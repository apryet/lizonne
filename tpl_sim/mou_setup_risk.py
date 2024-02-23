import os,shutil
import math 
import re
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
from pymarthe import MartheModel
from pymarthe.moptim import MartheOptim
from pymarthe.mfield import MartheField, MartheFieldSeries
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
import matplotlib.pyplot as plt
import pyemu
import helpers
encoding = 'latin-1'

# -------------------------------------------
# ------- Update model and post-proc -----------
# -------------------------------------------

print('Disabling initial steady state in .mart file ...')
martfile = 'Lizonne.mart'
with open(martfile, 'r', encoding=encoding) as f:
    mart_content = f.read()
    pattern = r'(.*)(\=Nombre maximal d\'itérations pour le pas de temps n.*)'
    repl = r'         0\2'
    mart_content = re.sub(pattern, repl, mart_content)

with open(martfile, 'w', encoding=encoding) as f:
    f.write(mart_content)

print('Update pastp file ...')
helpers.opt_pastp()

print('Reading MARTHE model ...')
# Spatial index has to be recomputed, reloading from file leads to error 
mm = MartheModel('Lizonne.rma', spatial_index = True)

# set time windows 
sim_start,sim_end = mm.get_time_window()

# load aquifer pumping data from MartheModel 
mm.load_prop('aqpump')
aqpump_df = mm.prop['aqpump'].data
aqpump_nodes = aqpump_df.node.unique()

# load river pumping data from MartheModel 
mm.load_prop('rivpump')
rivpump_df = mm.prop['rivpump'].data
rivpump_nodes = rivpump_df.node.unique()

# add reach number (aff_r) to rivpump
aff_r = MartheField('aff_r',mm.mlname+'.aff_r',mm, use_imask=False)
rivpump_df['aff_r'] = aff_r.data['value'][rivpump_df.node.values].astype(int)

# -------------------------------------------
# -- Compute maximum allowable pumping   ----
# -------------------------------------------
# => later used as constraint values and parameter bounds

# load steady state simulation 
mfs = MartheFieldSeries(mm=mm, chasim='chasim.out.cal')
mfs.load_field('CHARGE')

# get steady state simulated heads (option #1)
h0 = mfs.data['CHARGE'][0].data['value']

# get min of transient heads (option #2)
# => cells will be considered as confined if they are always confined
# => this leads to less restrictive constraints
isteps = list(mfs.data['CHARGE'].keys())[2:]
hs = np.stack([mfs.data['CHARGE'][i].data['value'] for i in isteps])
href=hs.min(axis=0)

# load geometry as 1D array
#top, botm = _top.flatten('F'), _hsubs.flatten('F')
top, botm = map(np.ravel,mm._get_top_bottom_arrays())

# identify confined cells 
conf = href > top

# compute aquifer thickness 
e = (top - botm)  # top - botm for confined cells
e[~conf] = (href  - botm)[~conf] # head - botm for unconfined cells

# compute minimum allowable head (used as optimization constraint)
hmin = top.copy() # confined 
hmin[~conf] = botm[~conf] + 1/3*e[~conf] # unconfined

# infer maximum allowable drawdown
smax = (href-top) # confined 
smax[~conf] = 2/3*e[~conf]  # unconfined

# infer maximum allowable pumping with Dupuit-Thiem (may be used as parameter bound)
permh = mm.prop['permh'].data['value'].flatten()
trans = permh*e
qmax = 2*np.pi*trans*smax/np.log(3) # m3/s

# aquifer pumping limits in df 
aqpump_lim_df = pd.DataFrame({
    'node':aqpump_nodes,
    'top': top[aqpump_nodes],
    'botm':botm[aqpump_nodes],
    'e':e[aqpump_nodes],
    'h0': h0[aqpump_nodes],
    'href': href[aqpump_nodes],
    'hmin': hmin[aqpump_nodes],
    'conf': conf[aqpump_nodes],
    'qmax': qmax[aqpump_nodes],
    }, index = [f'aqpump_{n:05d}' for n in aqpump_nodes])

aqpump_lim_df.to_csv('aqpump_lim.csv')

# --------------------------------------------
# -- append pumping wells to histo file
# --------------------------------------------
mm.build_modelgrid()
org_histo_file = open('Lizonne.histo.cal','r',encoding='ISO-8859-1')
new_histo_file = open('Lizonne.histo','w',encoding='ISO-8859-1')
ll = org_histo_file.readlines()
o = new_histo_file.writelines(ll[:-1])

for n in aqpump_nodes:
    xcc = mm.modelgrid.loc[n,'xcc']
    ycc = mm.modelgrid.loc[n,'ycc']
    l = mm.modelgrid.loc[n,'layer']
    # note l+1 : back to 1-based for layer id !
    o= new_histo_file.write(f'  /Charge       /HISTO/   =   /XCOO:X={xcc:>7.0f}Y={ycc:>7.0f}P={l+1:>7}; Name=aqpump_{n:05d}\n')


o = new_histo_file.writelines(ll[-1])

org_histo_file.close()
new_histo_file.close()

# -------------------------------------------
# ----- Run model with new set up -----------
# -------------------------------------------

print('Running Marthe model to update histo file...')
# run model once (so as to update dates_out in historiq.prn)
mm.run_model()

# write simulated values in syn_obs dir
syn_obs_dir = os.path.join('pest','syn_obs')
if os.path.exists(syn_obs_dir):
   shutil.rmtree(syn_obs_dir)

os.mkdir(syn_obs_dir)

helpers.opt_postproc(mm,syn_obs_dir)

# -------------------------------------------
# ---------- MartheOptim set up -------------
# -------------------------------------------

print('Setting up MartheOptim...')

# optimizatio over last simulation year
opt_start = pd.to_datetime(str(sim_end.year))

# typical PEST directory structure
dirs = {f'{k}_dir' : os.path.join('pest',k)
        for k in ['par', 'tpl', 'obs', 'ins', 'sim']}

# synthetic observations will be generated here 
dirs['obs_dir'] = syn_obs_dir

# clean and (re)-build directory tree
for k in ['par','tpl','ins','sim']:
    d = dirs[f'{k}_dir']
    if os.path.exists(d):
       shutil.rmtree(d)
    os.mkdir(d)

# initialize MartheOptim
mopt = MartheOptim(mm , name='opt_lizonne', **dirs)

# -------------------------------------------
# ------- Observations processing -----------
# -------------------------------------------

'''
3 types of "observations" are considered : 
    - Total pumping rate (objective #1, to maximize)
    - Total deficit at gaging  station (objective #2, to minimize)
    - Hydraulic heads at pumping wells (constraints)

Reference constraint values are not field observations but simulated outcomes. 
For convenience, the corresponding files are generated and stored 
in syn_obs dir to use mopt.add_obs() function. 

'''

# --- Total pumping (objective to maximize)

# ---  Pumping

# total pumping (objective #1, to maximize)
mopt.add_obs(data= pest_utils.read_simfile(os.path.join(mopt.obs_dir,'tot_pump.dat')), 
             obsnmes=['tot_pump'],
             locnme='tot_pump',
             weight= 0,
             datatype = 'volume', check_loc=False)

# total aquifer pumping (for informative purpose)
mopt.add_obs(data= pest_utils.read_simfile(os.path.join(mopt.obs_dir,'tot_aqpump.dat')),
             obsnmes=['tot_aqpump'],
             locnme='tot_aqpump',
             weight= 0,
             datatype = 'volume', check_loc=False)

# total river pumping (for informative purpose)
mopt.add_obs(data= pest_utils.read_simfile(os.path.join(mopt.obs_dir,'tot_rivpump.dat')), 
             obsnmes=['tot_rivpump'],
             locnme='tot_rivpump',
             weight= 0,
             datatype = 'volume', check_loc=False)


# ---  River deficits

# total deficit (objective #2, to minimize)
mopt.add_obs(data= pest_utils.read_simfile(os.path.join(mopt.obs_dir,f'deficit_tot.dat')), 
             obsnmes= ['deficit_tot'], 
             locnme = 'deficit_tot',
             obgnme = 'l_riv_deficit_m3', # l_ for "lower than" constraint
             weight= 0,
             datatype = 'volume', check_loc=False)

# deficit at stations, for informative purpose
station_locs = ['P8284010', 'P8215010', 'P7250001', 'P7270001']
for loc in station_locs:
    # add observation 
    mopt.add_obs(data= pest_utils.read_simfile(os.path.join(mopt.obs_dir,f'deficit_{loc}.dat')), 
                 obsnmes= [f'deficit_{loc}'], 
                 locnme = f'deficit_{loc}',
                 obgnme = 'l_riv_deficit_m3', # l_ for "lower than" constraint
                 weight= 1,
                 datatype = 'volume', check_loc=False)

# --- Observed levels at pumping wells  (constraints)

prn_df = marthe_utils.read_prn()

for loc in aqpump_lim_df.index:
    # pumping cell node
    aqpump_node = loc.split('_')[1]
    # write to syn_obs dir
    pest_utils.extract_prn(prn_df,loc, dates_out=prn_df.index, sim_dir=mopt.obs_dir)
    # add observation 
    obs_data = pest_utils.read_simfile(os.path.join(mopt.obs_dir,f'{loc}.dat'))
    mopt.add_obs(data= obs_data, 
                 locnme=loc, weight= 1,
                 # provide more explicit observation names 
                 obsnmes = [f'h_{aqpump_node}_n{i+1:03d}' for i in range(obs_data.shape[0])],
                 obgnme = f'g_pcellh_{aqpump_node}', # g_ for "greater than" constraint
                 datatype = 'head', check_loc=True)

# -------------------------------------------
# ---- Adjustable parameters (for UQ) -------
# -------------------------------------------

exec(open(os.path.join('..','master_ies','ies_parameterization.py')).read())

# --------- write template and parameter files of adjustable model parameters are handled by MartheOptim (dv are not).
mopt.write_parfile()
mopt.write_tplfile()

# ----------------------------------------------------------------
# --------- Climate model selection parameter --------------------
# ----------------------------------------------------------------
'''
# ----------------SKIPPED ---------------------------------: 
for the moment, climate model scenario is not considered as a parameter.
Instead, we keep up 3 quantiles in a population of simulation set account
for both inter-annual and climate scenario variability. 
but read from the directory name. 

cmid_df = pd.DataFrame({'parnme':['cm'],'defaultvalue':1})

# tpl file
cmid_tpl_file = os.path.join(dirs['tpl_dir'],'cm.tpl')
pest_utils.write_mlp_tplfile(cmid_tpl_file,cmid_df)

# par file 
cmid_par_file = os.path.join(dirs['par_dir'],'cm.dat')
pest_utils.write_mlp_parfile(cmid_par_file,cmid_df)
'''

# ----------------------------------------------------------------
# ---- Pumping parameterization (decision variables) -------------
# ----------------------------------------------------------------

'''

At the time of writting, PyMarthe does not handle factor-type parameters.
So we implement to level of parameterization, with overlay factor-type parameters for PEST
Pymarthe group pumping by cell, layer and istep. 
Factors (overlay parameters) are seen by PEST, grouped by :
    layer and istep for aquifer pumping
    reach id (aff_r) and istep for river pumping

writting of parameter values are handled by PyMarthe. The multiplication by factor is handled
in a preproc function from helpers, called in the forward run script. 

'''

# --- Aquifer pumping -------------------

# --- PyMarthe parameters (grouped by cell, layer and istep) 

# subset to summer months within opt period
ap_idx = pd.MultiIndex.from_frame(aqpump_df[['boundname', 'layer', 'istep']])
ap_dates = mm.mldates[ap_idx.get_level_values(2)]
# only summer values 
ap_kmi = ap_idx[ap_dates.month.isin(range(6,10)) & (ap_dates >= pd.to_datetime(opt_start))]

mopt.add_param(parname = 'aqpump',
               mobj = mm.prop['aqpump'],
               kmi = ap_kmi,
               value_col = 'value')


# write original parameter value file 
pest_utils.write_mlp_parfile(os.path.join(dirs['par_dir'],'aqpump.dat.org'), mopt.param['aqpump'].param_df)

# setup PEST parameters (factors), grouped by layer and istep
fac_parnmes = ap_kmi.to_frame().apply(lambda x: f'aqpumpfac__{x.layer}__{x.istep}',axis=1).unique()
aqpump_fac_df = pd.DataFrame({'parnme':fac_parnmes,'defaultvalue':1})

# tpl file
aqpump_tpl_file = os.path.join(dirs['tpl_dir'],'aqpumpfac.tpl')
pest_utils.write_mlp_tplfile(aqpump_tpl_file,aqpump_fac_df)

# par file 
aqpump_par_file = os.path.join(dirs['par_dir'],'aqpumpfac.dat')
pest_utils.write_mlp_parfile(aqpump_par_file,aqpump_fac_df)

# --- River pumping ----------------------

# --- PyMarthe parameters (grouped by cell, layer and istep)

# generate a multi-index from pumping df 
rp_idx = pd.MultiIndex.from_frame(rivpump_df[['boundname', 'aff_r', 'istep']])
rp_dates = mm.mldates[rp_idx.get_level_values(2)]

# subset to summer months within opt period 
rp_kmi = rp_idx[rp_dates.month.isin(range(6,10)) & (rp_dates >= pd.to_datetime(opt_start))]

mopt.add_param(parname = 'rivpump',
               mobj = mm.prop['rivpump'],
               kmi = rp_kmi,
               value_col = 'value')

# write original parameter value file 
pest_utils.write_mlp_parfile(os.path.join(dirs['par_dir'],'rivpump.dat.org'), mopt.param['rivpump'].param_df)

# --- setup PEST parameters (factors), grouped by reach and istep
rp_kmi_df = rp_kmi.to_frame()
rp_kmi_df['node'] = rp_kmi_df['boundname'].apply(lambda x: int(x.split('_')[1]))
rp_kmi_df['aff_r'] = aff_r.data['value'][rp_kmi_df.node.values].astype(int)

fac_parnmes = rp_kmi_df.apply(lambda x: f'rivpumpfac__{x.aff_r}__{x.istep}',axis=1).unique()
rivpump_fac_df = pd.DataFrame({'parnme':fac_parnmes,'defaultvalue':1})
rivpump_tpl_file = os.path.join(dirs['tpl_dir'],'rivpumpfac.tpl')
pest_utils.write_mlp_tplfile(rivpump_tpl_file, rivpump_fac_df)
rivpump_par_file = os.path.join(dirs['par_dir'],'rivpumpfac.dat')
pest_utils.write_mlp_parfile(rivpump_par_file,rivpump_fac_df)



# -----------------------------------------------------------------
# ------- if risk as an objective, introduce _risk_par  ------------
# -------------------------------------------------------------
risk_init = 0.50
risk_df = pd.DataFrame({'parnme':['_risk_'],'defaultvalue':risk_init})

# tpl file
risk_tpl_file = os.path.join(dirs['tpl_dir'],'risk.tpl')
pest_utils.write_mlp_tplfile(risk_tpl_file,risk_df)

# par file 
risk_par_file = os.path.join(dirs['par_dir'],'risk.dat')
pest_utils.write_mlp_parfile(risk_par_file,risk_df)

# -------------------------------------------
# ------- Generate PEST files  -----------
# -------------------------------------------

# -- Write PEST files
print('Writing PEST files ...')
mopt.write_insfile()

# copy simulated files from syn_obs 
for f in os.listdir(mopt.obs_dir):
    if f.endswith('.dat'):
        shutil.copy(os.path.join(mopt.obs_dir,f),os.path.join(mopt.sim_dir,f))

# -- Write configuration file
print('Writing configuration file ...')
mopt.write_config()

'''
# write forward run file (obsolete, further edits are necessary)
mopt.write_forward_run('forward_run.py.new',
                       f'{mopt.name}.config',
                       extra_py_imports=['helpers'],
                       preproc_functions = [helpers.lizonne_preproc],
                       extra_functions=[helpers.lizonne_postproc]
                       )
'''

# -- Build pst
print('Generating PEST Control File ...')

# collect PEST io files
tpl_files = mopt.collect_pest_files('tpl')+[os.path.join(dirs['tpl_dir'],'risk.tpl')]
par_files = mopt.collect_pest_files('par')+[os.path.join(dirs['par_dir'],'risk.dat')]
ins_files = mopt.collect_pest_files('ins')
out_files = mopt.collect_pest_files('out')

# generate pst from io files 
pst = pyemu.Pst.from_io_files(tpl_files, par_files, ins_files, out_files)

pst.model_command = 'python forward_run.py'
obs_df = mopt.get_obs_df(transformed=True)
pst.observation_data.loc[obs_df.index.str.lower()] = obs_df[pst.obs_fieldnames]

# -------------------------------------------------------
# ------- Set parameter  groups   ------------
# -------------------------------------------------------

par = pst.parameter_data

# decision variables
par.loc[par.parnme.str.startswith('aq'), 'pargp'] = 'aqpump'
par.loc[par.parnme.str.startswith('riv'),'pargp'] = 'rivpump'

# --- set parameter groups and update pst
# soil : parameter name as group name 
idx = par.pargp=='soil'
par.loc[idx,'pargp'] = par.loc[idx,'parnme'].apply(lambda x : '_'.join(x.split('_')[:-1]))

# emmca : parameter name + layer id as group name 
idx = par.parnme.str.startswith('emmca')
par.loc[idx,'pargp'] = par.loc[idx,'parnme'].apply(lambda x : '_'.join(x.split('_')[:2]))

# permh : parameter name + layer id as group name 
idx = par.parnme.str.startswith('permh')
par.loc[idx,'pargp'] = par.loc[idx,'parnme'].apply(lambda x : '_'.join(x.split('_')[:2]))

# emmli : parameter name + layer id as group name 
idx = par.parnme.str.startswith('emmca')
par.loc[idx,'pargp'] = par.loc[idx,'parnme'].apply(lambda x : '_'.join(x.split('_')[:2]))

# risk as an objective
par.loc['_risk_', 'pargp'] = 'risk'
par.loc['_risk_', 'partrans'] = 'none'
par.loc['_risk_', 'parval1'] = risk_init
par.loc['_risk_', 'parubnd'] = 0.99
par.loc['_risk_', 'parlbnd'] = 0.50

# update pgroups
pst.rectify_pgroups() 

# --------------------------------------------------------
# ---- get parameter ensemble (stack from IES)  ---------
# -------------------------------------------------------
'''
'''
# load IES final filtered parameter ensemble
ffpt_pe = pyemu.ParameterEnsemble.from_csv(pst=pst,filename='cal_lizonne.ffpt.par.csv')

'''
# NOTE : SECTION SKIPPED, cm_id no more considered as an adjustable parameter 

# ---- Add climate model selection parameter to the stack 
# load list of climate models 
cm_df = pd.read_excel('../data/DRIAS_Lizonne/clim_models.xlsx',index_col=0)
num_cm = cm_df.shape[0]

# append cm id to the index as a suffix ( real_id + clim_id)
pstack.index= pstack.index.get_level_values(1) + '_' +  pstack.index.get_level_values(0).astype(str)

# append columns referring to cm id 
pstack['cm'] = [ i for i in cm_df.index for j in range(ffpt_pe.shape[0])]
'''

# Append decision variables to the parameter ensemble 
dvpars = par.index[par.index.str.contains('fac') | par.index.str.contains('risk') ]
dv_pe = pd.DataFrame(1,index=ffpt_pe.index,columns=dvpars)
mou_pe = pyemu.ParameterEnsemble(pst, pd.concat((ffpt_pe._df,dv_pe),axis=1))
mou_pe.to_csv('cal_lizonne.par_stack.csv')

npar_reals = mou_pe.shape[0]
print(f'Size of the parameter stack : {npar_reals}')

# -------------------------------------------------------
# ------- Set parameter  values and bounds   ------------
# -------------------------------------------------------

# Set parameter value from base realization and set as fixed
par = pst.parameter_data

# Set parameter value from base realization and set as fixed
par.loc[mou_pe.columns,'parval1'] = mou_pe.loc['base'].values

# -- adjusting dv parameter settings
par.loc[dvpars,'parlbnd']=0
par.loc[dvpars,'parubnd']=10
par.loc[dvpars,'parval1']=0.
par.loc[dvpars,'partrans']='none'

# -- adjustable model  parameter bounds 
ies_pst = pyemu.Pst(os.path.join('..','master_ies',f'cal_lizonne.pst'))
ies_par = ies_pst.parameter_data 
par.loc[ies_par.index,'parlbnd'] = ies_par.parlbnd
par.loc[ies_par.index,'parubnd'] = ies_par.parubnd

# ---------------------------------------------
# ------- Set observation values   ------------
# ---------------------------------------------

obs = pst.observation_data
hmin = aqpump_lim_df.set_index('node')['hmin']

# set minimum allowable head values in pumping wells 
for n in hmin.index:
    idx = obs.obgnme.str.endswith(str(n))
    obs.loc[idx,'obsval'] = hmin.loc[n] 

# set 0-weight to head observations before opt_date
hobs = obs.loc[obs.index.str.startswith('h_')].copy()
hobs['rec_id'] = hobs.obsnme.str.extract(pat=r'n(\d+)').astype(int)-1
hobs['date'] = hobs.rec_id.apply(lambda x : prn_df.index[x])
obs.loc[hobs.loc[hobs.date<opt_start].index,'weight']=0

# set 0-weight to head of pumping cell offending constraints with 0 pumping.
#discard_locs = ['08656','18427','18428', '19029','19109','19505','20198','20514','20609','20976'] 

discard_locs = ['08656']
for loc in discard_locs:
    obs.loc[obs.obgnme=='g_pcellh_'+loc,'weight']=0

# ---------------------------------------------
# --- draw dv population   --------------------
# ---------------------------------------------

# NOTE : unsure whether pestpp considers values of calibrated parameter values 
# so it may not be necessary to provide them in the dv pop file. 

# dv population size (twice the number of dv)
num_dvpars = 2*dvpars.shape[0]
print(f'Number of dv : {num_dvpars}')
num_reals=2*116 # adjusted as a multiple of available cores (116)
print(f'Size of dv population : {num_reals}')

# set all parameters fixed and release dvpars 
par.loc[:,'partrans'] = 'fixed'
par.loc[dvpars,'partrans'] = 'none'

# lower the factor for aquifer pumping to avoid unsatisfied constraints (infeasible solutions)
par.loc[dvpars,'parubnd'] = 0.50
par.loc['_risk_', 'parubnd'] = 0.99

# generate initial dv ensemble by uniform draw and save to csv 
dvpop = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=num_reals)
dvpop.to_csv('initial_dvpop.csv')
par.loc[par.partrans=='fixed','partrans'] = 'none'

# tell PESTPP-MOU about the new file
pst.pestpp_options["mou_dv_population_file"] = 'initial_dvpop.csv'

# reset the decision variable upper bound
par.loc[dvpars,'parubnd'] = 10.
par.loc['_risk_', 'parubnd'] = 0.99

# ---------------------------------------------
# --- define dv and obs   --------------------
# ---------------------------------------------

#  decision variables
pst.pestpp_options['opt_dec_var_groups'] = ['aqpump','rivpump','risk']

# discard constraints #1 on deficit, they are now an  objective to minimize
obs.loc[obs.obsnme.str.startswith('deficit'),'weight']=0

# objectives 
#pst.pestpp_options["mou_objectives"] = ['deficit_tot','tot_pump']
pst.pestpp_options["mou_objectives"] = ['deficit_tot','tot_pump','_risk_']

# now define the direction of each objective and give them a non-zero weight:
obs = pst.observation_data
obs.loc['deficit_tot','obgnme'] = 'less_than_obj'
obs.loc['tot_pump', 'obgnme'] = 'greater_than_obj' 
obs.loc[['deficit_tot','tot_pump'],'weight'] = 1.0


# ---------------------------------------------
# -------- MOEA settings   --------------------
# ---------------------------------------------

pst.pestpp_options["mou_population_size"] = num_reals 
pst.pestpp_options["mou_save_population_every"] = 1

# risk as an objective
pst.pestpp_options['mou_risk_objective'] = True
pst.add_pi_equation(['_risk_'], pilbl = '_risk_', obs_group = 'greater_than')

pst.pestpp_options['opt_risk'] = 0.66
pst.pestpp_options['opt_par_stack'] = 'cal_lizonne.par_stack.csv'
#pst.pestpp_options['opt_obs_stack'] = 'cal_lizonne.obs_stack.csv'
pst.pestpp_options['opt_chance_points'] = 'single'

# NOTE : issue in pestpp-manual
pst.pestpp_options['opt_stack_size'] = npar_reals
pst.pestpp_options['opt_recalc_chance_every'] = 100

pst.pestpp_options['mou_generator'] = 'pso'
pst.pestpp_options['mou_env_selector'] = 'nsga'
pst.pestpp_options['mou_verbose_level'] = 4

# ---------------------------------------------
# -------- Write control files   --------------
# ---------------------------------------------

# ---- Write pst
pst.control_data.noptmax = 50
pst_name = 'mou_lizonne.pst'
pst.write(pst_name)

# --- write fac0 pst
par.loc[dvpars,'parval1']=0.
pst.control_data.noptmax = 0
pst_name = 'mou_lizonne_fac0.pst'
pst.write(pst_name)

# --- write fac1 pst
par.loc[dvpars,'parval1']=1.
pst.control_data.noptmax = 0
pst_name = 'mou_lizonne_fac1.pst'
pst.write(pst_name)

# ---------------------------------------------
# -------- run PEST with reference conditions -
# ---------------------------------------------
pyemu.os_utils.run('pestpp-mou mou_lizonne_fac0.pst')
pyemu.os_utils.run('pestpp-mou mou_lizonne_fac1.pst')

