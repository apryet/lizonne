import os, sys, shutil
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
from pymarthe import MartheModel
from pymarthe.moptim import MartheOptim
from pymarthe.mfield import MartheField
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
import matplotlib.pyplot as plt
import pyemu

# --- load model
print('Reading MARTHE model ...')
# Spatial index has to be recomputed, reloading from file leads to error 
mm = MartheModel('Lizonne.rma', spatial_index = True)

marthe_utils.set_tw(start=0,mm=mm)

# -------------------------------------------
# ---------- MartheOptim set up -------------
# -------------------------------------------

print('Setting up MartheOptim...')
dirs = {f'{k}_dir' : os.path.join('pest',k)
        for k in ['par', 'tpl', 'obs', 'ins', 'sim']}

# clean and (re)-build directory tree
for k in ['par','tpl','ins','sim']:
    d = dirs[f'{k}_dir']
    if os.path.exists(d):
       shutil.rmtree(d)
    os.mkdir(d)

# initialize MartheOptim
mopt = MartheOptim(mm , name='cal_lizonne', **dirs)

# -------------------------------------------
# ------- Observations processing -----------
# -------------------------------------------

# -- Add observations data
print('Adding observations ...')
obs_data_files = [f for f in os.listdir(mopt.obs_dir) if f.endswith('.dat')]

# weights 
sigma_head = 0.10 # m
sigma_rivflow = np.log10(1.1) # on log-transformed(Q), i.e. factor of 1.1
sigma_dic={'head':sigma_head,'rivflow':sigma_rivflow}

# transformation 
trans_dic={'head':'none','rivflow':np.log10}

for f in obs_data_files:
    # observation file name
    obsfile = os.path.join(mopt.obs_dir,f)
    locnme = os.path.split(obsfile)[-1].split('.')[0]
    # data type
    dt = 'rivflow' if locnme.startswith('P') else 'head'
    # -- add observations (with reduced time window)
    mopt.add_obs(data= obsfile, locnme=locnme, 
                 weight= 1./sigma_dic[dt],
                 datatype = dt, check_loc=True)

# -------------------------------------------
# ------- Parameter processing -----------
# -------------------------------------------

# ------ parameterization of soil properties ----

print('Parameterize soil properties ...')
mm.load_prop('soil')
ms = mm.prop['soil']
skmi = pest_utils.get_kmi(ms, keys=['soilprop', 'zone'], istep=0) # Constant in time
mopt.add_param(parname='soil', mobj=ms, kmi=skmi, parlbnd=0, parubnd=1e5)

# --- Parameterization of  perm_r ----
# -> zpc based on aff_r id

print('Parameterize `perm_r` field ...')
# set up izone for perm_r based on reach id  
mm.load_prop('aff_r',use_imask=False) # reach id 
aff_r = mm.prop['aff_r']
iperm_r = MartheField('iperm_r', 1, mm,use_imask=False)
iperm_r.data['value'] = -1 * aff_r.data['value'] # negative values for ZPCs

# initialize perm_r field from scratch 
perm_r = MartheField('perm_r', 0, mm,use_imask=False)
perm_r.data['value'][aff_r.data['value']!=0]=1e-5
mm.prop['perm_r'] = perm_r

# set up parameter
mopt.add_param(parname='perm_r', mobj=mm.prop['perm_r'],
               izone=iperm_r,
               trans='log10', btrans='lambda x: 10**x',
               parlbnd=1e-8, parubnd=1e-2, defaultvalue=1e-5,
               use_imask=False)


# --- Parameterization of  emmli  ---
# -- >zpc for all layers

emmli_value = 0.1 

print('Parameterize `emmli` field ...')
# initialize emmli field from scratch
mm.prop['emmli'] =  MartheField('emmli', emmli_value, mm)

mopt.add_param(parname='emmli', mobj=mm.prop['emmli'],
               izone=None, 
               trans='log10', btrans='lambda x: 10**x',
               parlbnd=1e-4, parubnd=1, defaultvalue=emmli_value)


# --- Set izone and pilot point spacings (permh and emmca)

# izone definition
izone = MartheField('izone', 1, mm, use_imask=True)

# set spacings depending on layer type 
spacing_ep_layer = 8*1000 # in meters
spacing_aq_layer = 4*1000 # in meters
ep_mask = mm.layers_infos.name.str.startswith('Eponte')
spacings = ep_mask.map({ True : spacing_ep_layer, False : spacing_aq_layer})

# --- Parameterization of  permh  ----
print('Parameterize `permh` field ...')

# pilot points set up
pp_permh = pp_utils.PilotPoints(izone)
for l,s in enumerate(spacings):
    pp_permh.add_spacing_pp(layer=l, zone=1, xspacing=s, yspacing=s)

# add parameters 
mopt.add_param(parname='permh', mobj=mm.prop['permh'],
               izone=izone, pp_data=pp_permh.to_pp_data(),
               trans='log10', btrans='lambda x: 10**x',
               parlbnd=1e-12, parubnd=1) 

# compute kriging factors
mopt.write_kriging_factors(pp_permh.extract_vgm_range(), parname='permh', cleanup=True, save_cov=True)

# --- Parameterization of  emmca  ---
print('Parameterize `emmca` field ...')
mm.load_prop('emmca')

# pilot points set up
pp_emmca = pp_utils.PilotPoints(izone)
for l,s in enumerate(spacings):
    pp_emmca.add_spacing_pp(layer=l, zone=1, xspacing=s, yspacing=s)

# add parameters 
mopt.add_param(parname='emmca', mobj=mm.prop['emmca'],
               izone=izone, pp_data=pp_emmca.to_pp_data(),
               trans='log10', btrans='lambda x: 10**x',
               parlbnd=1e-15, parubnd=1) 

# compute kriging factors
mopt.write_kriging_factors(pp_emmca.extract_vgm_range(), parname='emmca', cleanup=False, save_cov=True)

# -------------------------------------------
# ------- Generate PEST files  -----------
# -------------------------------------------

# --- Write new model property files  ----

# re-write properties (since we generated emmli, emmca, and perm_r from scratch) 
# values will be updated by set_data_from_parfile() with first PEST run, but masked
# values (in dmv list) will no longer be altered, this is the last opportunity !
mm.write_prop()

# -- Write PEST files
print('Writing PEST files ...')
mopt.write_insfile()
mopt.write_simfile()
mopt.write_parfile()
mopt.write_tplfile()

# -- Write configuration file
print('Writing configuration file ...')
mopt.write_config()

# -- Build pst
print('Writing PEST Control File ...')
pst = mopt.build_pst(add_reg0= False,
                     write_pst= False,
                     write_config= True,
                     write_fr= True,
                     noptmax=0)

# temporary fix (will be implemented in build_pst()
pst.rectify_pgroups()

# --------------------------------------------------------------
# ------- Set parameter initial values and bounds   ------------
# --------------------------------------------------------------

par = pst.parameter_data

# -- adjusting parameter bounds
# lower and upper bounds representative sigma range of 4

# soil parameters
par.loc[par.parnme.str.startswith('cap_sol_progr'),'parlbnd']=0       # mm
par.loc[par.parnme.str.startswith('cap_sol_progr'),'parval1']=269.55  # mm
par.loc[par.parnme.str.startswith('cap_sol_progr'),'parubnd']=1000    # mm

par.loc[par.parnme.str.startswith('equ_ruis_perc'),'parlbnd']=0        # mm
par.loc[par.parnme.str.startswith('equ_ruis_perc'),'parval1']=14.26    # mm
par.loc[par.parnme.str.startswith('equ_ruis_perc'),'parubnd']=1000     # mm

par.loc[par.parnme.str.startswith('t_demi_percol'),'parlbnd']=0      # days
par.loc[par.parnme.str.startswith('t_demi_percol'),'parval1']=0.767  # days
par.loc[par.parnme.str.startswith('t_demi_percol'),'parubnd']=100    # days

# unconfined storage (emmli)
par.loc[par.parnme.str.startswith('emmli'),'parlbnd']=-4   # [-]
par.loc[par.parnme.str.startswith('emmli'),'parubnd']=0      # [-]

# confined storage (emmca, log-transformed)
par.loc[par.parnme.str.startswith('emmca'),'parlbnd']=-12      # [-]
par.loc[par.parnme.str.startswith('emmca'),'parubnd']=-3      # [-]

# river network conductivity (uniform)
par.loc[par.parnme.str.startswith('perm_r'),'parlbnd'] = -10    # m/s
par.loc[par.parnme.str.startswith('perm_r'),'parubnd'] = -2     # m/s

# hydraulic conductivity - EPCS
par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l01') ,'parlbnd']=-12  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l01') ,'parval1']=-8  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l01') ,'parubnd']=-2  # log10(m/s)

# hydraulic conductivity - COST
par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l02') ,'parlbnd']=-9  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l02') ,'parval1']=-5  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l02') ,'parubnd']=-1  # log10(m/s)

# hydraulic conductivity - EPTU
par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l03') ,'parlbnd']=-12  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l03') ,'parval1']=-8  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l03') ,'parubnd']=-5  # log10(m/s)

# hydraulic conductivity - TURO
par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l04') ,'parlbnd']=-9  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l04') ,'parval1']=-5  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l04') ,'parubnd']=-1  # log10(m/s)

# hydraulic conductivity - EPCE
par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l05') ,'parlbnd']=-12  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l05') ,'parval1']=-8  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l05') ,'parubnd']=-5  # log10(m/s)

# hydraulic conductivity - CE
par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l06') ,'parlbnd']=-9  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l06') ,'parval1']=-5  # log10(m/s)

par.loc[par.parnme.str.startswith('permh') & \
        par.parnme.str.contains('l06') ,'parubnd']=-1  # log10(m/s)

# -----------------------------------------------------------------
# generate parameter covariance matrix 
# -----------------------------------------------------------------

print('Build prior covariance matrix ...')
# -- Build full covariance matrix from parameter bounds
par_cov = pyemu.Cov.from_parameter_data(pst, sigma_range=4)

# -- Load and replace off diagonal terms from pilot points geostructures
for pp_cov_file in os.listdir(mopt.par_dir):
    if pp_cov_file.endswith('.jcb'): # par_dir should be clear of parasitic .jcb
        fn = os.path.join(mopt.par_dir, pp_cov_file)
        # load correlation matrix of pilot points
        corr = pyemu.Cov.from_binary(fn)
        # convert correlation matrix to covariance matrix
        sigma_vec = par_cov.get(corr.names,corr.names).get_diagonal_vector().sqrt
        cov = corr.hadamard_product(pyemu.Matrix(x=np.outer(sigma_vec.x,sigma_vec.x),
                   row_names=sigma_vec.names,
                   col_names=sigma_vec.names))
        par_cov.replace(cov)

# -- Write final prior parameter covariance matrix
par_cov.to_binary('prior.cov')

'''
# --  generate parmeter ensemble for IES
num_reals = 496
prior_pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov=par_cov, num_reals=num_reals)
prior_pe.add_base()
prior_pe.enforce()
prior_pe.to_binary('prior_pe.jcb')

# plot selection of parameters
#prior_pe._df.loc[:,].plot(kind="hist",bins=20,alpha=0.5)

'''

# -----------------------------------------------------------------
# generate observation covariance matrix 
# -----------------------------------------------------------------

# get observation dataframe with full data (observation dates and types) 
obs_df = mopt.get_obs_df()

# take a peek to observation auto-correlation plots
'''
for locnme in obs_df.locnme.unique(): 
    try :
        # NOTE check lag unit, not sure wheter it remains relevant for irregular dates
        ts = pd.Series(obs_df.loc[obs_df.locnme==locnme,'obsval'].values,
                       index=obs_df.loc[obs_df.locnme==locnme,'date'].values)
        fig,ax = plt.subplots(figsize=(5,4))
        pd.plotting.autocorrelation_plot(ts,ax=ax)
        fig.savefig(os.path.join('figs','obs_autocorr',f'{locnme}.png'),dpi=300)
    except :
        continue
'''

''''
# => after a rough anlaysis, we get auto-corr. of ~100 days for both Q and heads
v = pyemu.geostats.ExpVario(a=33,contribution=1.0)

# get full (diagonal) observation covariance matrix from obs weights 
obs_cov = pyemu.Cov.from_observation_data(pst)

# -- Load and replace off diagonal terms from autocorrelation
print('Computing off-diag terms of observation covariance matrix...')
for locnme in obs_df.locnme.unique():
    # get observation names for current loc
    onmes = obs_df.loc[obs_df.locnme==locnme].index.to_list()
    # get correlation matrix from variogram model
    dates=pd.to_datetime(obs_df.loc[obs_df.locnme==locnme,'date'])
    x = (dates - dates[0]).astype('timedelta64[D]').values  # days from first obs.
    y = np.zeros_like(x)  # dummy constant vector for Vario2d methods
    corr = v.covariance_matrix(x,y,names=onmes)
    # get covariance matrix from correlation matrix scaled by the sigma_i*sigma_j matrix
    # we get this matrix with the outer product of the sigma_vec by itself.
    # hadamard product corresponds to np.multiply (matrix term-by-term product)
    sigma_vec = obs_cov.get(onmes,onmes).get_diagonal_vector().sqrt
    cov = corr.hadamard_product(
            pyemu.Matrix(x=np.outer(sigma_vec.x,sigma_vec.x),
                         row_names=onmes,
                         col_names=onmes)
            )
    obs_cov.replace(cov)

obs_cov.to_binary('obs_cov.jcb')

# generate observation ensemble for IES
oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst=pst, 
                                                num_reals=num_reals,
                                                cov=obs_cov) 

oe.add_base()
oe.to_binary('oe.jcb')
'''

# -----------------------------------------------------------------
# -- Set up some weighting and regularization adjustments 
# -----------------------------------------------------------------

# load desired contribution to phi for each obs groups (rough and subjective)
obsgp_contrib = pd.read_excel('obsgp_contrib.xlsx')
balanced_groups = {g:w for g,w in zip(obsgp_contrib.name,obsgp_contrib.contrib)}

# get res file from last calibration 
pst.set_res('../master_hp_2/caleval_lizonne.rei')

# adjust weights
pst.adjust_weights(obsgrp_dict=balanced_groups)

# target objective function value 
pst.reg_data.phimlim = 0.5*obsgp_contrib['contrib'].sum()
pst.reg_data.phimaccept = 1.05*pst.reg_data.phimlim
pst.reg_data.fracphim = 0.2

# zero-order Tikhonov reg 
pyemu.helpers.zero_order_tikhonov(pst)

# 1st-order Tikhonov reg 
pyemu.helpers.first_order_pearson_tikhonov(pst,par_cov,reset=False,abs_drop_tol=0.2)

# -----------------------------------------------------------------
# -- Set up some additional ++ options
# -----------------------------------------------------------------

print('Set some PEST ++ options ...')

# -- Set parameter prior covariance matrix
#pst.pestpp_options["parcov"] = "prior.cov"

# -- Derivative computation options 
pst.parameter_groups.loc["soil","inctyp"] = 'relative'
pst.parameter_groups.loc["soil","derinc"] = 0.25
pst.parameter_groups.loc["soil","derinclb"] = 0.01

# absolute increment type for log-transformed parameters
pst.parameter_groups.loc[['emmli_zpc','permh_pp','emmca_pp','perm_r_zpc'],'inctyp']='absolute'
pst.parameter_groups.loc[['emmli_zpc','permh_pp','emmca_pp','perm_r_zpc'],'derinc']=0.1 # factor of 25%

# Parameter upgrade limit options 
par.parchglim='absolute(1)' # for pest_hp only
par.loc[par.pargp=='emmli_zpc','parchglim'] = 'factor'
par.loc[par.pargp=='soil','parchglim'] = 'factor'
pst.control_data.facparmax = 5

# -- Set PEST PP GLM options 
pst.pestpp_options['svd_pack'] = 'redsvd'
pst.pestpp_options['uncertainty'] = 'false'

# slower but more robust (avoids 0-length parameter upgrade vector after bounds enforcement)


# -- Set overdue rescheduling factor to twice the average model run
pst.pestpp_options['overdue_resched_fac'] = 2
pst.pestpp_options['panther_agent_no_ping_timeout_secs'] = 36000

# -- IES options 
'''
pst.pestpp_options["ies_num_reals"] = num_reals
pst.pestpp_options['ies_parameter_ensemble'] = 'prior_pe.jcb'
pst.pestpp_options["ies_observation_ensemble"] = "oe.jcb"
pst.pestpp_options["ies_bad_phi_sigma"] = 2.0
'''


# -----------------------------------------------------------------
# -- Write final PEST control file - End of the story
# -----------------------------------------------------------------

# -- Write .pst file
pst.control_data.noptmax=15

pst.write(mopt.name + '.pst')

# ugly trick to append absparmax value
with open(mopt.name + '.pst') as f :
    ls = f.readlines()

ls[6] = ls[6].replace('\n','   absparmax(1)=2\n')

with open(mopt.name + '.pst','w') as f:
    f.writelines(ls)

