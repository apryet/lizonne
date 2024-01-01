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

# -------------------------------------------
# ---------- load  MartheModel  -------------
# -------------------------------------------

# --- load model
print('Reading MARTHE model ...')
# Spatial index has to be recomputed, reloading from file leads to error 
mm = MartheModel('Lizonne.rma', spatial_index = True)

# read model output (historiq.prn should be up to date !!!)
#mm.run_model() # run model if .pastp has been modified
prn_df = marthe_utils.read_prn('historiq.prn')

# -------------------------------------------
# ---------- MartheOptim set up -------------
# -------------------------------------------

# set optimization start
# sim starts on 2010-07-31, ends on 2019-07-31
# 18 months warm up
cal_start = pd.to_datetime('2012-07-31')

print('Setting up MartheOptim...')
# clean and (re)-build directory tree
dirs = {f'{k}_dir' : os.path.join('pest',k)
        for k in ['par', 'tpl', 'obs', 'ins', 'sim']}

for k in ['par','tpl','ins','sim']:
    d = dirs[f'{k}_dir']
    if os.path.exists(d):
       shutil.rmtree(d)
    os.mkdir(d)

# initialize MartheOptim
mopt = MartheOptim(mm, name='cal_lizonne', **dirs)

# -------------------------------------------
# ------- Observations processing -----------
# -------------------------------------------

# -- Add observations data
print('Adding observations ...')
obs_data_files = [f for f in os.listdir(mopt.obs_dir) if f.endswith('.dat')]

# standard deviations based on a decent model-to-measurement misfit at ±2*sigma
sigma_head = 1.00 # m
sigma_rivflow = np.log10(1.15) # on log-transformed(Q), i.e. factor of 1.15
sigma_dic={'head':sigma_head,'rivflow':sigma_rivflow}

# transformation 
qtrans = 'lambda  x: np.log10(max(x,1e-3))' # minimum flow value of 1 L/s
trans_dic={'head':'none','rivflow':qtrans}


for f in obs_data_files:
    # observation file name
    obsfile = os.path.join(mopt.obs_dir,f)
    locnme = os.path.split(obsfile)[-1].split('.')[0]
    # data type
    dt = 'rivflow' if locnme.startswith('P') else 'head'
    # read observations data
    odf = marthe_utils.read_obsfile(obsfile)
    # keep only records over calibration period (after warm up)
    odf = odf.loc[odf.index>cal_start]
    if odf.shape[0]>0 :
        mopt.add_obs(data= odf, locnme=locnme, 
                     weight= 1./sigma_dic[dt],
                     datatype = dt, trans=trans_dic[dt],
                     check_loc=True)
        # remove duplicated sim values
        obs_df = mopt.obs[locnme].obs_df
        sim_df = prn_df.loc[:,(slice(None),locnme)]
        obs_df['simval'] = sim_df.loc[obs_df.date].iloc[:,0].values
        obs_df.drop_duplicates('simval',keep='first',inplace=True)


# add fluctuations from median for head observations 
hobs_locs = [loc for loc in mopt.obs.keys() if not loc.startswith('P')]
mopt.add_fluc(locnme=hobs_locs,tag='m',on='median',weight=1./sigma_dic['head']) # obsgmne <12 char

# -------------------------------------------
# ------- Parameter processing -----------
# -------------------------------------------

exec(open('ies_parameterization.py').read())

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
pest_utils.run_from_config(f'{mopt.name}.config',run_model=False)

# -- Build pst
print('Writing PEST Control File ...')
pst = mopt.build_pst(add_reg0= False,
                     write_pst= False,
                     write_config= True,
                     write_fr= True,
                     noptmax=0)

# --------------------------------------------------------------
# ------- Set parameter initial values and bounds   ------------
# --------------------------------------------------------------

par = pst.parameter_data

# more granularity in parameter groups 

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

pst.rectify_pgroups()
                                                                               
# -- adjusting parameter bounds
# lower and upper bounds representative sigma range of 4

# soil parameters (see sec. 3.4.9 notice Gardenia)
par.loc[par.parnme.str.startswith('cap_sol_progr'),'parlbnd']=100       # mm
par.loc[par.parnme.str.startswith('cap_sol_progr'),'parval1']=250  # mm
par.loc[par.parnme.str.startswith('cap_sol_progr'),'parubnd']=400    # mm

par.loc[par.parnme.str.startswith('equ_ruis_perc'),'parlbnd']=100        # mm
par.loc[par.parnme.str.startswith('equ_ruis_perc'),'parval1']=350    # mm
par.loc[par.parnme.str.startswith('equ_ruis_perc'),'parubnd']=600     # mm

par.loc[par.parnme.str.startswith('t_demi_percol'),'parlbnd']= 1   # days
par.loc[par.parnme.str.startswith('t_demi_percol'),'parval1']= 15  # days
par.loc[par.parnme.str.startswith('t_demi_percol'),'parubnd']= 30  # days

# unconfined storage (emmli)
par.loc[par.parnme.str.startswith('emmli'),'parlbnd']=-4   # [-]
par.loc[par.parnme.str.startswith('emmli'),'parval1']=-2   # [-]
par.loc[par.parnme.str.startswith('emmli'),'parubnd']=0      # [-]

# confined storage (emmca, log-transformed)
par.loc[par.parnme.str.startswith('emmca'),'parlbnd']=-12      # [-]
par.loc[par.parnme.str.startswith('emmca'),'parval1']=-7      # [-]
par.loc[par.parnme.str.startswith('emmca'),'parubnd']=-2      # [-]

# river network conductivity (uniform)
par.loc[par.parnme.str.startswith('perm_r'),'parlbnd'] = -10    # m/s
par.loc[par.parnme.str.startswith('perm_r'),'parval1'] = -6.5    # m/s
par.loc[par.parnme.str.startswith('perm_r'),'parubnd'] = -3     # m/s

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
        par.parnme.str.contains('l05') ,'parubnd']=-3  # log10(m/s)

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
par_cov.to_coo('prior_pcov.jcb')

# --  generate parameter ensemble for IES
num_reals = 2*116 # multiple of the number of available CPU cores
prior_pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov=par_cov, num_reals=num_reals)
prior_pe.enforce()
prior_pe.to_binary('prior_pe.jcb')

# plot selection of parameters
#prior_pe._df.loc[:,].plot(kind="hist",bins=20,alpha=0.5)


# -----------------------------------------------------------------
# adjust parameter bounds after drawing  
# -----------------------------------------------------------------

# river network conductivity (uniform)
#par.loc[par.parnme.str.startswith('perm_r'),'parlbnd'] = -7    # m/s
#par.loc[par.parnme.str.startswith('perm_r'),'parubnd'] = 0    # m/s

# -----------------------------------------------------------------
# generate observation covariance matrix 
# -----------------------------------------------------------------

# get observation dataframe with full data (observation dates and types) 
obs_df = mopt.get_obs_df()

'''
# take a peek to observation auto-correlation plots
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

obs_cov.to_coo('obs_cov.jcb')

'''
## reload saved obs. cov mat
obs_cov = pyemu.Cov.from_binary('obs_cov.jcb')

# generate observation ensemble for IES
oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst=pst, 
                                                num_reals=num_reals,
                                                cov=obs_cov) 

oe.add_base()
oe.to_binary('oe.jcb')

# -----------------------------------------------------------------
# -- Adjust observation group name and set up weighting  adjustments 
# -----------------------------------------------------------------

obs = pst.observation_data

# reduce the weight of fluctuations to the third of absolute values
obs.loc[obs.obgnme.str.endswith('mf'),'weight'] = obs.loc[obs.obgnme.str.endswith('mf'),'weight'].div(2)

# load xls file with desired contribution of obs group in PHI ("contrib" column)
obsgp_data_file = 'obsgp_data.xlsx'
obsgp_data = pd.read_excel(obsgp_data_file,index_col='name')

# number of non-zero weighted observations per group (nnz)
obsgp_nnz = obs.loc[obs.weight>0].groupby('obgnme').count()['obsval']

# original obs loc (not fluctuations)
locs = obsgp_nnz.loc[~obsgp_nnz.index.str.endswith('mf')].index
# get nnz
obsgp_data.loc[locs,'nnz'] = obsgp_nnz.loc[locs].astype(int)

# compute phi factor from user-defined contribution provided in obsgp_data_file 
sum_contrib = obsgp_data.loc[obsgp_data.nnz>0,'contrib'].sum()
obsgp_data['phi_factor'] = obsgp_data.contrib/sum_contrib

# write back excel file  (we have added nnz and phi_factor columns)
obsgp_data.to_excel(obsgp_data_file)

# column of obsgp_data to consider as target phi contribution
phi_contrib = 'phi_factor'

# write ies_phi_factor file as csv 
# locs are obs group "tags" in phi factors file, so fluctuations will also be considered
obsgp_data.loc[locs,phi_contrib].to_csv('ies_phi_factor.csv', index=True, header=False)

# -----------------------------------------------------------------
# -- some more ++ options
# -----------------------------------------------------------------

print('Set some PEST ++ options ...')

# Parameter upgrade limit options 
par.parchglim='absolute(1)' # for pest_hp only
par.loc[par.pargp=='emmli_zpc','parchglim'] = 'factor'
par.loc[par.pargp=='soil','parchglim'] = 'factor'
pst.control_data.facparmax = 5

# -- Set PEST SVD options 
pst.pestpp_options['svd_pack'] = 'redsvd'
pst.svd_data.maxsing = pst.npar
pst.svd_data.eigthresh = 1e-6

# -- IES options 
pst.pestpp_options['ies_num_reals'] = num_reals

# par and obs cov mats
pst.pestpp_options['ies_use_empirical_prior'] = False
pst.pestpp_options['parcov'] = 'prior_pcov.jcb'
pst.pestpp_options['obscov'] = 'obs_cov.jcb'

# par and obs ensembles 
pst.pestpp_options['ies_parameter_ensemble'] = 'prior_pe.jcb'
pst.pestpp_options['ies_observation_ensemble'] = 'oe.jcb'

# regularization (see sec. 9.2.3 of pestpp manual)
pst.pestpp_options['ies_use_approx'] = False # expensive form, Eq. 18 from Chen & Oliver (2013)
#pst.pestpp_options['ies_use_approx'] = True # simplified form, Eq. 19 from Chen & Oliver (2013)
pst.pestpp_options['ies_reg_factor'] = 0.05 # very light reg, priority to obs. fit. 

# (automatic) adaptive localization
pst.pestpp_options['ies_autoadaloc'] = False # not effective for this problem
pst.pestpp_options['ies_autoadaloc_sigma_dist'] = 1
pst.pestpp_options['ies_num_threads'] = 20

# Mean-Update Iterations (not effective for this problem)
pst.pestpp_options['ies_n_iter_mean'] = 0

# lambda testing 
# check that : subset_size*n_lambda_mults*n_lambda_scale_fac < n CPU cores
pst.pestpp_options['ies_subset_size'] =  12 
pst.pestpp_options['ies_subset_how'] = 'random'
pst.pestpp_options['ies_lambda_mults'] = [0.1, 1.0, 10.0]
pst.pestpp_options['lambda_scale_fac'] = [0.75, 1.0, 1.25]
pst.pestpp_options['ies_initial_lambda'] = 10000 # (very) high starting value
pst.pestpp_options['ies_accept_phi_fac'] = 1.05
pst.pestpp_options['ies_lambda_inc_fac'] = 1. # so as to avoid extreme values 
pst.pestpp_options['ies_lambda_dec_fac'] = 0.75
pst.pestpp_options['ies_save_lambda_ensembles'] = False
pst.pestpp_options['ies_enforce_bounds'] = True

# reweighting  
pst.pestpp_options['ies_phi_factor_file'] = 'ies_phi_factor.csv'

# bad phi
pst.pestpp_options['ies_bad_phi'] = 1e8  # drop realizations above threshold 
pst.pestpp_options['ies_bad_phi_sigma'] = -0.90 # drop realizations above the 90-th decile
pst.pestpp_options['ies_drop_conflicts'] = False  # not supported with full obs. cov

# run failure
pst.pestpp_options['max_run_fail'] = 1
pst.pestpp_options['overdue_giveup_minutes'] = 100
pst.pestpp_options['overdue_giveup_fac'] = 3

# other options 
pst.pestpp_options['ies_verbose_level'] = 1
pst.control_data.jcosaveitn='jcosaveitn'
pst.pestpp_options['panther_agent_no_ping_timeout_secs'] = 36000

# -----------------------------------------------------------------
# -- Write final PEST control file - End of the story
# -----------------------------------------------------------------

# -- Write .pst file
pst.control_data.noptmax=10
pst.write(mopt.name + '.pst', version=2)

pst.control_data.noptmax=0
pst.write(mopt.name + '_base.pst', version=2)

