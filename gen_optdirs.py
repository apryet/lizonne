import os, sys, shutil
import pandas as pd
import numpy as np

from pymarthe import MartheModel
from pymarthe.mfield import MartheField, MartheFieldSeries
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils


print('OPTIMIZATION ROADMAP:')
print('--------------------------------------------------------------------------------')
print('DONE: pproc_ies from master_ies copies model files + chasim.out.cal to sim dir')
print('--------------------------------------------------------------------------------')
print('DONE: sim_setup.py')
print('--------------------------------------------------------------------------------')
print('DONE: run clim_analysis(_ens).py.') # _ens if multiple scenarios are considered
print('--------------------------------------------------------------------------------')
print('NOW:  generation of simulation dirs for PESTPP-MOU')
print('--------------------------------------------------------------------------------')

# ===========================================================
# -- generate simulation dirs from template 
# ===========================================================

# template dir
'''
The template sim is run over the whole simulation period (1950-2100) at the daily time step 
for Gardenia, weekly time step for hydrodynamic. 
Withdrawals at median historical levels 
Heads saved on August, 1st
'''
tpl_dir = 'tpl_sim'

print('----------------------------------------------------------------------------')
print(f'INFO :  {tpl_dir} is the template simulation dir.')
print('WARNING #1: make sure to have updated modelfiles (Lizonne.*) + chasim.out.cal  !!!')
print('----------------------------------------------------------------------------')
print('WARNING #2: make sure to have run clim_analysis_ens.py !!!')
print('----------------------------------------------------------------------------')

# get simulation years from table generated by clim_analysis.py
simyrs_file = 'simyrs.csv'
print(f'Reading {simyrs_file}....')
qyrs = pd.read_csv(simyrs_file,index_col=[0,1])
qyrs['quant'] = qyrs.index.get_level_values(1)
qyrs['sim_dir'] = qyrs.apply(lambda x: f'sim_{x.quant}_{x.year}_{x.cm}' ,axis=1)

# keep only median years
#qyrs = qyrs.loc[(slice(None),['Q10','Q90']),:]

# generate simulation dir for each year of interest 
# NOTE : pastp is now updated on runtime during preproc
for cm_id, year, sim_dir in zip(qyrs.cm, qyrs.year, qyrs.sim_dir):
    print(f'Generating {sim_dir}...')
    # remove if exists 
    if os.path.exists(sim_dir):
           shutil.rmtree(sim_dir)
    # copy from template 
    shutil.copytree(tpl_dir,sim_dir)
    # add climate model selection file 
    cm_df = pd.DataFrame({'parnme':['cm'],'defaultvalue':cm_id})
    cm_par_file = os.path.join(sim_dir,'cm.dat')
    pest_utils.write_mlp_parfile(cm_par_file,cm_df)
    # update initial heads 
    fullsim_dir = os.path.join('sims',f'sim_{cm_id:02d}')
    mm = MartheModel(os.path.join(fullsim_dir,'Lizonne.rma'))
    chasim  = os.path.join(fullsim_dir,'chasim.out.full')
    mfs = MartheFieldSeries(mm=mm, chasim=chasim)
    mfs.load_field('CHARGE')
    # write initial heads
    sim_start = pd.to_datetime(str(year-2)+'-08-01') 
    istep = int(np.argwhere(mm.mldates==sim_start))
    hinit= mfs.data['CHARGE'][istep]
    hinit.write_data('Lizonne.charg')
    os.chdir(sim_dir)
    os.system('python mou_setup.py > mou_setup.log &')
    os.chdir('..')
    

