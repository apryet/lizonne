import os,sys
import shutil
import pandas as pd 

simyrs = pd.read_csv('sim_yrs.csv')

for i in range(1,ncms+1):
    sim_dir = f'sim_{i:02d}'
    print(f'Replicating {tpl_dir} into {sim_dir} and execute run_fullsim.py')
    if os.path.exists(sim_dir):
       shutil.rmtree(sim_dir)
    shutil.copytree(tpl_dir,sim_dir)
    os.chdir(sim_dir)
    os.system('python run_fullsim.py > run_fullsim.log &')
    os.chdir('..')

