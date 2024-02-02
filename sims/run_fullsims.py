import os,sys
import shutil

tpl_dir = os.path.join('..','tpl_sim')

# number of climate models referenced in ../data/clim_models.xlsx
ncms=17

for i in range(1,ncms+1):
    sim_dir = f'sim_{i:02d}'
    print(f'Replicating {tpl_dir} into {sim_dir} and execute run_fullsim.py')
    if os.path.exists(sim_dir):
       shutil.rmtree(sim_dir)
    shutil.copytree(tpl_dir,sim_dir)
    os.chdir(sim_dir)
    os.system('python run_fullsim.py > run_fullsim.log &')
    os.chdir('..')

