import os 
import argparse
import pyemu

parser = argparse.ArgumentParser(description='PEST parallel run launcher')
parser.add_argument('--algo', type = str, default='ies')
parser.add_argument('--pst', type = str, default='cal_lizonne.pst')
parser.add_argument('--tpl_dir', type = str, default='cal')
algo, pst, tpl_dir  = [v for _, v in parser.parse_args()._get_kwargs()]

# clear worker directories
os.system('rm -r workers/*')

# start master + loopt workers
pyemu.helpers.start_workers(tpl_dir,f'pestpp-{algo}',pst,num_workers=12,
                              worker_root= 'workers',cleanup=False,
                                master_dir=f'master_{algo}')

