import argparse
import os 
import subprocess

parser = argparse.ArgumentParser(description='PEST parallel run launcher')
parser.add_argument('--tpl_dir', type = str, default='cal_lizonne_v5')
tpl_dir,  = [v for _, v in parser.parse_args()._get_kwargs()]

print(tpl_dir)

# compress template dir
os.system(f'tar --exclude="*.jcb"  -czvf {tpl_dir}.tar.gz {tpl_dir}')

# start obelix
os.system(f'scp {tpl_dir}.tar.gz apryet@obelix.ensegid.fr:/srv/common/eauxscars/')

# start numerix
os.system(f'scp {tpl_dir}.tar.gz apryet@numerix.ensegid.fr:/srv/common/eauxscars/')

