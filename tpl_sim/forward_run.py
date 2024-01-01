

print('███████╗ ██████╗ ██████╗ ██╗    ██╗ █████╗ ██████╗ ██████╗     ██████╗ ██╗   ██╗███╗   ██╗')
print('██╔════╝██╔═══██╗██╔══██╗██║    ██║██╔══██╗██╔══██╗██╔══██╗    ██╔══██╗██║   ██║████╗  ██║')
print('█████╗  ██║   ██║██████╔╝██║ █╗ ██║███████║██████╔╝██║  ██║    ██████╔╝██║   ██║██╔██╗ ██║')
print('██╔══╝  ██║   ██║██╔══██╗██║███╗██║██╔══██║██╔══██╗██║  ██║    ██╔══██╗██║   ██║██║╚██╗██║')
print('██║     ╚██████╔╝██║  ██║╚███╔███╔╝██║  ██║██║  ██║██████╔╝    ██║  ██║╚██████╔╝██║ ╚████║')
print('╚═╝      ╚═════╝ ╚═╝  ╚═╝ ╚══╝╚══╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝')


# -- Import basic python modules
import os
import pyemu
import multiprocessing as mp
import numpy as np
import pandas as pd
import pymarthe
from pymarthe.utils import marthe_utils, pest_utils

import helpers

def main():
	# -- Run preproc functions
	helpers.opt_parfacmul()
	# -- Run model from .config file
	pymarthe.utils.pest_utils.run_from_config("opt_lizonne.config")
	# -- Run extra functions
	helpers.opt_postproc()

# -- Launch forward run
if __name__ == '__main__':
	mp.freeze_support()
	main()

