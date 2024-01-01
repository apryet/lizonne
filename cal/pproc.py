import os, sys
import shutil
import pandas as pd
import numpy as np
import pyemu
from pymarthe import MartheModel
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
from pymarthe.mfield import MartheField, MartheFieldSeries
from pymarthe.helpers.postprocessing import PestPostProcessing
import matplotlib.pyplot as plt

# --- load model
print('Reading MARTHE model ...')
# Spatial index has to be recomputed, reloading from file leads to error 
configfile = 'cal_lizonne.config'
mm = MartheModel('Lizonne.rma', spatial_index = False)

# setup and clean output directories 
outdir = {'sim_obs':os.path.join('pproc', 'sim_obs'),
           'props':os.path.join('pproc', 'props'),
           'heads':os.path.join('pproc', 'heads')
           }
for d in outdir.values():
    if os.path.exists(d):
        shutil.rmtree(d)
    os.mkdir(d)

_, _, locs = pest_utils.read_config(configfile)

# keep only actual locations (not fluctuations)
loc_list = [d['locnme'] for d in locs if not d['locnme'].endswith('mf') ]

# -- extraction des données simulées
prn_df = marthe_utils.read_prn('historiq.prn')

# -- Préparation de la figure
plt.rc('font', family='serif', size=11)
# -- sim vs obs  
for loc in loc_list:
    # Simulated
    hsim_df = prn_df.loc[:,(slice(None), loc)]
    tmin, tmax = hsim_df.index.min(), hsim_df.index.max()
    ax = hsim_df.plot(figsize= (12,8), c='red', alpha=0.4, lw=0.8)
    # Observed
    hobs_df = marthe_utils.read_obsfile(
                        os.path.join('pest', 'obs', f'{loc}.dat')
                            )[tmin:tmax].plot(ax=ax, ls='none', mew=0.7, legend=False,
                                              marker='+', alpha=0.4, ms=5, color='blue')
    # Ajout du nom du point d'observation
    ax.text(0.01, 0.95, loc, ha='left', va='center',
            transform = ax.transAxes, weight='bold')
    ax.get_figure().savefig(os.path.join(outdir['sim_obs'], f'{loc}.pdf'),dpi=300)


#  hydraulic properties
props  = ['permh','emmca','emmli']
for prop in props:
    mm.load_prop(prop)
    print(f'---> writing {prop} field')
    for ilay in range(mm.nlay):
        ax = mm.prop[prop].plot(layer=ilay, log=True)
        ax.set_title(f'layer n°{ilay+1}', fontsize=12)
        plt.savefig(os.path.join(outdir['props'], f'{prop}_l{ilay+1}.png'))


# -- heads (when available)
print('---> writing heads')
mfs = MartheFieldSeries(mm=mm, chasim='chasim.out')
mfs.load_field('CHARGE')
for istep, mf in mfs.data['CHARGE'].items():
    for ilay in range(mm.nlay):
        ax = mf.plot(layer=ilay, log=False)
        ax.set_title(f'layer n°{ilay+1}', fontsize=12)
        plt.savefig(os.path.join(outdir['heads'], f'head_l{ilay+1}.png'))

