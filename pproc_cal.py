'''
cal post-proc of a single realization
- plot sim-obs tseries 
- plot fields hydraulic properties 
- plot fields of hydraulic heads
- compute days below DOE and save to xlsx 
- plot indicators 
'''

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
plt.rc('font', family='serif', size=11)

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

# -- extract simulation data 
prn_df = marthe_utils.read_prn('historiq.prn')

# -- plot sim-obs time series at each loc 
plt.rc('font', family='serif', size=11)
# -- sim vs obs  
for loc in loc_list:
    # Simulated
    try :
        hsim_df = prn_df.loc[:,(slice(None), loc)]
    except : 
        print(f'{loc} not found in historiq.prn, skipping.')
        continue
    tmin, tmax = hsim_df.index.min(), hsim_df.index.max()
    ax = hsim_df.plot(figsize= (12,8), c='red', alpha=0.4, lw=0.8)
    # Observed
    try : 
        hobs_df = marthe_utils.read_obsfile(
                        os.path.join('pest', 'obs', f'{loc}.dat')
                            )[tmin:tmax].plot(ax=ax, ls='none', mew=0.7, legend=False,
                                              marker='+', alpha=0.4, ms=5, color='blue')
    except : 
        print(f'observation data for {loc} not found, skipping.')
    # Ajout du nom du point d'observation
    _ = ax.text(0.01, 0.95, loc, ha='left', va='center',
            transform = ax.transAxes, weight='bold')
    ax.get_figure().savefig(os.path.join(outdir['sim_obs'], f'{loc}.pdf'),dpi=300)


# -- plot hydraulic properties
props  = ['permh','emmca','emmli']
for prop in props:
    mm.load_prop(prop)
    print(f'---> writing {prop} field')
    for ilay in range(mm.nlay):
        ax = mm.prop[prop].plot(layer=ilay, log=True)
        ax.set_title(f'layer n°{ilay+1}', fontsize=12)
        plt.savefig(os.path.join(outdir['props'], f'{prop}_l{ilay+1}.png'))


# plot head fields (when available)
print('---> writing heads')
mfs = MartheFieldSeries(mm=mm, chasim='chasim.out')
mfs.load_field('CHARGE')
for istep, mf in mfs.data['CHARGE'].items():
    for ilay in range(mm.nlay):
        ax = mf.plot(layer=ilay, log=False)
        ax.set_title(f'layer n°{ilay+1}', fontsize=12)
        plt.savefig(os.path.join(outdir['heads'], f'head_l{ilay+1}.png'))



# --------------  compute and plot indicators 


# load DOE data
doe_file = os.path.join('..','data','DOE','DOE.xlsx')
doe_df = pd.read_excel(doe_file,index_col='id')

# list of stations 
loc_list = ['P7250001','P7270001','P8215010','P8284010']

# load sim data for loc_list 
prn_df = marthe_utils.read_prn('historiq.prn')
# extract discharge data
qsim_df = prn_df.loc[:,(slice(None), loc_list)].loc[:,'Débit_Rivi']

# load obs data for loc_list
qobs_df_list = []
for locnme in loc_list :  
    obsfile = os.path.join('pest', 'obs', f'{locnme}.dat')
    loc_df = marthe_utils.read_obsfile(obsfile)
    loc_df.rename(columns={'value':locnme},inplace=True)
    qobs_df_list.append(loc_df)

qobs_df = pd.concat(qobs_df_list,axis=1)   


# computes number of days under alert.  
def get_qriv_alerts(q_df):
    '''
    Compute days of alerte (alterte, renforce, crise) from q_df
    '''
    # identify periods with flow discharge below threshold levels
    qcrise_df    = q_df.le(doe_df['crise'])
    qrenforce_df = q_df.le(doe_df['renforce'])*(~qcrise_df)
    qalerte_df   = q_df.le(doe_df['alerte'])*(~qrenforce_df)*(~qcrise_df)
    # concatenate all types of alerts
    qalerts_df = pd.concat( [qalerte_df, qrenforce_df, qcrise_df],axis=1,
                         keys=['alerte','renforce','crise'],
                         names=['seuil','station']
                         )
    # summary of days with alerts (number of time steps where threshold values are met)
    qalerts_sum = qalerts_df.astype(int).sum().unstack()
    return(qalerts_df,qalerts_sum)



# subset to intersection with available obs values (for sim-obs comparative purpose)
idx = qobs_df.index.intersection(qsim_df.index)
# keep calibration period only (remove initial warmup)
cal_start = pd.to_datetime('2012-07-31')
idx = idx[idx > cal_start]

ss_obs_df = qobs_df.loc[idx]
ss_sim_df = qsim_df.loc[idx]

# apply function to sim and obs values 
sim_qalerts_df, sim_qalerts_sum = get_qriv_alerts(ss_sim_df)
obs_qalerts_df, obs_qalerts_sum = get_qriv_alerts(ss_obs_df)

qalerts_sum = pd.concat([sim_qalerts_sum,obs_qalerts_sum], keys=['sim','obs'])
simyrs = round((ss_sim_df.index.max()-ss_sim_df.index.min()).days/365.25)
qalerts_pyr = qalerts_sum.div(simyrs).round()
qalerts_pyr.unstack().to_excel(os.path.join('pproc','doe_simobs.xlsx'))

# plot bar plots of days 
fig,axs=plt.subplots(1, 4, sharey=True, figsize=(8,4))
for ax,loc in zip(axs.ravel(),loc_list):
    qalerts_pyr[loc].unstack().T.plot.bar(color=['red','blue'],ax=ax,legend=False)
    ax.set_title(loc)
    ax.set_xlabel('')

axs[0].legend()
axs[0].set_ylabel('Number of days per year')
fig.tight_layout()
fig.savefig(os.path.join('pproc','doe_simobs.pdf'))


# --- plot q records with alerts over simulated period 
sim_qalerts_df, sim_qalerts_sum = get_qriv_alerts(qsim_df)
loc_list = ['P7250001','P7270001','P8215010','P8284010']

def plot_indic(yscale='log'):
    fig,axs=plt.subplots(2,2,sharex=True,figsize=(14,8))
    for lax, locnme in zip(axs.ravel(),loc_list) :  
        # plot setup
        loc_long_name = doe_df.loc[locnme,'nom']
        lax.set_title(f'Débit à la station \"{loc_long_name}\" ({locnme})')
        bax = lax.twinx()
        # plot sim 
        lax.plot(qsim_df[locnme],c='red',alpha=0.8, lw=0.8,label='Sim.')
        # plot obs 
        lax.plot(qobs_df[locnme], ls='none', mew=0.7, marker='+', alpha=0.3, ms=2, color='blue',label='obs')
        lax.set_ylabel('Débit [m$^3$/s]')
        lax.set_xlim(qsim_df[locnme].index.min(),qsim_df[locnme].index.max())
        if yscale=='log':
            lax.set_yscale('log')
        # plot critical levels 
        lax.axhline(y=doe_df.loc[locnme]['alerte'], color='gold', lw=0.8, ls='--')
        lax.axhline(y=doe_df.loc[locnme]['renforce'], color='orange',lw=0.8, ls='--')
        lax.axhline(y=doe_df.loc[locnme]['crise'], color='red',lw=0.8, ls='--')
        # plot bars in background
        bax.bar(sim_qalerts_df.index, sim_qalerts_df.loc[:,('alerte',locnme)].astype(int),
                color='gold',width=1,align='center',alpha=0.5,label='Alerte')
        bax.bar(sim_qalerts_df.index,sim_qalerts_df.loc[:,('renforce',locnme)].astype(int),
                color='orange',align='center',width=1, alpha=0.5,label='Alerte renforcée')
        bax.bar(sim_qalerts_df.index,sim_qalerts_df.loc[:,('crise',locnme)].astype(int),
                color='red',align='center',width=1, alpha=0.5,label='Crise')
        bax.set_ylim(0,1)
        bax.set_yticks([])
        bax.legend()
    fig.tight_layout()
    plt.savefig(os.path.join('pproc', f'doe_simobs_{yscale}.png'),dpi=300)


plot_indic(yscale='log')
plot_indic(yscale='linear')
