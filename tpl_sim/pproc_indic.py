import os, sys, shutil

import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pymarthe import MartheModel
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
from pymarthe.mfield import MartheField, MartheFieldSeries

FFMT = lambda x: "{0:<20.10E} ".format(float(x))
encoding = 'latin-1'

configfile = 'cal_lizonne.config'
_, _, locs = pest_utils.read_config(configfile)
loc_list = [d['locnme'] for d in locs]


# load DOE data
doe_file = os.path.join('..','data','DOE','DOE.xlsx')
doe_df = pd.read_excel(doe_file,index_col='id')

# list of stations 
loc_list = ['P7250001','P7270001','P8215010','P8284010']

# load sim data for loc_list 
prn_df = marthe_utils.read_prn('historiq.prn')
# extract discharge data
qsim_df = prn_df.loc[:,(slice(None), loc_list)].loc[:,'Débit_Rivi']
# remove  steady state time step + 1st time step
# this is required because there is less than 10 days between the first and second time step !
qsim_df = qsim_df.iloc[2:,:]
qsim_df = qsim_df.asfreq('7D')



# load obs data for loc_list
qobs_df_list = []
for locnme in loc_list :  
    obsfile = os.path.join('pest', 'obs', f'{locnme}.dat')
    loc_df = marthe_utils.read_obsfile(obsfile)
    loc_df.rename(columns={'value':locnme},inplace=True)
    qobs_df_list.append(loc_df)

qobs_df = pd.concat(qobs_df_list,axis=1)   


# function to compute and write number of days under alert.  
def get_qriv_alerts(q_df,filename):
    '''
    Compute days of alerte (alterte, renforce, crise) from q_df
    Writes to filename 
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



# apply function to all simulated values (for plotting purpose)
qalerts_df, qalerts_sum = get_qriv_alerts(qsim_df,'qalerts_sum.dat')

# -- apply function to subsets
# process obs to be comparable to simulated values
# 10-day average, labeled at end of period, supposed to correspond to discharge data in the prn file 
date_start = pd.to_datetime('2003-07-31')
qobs_agg_df = qobs_df.groupby(pd.Grouper(freq='7D',label='right',origin=date_start)).agg('mean')
idx = qobs_agg_df.index.intersection(qsim_df.index)

# subset to intersection
ss_obs_df = qobs_agg_df.loc[idx]
ss_sim_df = qsim_df.loc[idx]

# apply function to sim and obs values 
sim_qalerts_df, sim_qalerts_sum = get_qriv_alerts(ss_sim_df,'sim_qalerts_sum.dat')
obs_qalerts_df, obs_qalerts_sum = get_qriv_alerts(ss_obs_df,'obs_qalerts_sum.dat')


# load sim data for loc_list 
prncal_df = marthe_utils.read_prn('historiq.prn.cal')
# extract discharge data
qsimcal_df = prn_df.loc[:,(slice(None), loc_list)].loc[:,'Débit_Rivi']


# plot 
for locnme in loc_list :  
    # plot setup
    plt.rc('font', family='serif', size=11)
    fig,lax = plt.subplots(figsize=(8,4))
    loc_long_name = doe_df.loc[locnme,'nom']
    lax.set_title(f'Débit à la station \"{loc_long_name}\" ({locnme})')
    bax = lax.twinx()
    # plot sim 
    lax.plot(qsim_df[locnme],c='red',alpha=0.8, lw=0.8,label='Sim.')
    lax.plot(qsimcal_df[locnme],ls='none', marker='+',c='black',alpha=0.8, lw=0.8,label='Sim.')
    # plot obs 
    lax.plot(qobs_df[locnme], ls='none', mew=0.7, marker='+', alpha=0.4, ms=5, color='blue',label='obs')
    lax.set_ylabel('Débit [m$^3$/s]')
    lax.set_yscale('log')
    # plot critical levels 
    lax.axhline(y=doe_df.loc[locnme]['alerte'], color='gold', lw=0.8, ls='--')
    lax.axhline(y=doe_df.loc[locnme]['renforce'], color='orange',lw=0.8, ls='--')
    lax.axhline(y=doe_df.loc[locnme]['crise'], color='red',lw=0.8, ls='--')
    # plot bars in background (bar widths of 7 days)
    bax.bar(qalerts_df.index, qalerts_df.loc[:,('alerte',locnme)].astype(int),
            color='gold',width=7,align='center',alpha=0.5,label='Alerte')
    bax.bar(qalerts_df.index,qalerts_df.loc[:,('renforce',locnme)].astype(int),
            color='orange',align='center',width=7, alpha=0.5,label='Alerte renforcée')
    bax.bar(qalerts_df.index,qalerts_df.loc[:,('crise',locnme)].astype(int),
            color='red',align='center',width=7, alpha=0.5,label='Crise')
    bax.set_ylim(0,1)
    bax.set_yticks([])
    bax.legend()
    fig.tight_layout()
    plt.savefig(os.path.join('pproc','indic', f'{locnme}.png'),dpi=300)


