import os, sys
import shutil
import platform
import subprocess as sp
import pandas as pd
import geopandas as gpd
import numpy as np
import pyemu
from pymarthe import MartheModel
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
from pymarthe.mfield import MartheField, MartheFieldSeries
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib 
import kneed 


plt.rc('font', family='serif', size=11)


# ----------------------------------------------------
# pproc mou  : last generation pareto with fac1 config
# ----------------------------------------------------


def plot_pareto(master_dir, gen, label, color, marker, ax, is_feasible=False):
    pasum_df = pd.read_csv(os.path.join(master_dir,'mou_lizonne.pareto.archive.summary.csv'))
    if is_feasible:
        feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
    else:  
        feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1, axis=1),:]
    df = feas_front_df.loc[feas_front_df.generation==gen,:]
    # non-dominated realizations (pareto)
    pax = ax.scatter(df.loc[:,'deficit_tot'],df.loc[:,'tot_pump'],
               marker=marker,s=60,color=color,label=label,alpha=0.7)
    # knee point (maximum curvature along pareto front)
    sdf = df.sort_values('deficit_tot')
    kn = kneed.KneeLocator(
            sdf.loc[:,'deficit_tot'],
            sdf.loc[:,'tot_pump'],
            curve='concave',
            direction='increasing',
            interp_method='interp1d',
        )
    #kax = ax.scatter(kn.knee,kn.knee_y,marker='D', edgecolor='black', color="none",s=60, label='Knee point')
    kax=None
    # historical reference (fac=1)
    fac1_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac1.base.rei'))
    fac1_pump = fac1_rei.loc['tot_pump','modelled']
    fac1_deficit = fac1_rei.loc['deficit_tot','modelled']
    f1ax = ax.scatter(fac1_deficit,fac1_pump,marker="+",color=color,s=80,label='Factor=1')
    return(pax, kax, f1ax) 


marker_dic = {'Q10':'^','Q50':'o','Q90':'v'}

#master_dirs = ['master_sim_Q10_1987_1','master_sim_Q50_2002_3','master_sim_Q50_2090_4','master_sim_Q90_1977_3']
master_dirs = ['master_sim_Q50_2002_3','master_sim_Q50_2090_4']
master_dirs = ['master_sim_Q10_1987_1','master_sim_Q50_2002_3','master_sim_Q90_1977_3']
master_dirs = ['master_sim_Q10_1987_1','master_sim_Q50_2002_3','master_sim_Q90_1977_3','master_sim_Q10_2086_2','master_sim_Q50_2090_4','master_sim_Q90_2092_4']

fig,ax = plt.subplots(1,1,figsize=(5,5))
for master_dir in master_dirs:
    year = int(master_dir.split('_')[-2])
    qt = master_dir.split('_')[-3]
    prefix = 'hist.' if year < 2005 else 'fut.'
    color = 'grey' if year < 2005 else 'tan'
    marker=marker_dic[qt]
    pax, kax, f1ax = plot_pareto(master_dir, gen=9, label=f'{prefix}-{qt}', color=color, marker=marker, ax=ax)


q90_label = Line2D([0], [0], label='Q90', marker=None, linestyle= '')
q50_label = Line2D([0], [0], label='Q50', marker=None, linestyle= '')
q10_label = Line2D([0], [0], label='Q10', marker=None, linestyle= '')

hist_q10_marker = Line2D([0], [0], label='', marker='^',
    markersize=10, markeredgecolor='black', markerfacecolor='grey', linestyle='')
hist_q50_marker = Line2D([0], [0], label='', marker='o',
    markersize=10, markeredgecolor='black', markerfacecolor='grey', linestyle='')
hist_q90_marker = Line2D([0], [0], label='', marker='v',
    markersize=10, markeredgecolor='black', markerfacecolor='grey', linestyle='')

fut_q10_marker = Line2D([0], [0], label='', marker='^',
    markersize=10, markeredgecolor='black', markerfacecolor='tan', linestyle='')
fut_q50_marker = Line2D([0], [0], label='', marker='o',
    markersize=10, markeredgecolor='black', markerfacecolor='tan', linestyle='')
fut_q90_marker = Line2D([0], [0], label='', marker='v',
    markersize=10, markeredgecolor='black', markerfacecolor='tan', linestyle='')


ax.legend(handles=[q90_label, q50_label, q10_label,
                    hist_q90_marker,hist_q50_marker,hist_q10_marker,
                  fut_q90_marker,fut_q50_marker,fut_q10_marker],
          columnspacing=1., handletextpad=0, handlelength=1.,borderaxespad=0.5,
          title = '              Hist. Fut. ',
     loc='lower right', ncols=3)

ax.set_ylim([0,2e7])
ax.set_xlabel('Total river deficit [m$^3$]')
ax.set_ylabel('Total pumping [m$^3$]')
fig.tight_layout()
fig.savefig(os.path.join('figs',f'comp_pareto.pdf'),dpi=300)



