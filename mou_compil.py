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
               marker=marker,s=60,color=color,label=label,alpha=0.8)
    # knee point (maximum curvature along pareto front)
    sdf = df.sort_values('deficit_tot')
    kn = kneed.KneeLocator(
            sdf.loc[:,'deficit_tot'],
            sdf.loc[:,'tot_pump'],
            curve='concave',
            direction='increasing',
            interp_method='interp1d',
        )
    kax = ax.scatter(kn.knee,kn.knee_y,marker=marker, edgecolor='darkred', alpha=0.8, color=color,s=60, label='Knee point')
    #kax=None
    # no pumping (fac=0)
    fac0_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac0.base.rei'))
    fac0_pump = fac0_rei.loc['tot_pump','modelled']
    fac0_deficit = fac0_rei.loc['deficit_tot','modelled']
    f0ax = ax.scatter(fac0_deficit,fac0_pump,marker=marker,edgecolor='black',color=color,s=60,label='Factor=0')
    # historical reference (fac=1)
    fac1_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac1.base.rei'))
    fac1_pump = fac1_rei.loc['tot_pump','modelled']
    fac1_deficit = fac1_rei.loc['deficit_tot','modelled']
    f1ax = ax.scatter(fac1_deficit,fac1_pump,marker=marker,edgecolor='black',color=color,alpha=0.8,s=60,label='Factor=1')
    return(pax, kax, f0ax, f1ax) 


marker_dic = {'Q5':'^','Q50':'o','Q95':'v'}


master_dirs = ['master_sim_Q50_1978_5','master_sim_Q50_2076_1','master_sim_Q5_1981_2','master_sim_Q5_2079_11']

fig,ax = plt.subplots(1,1,figsize=(5,5))
for master_dir in master_dirs:
    year = int(master_dir.split('_')[-2])
    qt = master_dir.split('_')[-3]
    prefix = 'hist.' if year < 2005 else 'fut.'
    color = 'grey' if year < 2005 else 'tan'
    marker=marker_dic[qt]
    pax, kax, f0ax, f1ax = plot_pareto(master_dir, gen=9, label=f'{prefix}-{qt}', color=color, marker=marker, ax=ax)



# lines 

fac0_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac0.base.rei'))
fac0_pump = fac0_rei.loc['tot_pump','modelled']
fac0_deficit = fac0_rei.loc['deficit_tot','modelled']
f0lax = ax.axhline(fac0_pump,color='black',lw=0.5,ls='--')
# historical reference (fac=1)
fac1_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac1.base.rei'))
fac1_pump = fac1_rei.loc['tot_pump','modelled']
fac1_deficit = fac1_rei.loc['deficit_tot','modelled']
f1lax = ax.axhline(fac1_pump,color='black',lw=0.5,ls='--')

ax.text(1e5,fac0_pump+2e5,'$f=0$',c='black')
ax.text(1e5,fac1_pump+2e5,'$f=1$',c='black')

# legend
q95_label = Line2D([0], [0], label='Q95', marker=None, linestyle= '')
q50_label = Line2D([0], [0], label='Q50', marker=None, linestyle= '')
q05_label = Line2D([0], [0], label='Q05', marker=None, linestyle= '')

hist_q05_marker = Line2D([0], [0], label='', marker='^',
    markersize=10, markeredgecolor='black', markerfacecolor='grey', linestyle='')
hist_q50_marker = Line2D([0], [0], label='', marker='o',
    markersize=10, markeredgecolor='black', markerfacecolor='grey', linestyle='')
hist_q95_marker = Line2D([0], [0], label='', marker='v',
    markersize=10, markeredgecolor='black', markerfacecolor='grey', linestyle='')

fut_q05_marker = Line2D([0], [0], label='', marker='^',
    markersize=10, markeredgecolor='black', markerfacecolor='tan', linestyle='')
fut_q50_marker = Line2D([0], [0], label='', marker='o',
    markersize=10, markeredgecolor='black', markerfacecolor='tan', linestyle='')
fut_q95_marker = Line2D([0], [0], label='', marker='v',
    markersize=10, markeredgecolor='black', markerfacecolor='tan', linestyle='')


ax.legend(handles=[q95_label, q50_label, q05_label,
                    hist_q95_marker,hist_q50_marker,hist_q05_marker,
                  fut_q95_marker,fut_q50_marker,fut_q05_marker],
          columnspacing=1., handletextpad=0, handlelength=1.,borderaxespad=0.5,
          title = '              Hist. Fut. ',
     loc='lower right', ncols=3)

ax.set_ylim([0,3e7])
ax.set_xlim([0,1.2e7])
ax.set_xlabel('Total river deficit [m$^3$]')
ax.set_ylabel('Total pumping [m$^3$]')
fig.tight_layout()
fig.savefig(os.path.join('figs',f'comp_pareto.pdf'),dpi=300)

# --------------------------------------------------
# compare pre-defined risk, with risk as an objective 

'''
fig,ax = plt.subplots(1,1,figsize=(6,5))
pasum_df = pd.read_csv(os.path.join('pproc_master_simrisk_Q50_1986_7','mou_lizonne.pareto.archive.summary.csv'))
feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
gen=feas_front_df.generation.max()
df = feas_front_df.loc[feas_front_df.generation==gen,:]
# non-dominated realizations (pareto)
pax = ax.scatter(df.loc[:,'deficit_tot'],df.loc[:,'tot_pump'],c=df.loc[:,'_risk_'],
                 cmap='seismic_r',vmin=0,vmax=1, marker='o',ec='black',s=60,label='risk as an objective',alpha=1)

fig.colorbar(pax, label='Reliability')

pasum_df = pd.read_csv(os.path.join('pproc_master_sim_Q50_1986_7','mou_lizonne.pareto.archive.summary.csv'))
feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
gen=feas_front_df.generation.max()
df = feas_front_df.loc[feas_front_df.generation==gen,:]
feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
pax = ax.scatter(df.loc[:,'deficit_tot'],df.loc[:,'tot_pump'],c='none',
                  marker='o',ec='darkgreen',s=60,label='opt_risk=0.66',alpha=1)

ax.set_xlabel('Total river deficit [m$^3$]')
ax.set_ylabel('Total pumping [m$^3$]')
fig.tight_layout()
fig.savefig(os.path.join('figs',f'pareto_risk_as_objective.png'),dpi=300)

'''
