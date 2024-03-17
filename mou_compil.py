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

plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54
mdcol_width = 14/2.54
dbcol_width = 19/2.54

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
    pax = ax.scatter(df.loc[:,'deficit_tot']*1e-6,df.loc[:,'tot_pump']*1e-6,
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
    kax = ax.scatter(kn.knee*1e-6,kn.knee_y*1e-6,marker=marker, edgecolor='darkred', alpha=0.8, color=color,s=60, label='Knee point')
    # no pumping (fac=0)
    fac0_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac0.base.rei'))
    fac0_pump = fac0_rei.loc['tot_pump','modelled']*1e-6
    fac0_deficit = fac0_rei.loc['deficit_tot','modelled']*1e-6
    f0ax = ax.scatter(fac0_deficit,fac0_pump,marker=marker,edgecolor='black',color=color,s=60,label='Factor=0')
    # historical reference (fac=1)
    fac1_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac1.base.rei'))
    fac1_pump = fac1_rei.loc['tot_pump','modelled']*1e-6
    fac1_deficit = fac1_rei.loc['deficit_tot','modelled']*1e-6
    f1ax = ax.scatter(fac1_deficit,fac1_pump,marker=marker,edgecolor='black',color=color,alpha=0.8,s=60,label='Factor=1')


marker_dic = {'Q5':'^','Q50':'o','Q95':'v'}

#master_dirs = ['master_sim_Q50_1978_5','master_sim_Q50_2076_1','master_sim_Q5_1981_2','master_sim_Q5_2079_11']

#master_dirs = ['sim_Q50_1978_5','sim_Q50_2076_1','sim_Q5_1981_2','sim_Q5_2079_11','sim_Q95_2003_5','sim_Q95_2069_17']

master_dirs=['master_sim_Q5_1981_2','master_sim_Q5_2079_11','master_sim_Q50_1978_5','master_sim_Q50_2076_1','master_sim_Q95_2003_5','master_sim_Q95_2069_17']

fig,ax = plt.subplots(1,1,figsize=(sgcol_width,sgcol_width))
for master_dir in master_dirs:
    year = int(master_dir.split('_')[-2])
    qt = master_dir.split('_')[-3]
    prefix = 'hist.' if year < 2005 else 'fut.'
    color = 'grey' if year < 2005 else 'tan'
    marker=marker_dic[qt]
    plot_pareto(master_dir, gen=10, label=f'{prefix}-{qt}', color=color, marker=marker, ax=ax)



# lines 

fac0_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac0.base.rei'))
fac0_pump = fac0_rei.loc['tot_pump','modelled']*1e-6
fac0_deficit = fac0_rei.loc['deficit_tot','modelled']*1e-6
f0lax = ax.axhline(fac0_pump,color='black',lw=0.5,ls='--')
# historical reference (fac=1)
fac1_rei = pyemu.pst_utils.read_resfile(os.path.join(master_dir,'mou_lizonne_fac1.base.rei'))
fac1_pump = fac1_rei.loc['tot_pump','modelled']*1e-6
fac1_deficit = fac1_rei.loc['deficit_tot','modelled']*1e-6
f1lax = ax.axhline(fac1_pump,color='black',lw=0.5,ls='--')

ax.text(4.3,fac0_pump+0.2,'$f=0$',c='black')
ax.text(6.7,fac1_pump+0.2,'$f=1$',c='black')

# legend
q95_label = Line2D([0], [0], label='Q95', markersize=1, marker=None, linestyle= '')
q50_label = Line2D([0], [0], label='Q50', markersize=1, marker=None, linestyle= '')
q05_label = Line2D([0], [0], label='Q05', markersize=1, marker=None, linestyle= '')

hist_q05_marker = Line2D([0], [0], label='', marker='^',
    markersize=8, markeredgecolor='black', markerfacecolor='grey', linestyle='')
hist_q50_marker = Line2D([0], [0], label='', marker='o',
    markersize=8, markeredgecolor='black', markerfacecolor='grey', linestyle='')
hist_q95_marker = Line2D([0], [0], label='', marker='v',
    markersize=8, markeredgecolor='black', markerfacecolor='grey', linestyle='')

fut_q05_marker = Line2D([0], [0], label='', marker='^',
    markersize=8, markeredgecolor='black', markerfacecolor='tan', linestyle='')
fut_q50_marker = Line2D([0], [0], label='', marker='o',
    markersize=8, markeredgecolor='black', markerfacecolor='tan', linestyle='')
fut_q95_marker = Line2D([0], [0], label='', marker='v',
    markersize=8, markeredgecolor='black', markerfacecolor='tan', linestyle='')


lg = ax.legend(handles=[q95_label, q50_label, q05_label,
                    hist_q95_marker,hist_q50_marker,hist_q05_marker,
                  fut_q95_marker,fut_q50_marker,fut_q05_marker],
          columnspacing=1., handlelength=1.,
          title = '        Hist.  Fut.',
     loc='lower right',fontsize=7, ncols=3)

plt.setp(lg.get_title(),fontsize=8)

ax.set_ylim([0,32])
ax.set_xlim([0,9])
ax.set_xlabel('Total River Deficit [Mm$^3$]')
ax.set_ylabel('Total Pumping [Mm$^3$]')
fig.tight_layout()
fig.savefig(os.path.join('figs',f'comp_pareto.pdf'),dpi=300)

# --------------------------------------------------
# compare pre-defined risk, with risk as an objective 

fig,ax = plt.subplots(1,1,figsize=(sgcol_width,sgcol_width))
pasum_df = pd.read_csv(os.path.join('master_simrisk_Q50_1978_5','mou_lizonne.pareto.archive.summary.csv'))
feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
gen=feas_front_df.generation.max()
df = feas_front_df.loc[feas_front_df.generation==gen,:]
# non-dominated realizations (pareto)
pax = ax.scatter(df.loc[:,'deficit_tot']*1e-6,df.loc[:,'tot_pump']*1e-6,c=df.loc[:,'_risk_'],
                 cmap='seismic_r',vmin=0,vmax=1, marker='o',ec='black',s=60,alpha=1)

cax = fig.colorbar(pax, label='Reliability',orientation='horizontal',
             cax=ax.inset_axes((0.4, 0.30, 0.5, 0.05)))

pasum_df = pd.read_csv(os.path.join('master_sim_Q50_1978_5','mou_lizonne.pareto.archive.summary.csv'))
feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
gen=feas_front_df.generation.max()
df = feas_front_df.loc[feas_front_df.generation==gen,:]
feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
pax = ax.scatter(df.loc[:,'deficit_tot']*1e-6,df.loc[:,'tot_pump']*1e-6,c='grey',
                  marker='o',s=60,alpha=0.9, label='Specified reliability at 0.66')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles,labels=labels,loc='upper right',bbox_to_anchor=(1,0.125),fontsize=8)

ax.set_xlabel('Total River Deficit [Mm$^3$]')
ax.set_ylabel('Total Pumping [Mm$^3$]')
ax.set_ylim([0,32])
ax.set_xlim([0,9])
fig.tight_layout()
fig.savefig(os.path.join('figs',f'pareto_risk_as_objective.pdf'),dpi=300)

