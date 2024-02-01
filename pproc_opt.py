# ---- Import usefull modules
import os, sys
import shutil
import platform
import subprocess as sp
import pandas as pd
import numpy as np
import pyemu
from pymarthe import MartheModel
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
from pymarthe.mfield import MartheField, MartheFieldSeries
import matplotlib.pyplot as plt
import matplotlib 


# clean and (re)-build directory tree
shutil.rmtree('pproc')
os.mkdir('pproc')
for d in ['haqpump','qriv','heads']:
    os.mkdir(os.path.join('pproc',d))

mm = MartheModel('Lizonne.rma', spatial_index = True)

# ---------------------------------------------
# plot pumping factors 
# ---------------------------------------------

# --- aquifer 
fackmi, aqfacvals = pest_utils.parse_mlp_parfile(os.path.join('pest','par','aqpumpfac.dat'),
                                               keys=['boundname','layer','istep'],
                                               value_col=1, btrans='none')



aqfacvals.index = pd.MultiIndex.from_frame(pd.DataFrame({
    'layer':fackmi.get_level_values('layer'),
    'date':mm.mldates[fackmi.get_level_values('istep')].strftime('%Y-%m')
    }))


# plot 
fig,ax = plt.subplots(1,1,figsize=(10,6))
aqfacvals.unstack().T.plot(ax=ax,kind='bar')
ax.axhline(y=1,color='grey', ls='--')
fig.savefig(os.path.join('pproc','aqfacvals.png'),dpi=300)


# --- surface
fackmi, rivfacvals = pest_utils.parse_mlp_parfile(os.path.join('pest','par','rivpumpfac.dat'),
                                               keys=['prefix','aff_r','istep'], 
                                               value_col=1, btrans='none')


rivfacvals.index = pd.MultiIndex.from_frame(pd.DataFrame({
    'aff_r':fackmi.get_level_values('aff_r'),
    'date':mm.mldates[fackmi.get_level_values('istep')].strftime('%Y-%m')
    }))


import geopandas as gpd
simriv_shp = os.path.join('gis','sim_riv.shp')
simriv_gdf = gpd.read_file(simriv_shp)
simriv_gdf.set_index(simriv_gdf.val.astype(int),inplace=True)
simriv_gdf = simriv_gdf.merge(rivfacvals.unstack(),left_index=True, right_index=True)
plot_dates = rivfacvals.index.get_level_values('date').unique()


plt.rc('font', family='serif', size=11)
fig,axs = plt.subplots(3,4,figsize=(10,6))
for d,ax in zip(plot_dates,axs.ravel()):
    simriv_gdf.plot(ax=ax,column=d,vmin=0,vmax=10)
    ax.set_title(d)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.2, 0.02, 0.60])
fig.colorbar(ax.get_children()[0],cax=cbar_ax,label='factor')
fig.savefig(os.path.join('pproc','qrifac.png'),dpi=300)

# ---------------------------------------------
# plot constraints 
# ---------------------------------------------

# load sim data  
prn_df = marthe_utils.read_prn('historiq.prn')

# ---- plot water levels in pumping cells 

# extract discharge data
aqpump_lim_df = pd.read_csv('aqpump_lim.csv',index_col=0)
aqpump_names = aqpump_lim_df.index

# -- Préparation de la figure
plt.rc('font', family='serif', size=11)
# -- Construction des graphiques pour chaque points d'observation
for loc in aqpump_names:
    fig,ax = plt.subplots(figsize=(8,4))
    # plot sim
    hsim_df = prn_df.loc[:,(slice(None), loc)]
    ax = hsim_df.plot(ax=ax,figsize= (12,8), c='red', alpha=0.4, lw=0.8)
    # plot constraint value 
    ax.axhline(y=aqpump_lim_df.loc[loc]['hmin'], color='red',lw=0.8, ls='--')
    # print loc name 
    ax.text(0.01, 0.95, loc, ha='left', va='center',
            transform = ax.transAxes, weight='bold')
    fig.savefig(os.path.join('pproc','haqpump', f'{loc}.png'))



# ---- plot river discharge at gaging stations  

loc_list = ['P8284010', 'P8215010', 'P7250001', 'P7270001']

# extract discharge data
qsim_df = prn_df.loc[:,(slice(None), loc_list)].loc[:,'Débit_Rivi']

doe_file = os.path.join('..','data','DOE','DOE.xlsx')
doe_df = pd.read_excel(doe_file,index_col='id')
loc_list = doe_df.index 

# --- compute days of alerte (alterte, renforce, crise) from qsim_df
#  periods with flow discharge below threshold levels
qcrise_df    = qsim_df.le(doe_df['crise'])
qrenforce_df = qsim_df.le(doe_df['renforce'])*(~qcrise_df)
qalerte_df   = qsim_df.le(doe_df['alerte'])*(~qrenforce_df)*(~qcrise_df)
# concatenate all types of alerts
qalerts_df = pd.concat( [qalerte_df, qrenforce_df, qcrise_df],axis=1,
                     keys=['alerte','renforce','crise'],
                     names=['seuil','station']
                     )
# summary of days with alerts (number of time steps where threshold values are met)
qalerts_sum = qalerts_df.astype(int).sum().unstack()

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
    lax.legend(loc='upper left')
    # plot critical levels 
    lax.axhline(y=doe_df.loc[locnme]['alerte'], color='gold', lw=0.8, ls='--')
    lax.axhline(y=doe_df.loc[locnme]['renforce'], color='orange',lw=0.8, ls='--')
    lax.axhline(y=doe_df.loc[locnme]['crise'], color='red',lw=0.8, ls='--')
    # plot bars in background (bar widths of 1 days)
    bax.bar(qalerts_df.index, qalerts_df.loc[:,('alerte',locnme)].astype(int),
            color='gold',width=1,align='center',alpha=0.5,label='Alerte')
    bax.bar(qalerts_df.index,qalerts_df.loc[:,('renforce',locnme)].astype(int),
            color='orange',align='center',width=1, alpha=0.5,label='Alerte renforcée')
    bax.bar(qalerts_df.index,qalerts_df.loc[:,('crise',locnme)].astype(int),
            color='red',align='center',width=1, alpha=0.5,label='Crise')
    bax.set_ylim(0,1)
    bax.set_yticks([])
    bax.legend(loc='upper right')
    fig.tight_layout()
    plt.savefig(os.path.join('pproc','qriv', f'{locnme}.png'),dpi=300)

# ---------------------------------------------
# plot heads
# ---------------------------------------------

# -- Heads
print('---> Heads')
mfs = MartheFieldSeries(mm=mm, chasim='chasim.out')
mfs.load_field('CHARGE')
for istep, mf in mfs.data['CHARGE'].items():
    mm.mldates[istep].strftime
    folder = os.path.join('pproc','heads', f'head_i{istep}')
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)
    for ilay in range(mm.nlay):
        ax = mf.plot(layer=ilay, log=False)
        ax.set_title(f'layer n°{ilay+1}', fontsize=12)
        plt.savefig(os.path.join(folder, f'head_l{ilay+1}.png'))

