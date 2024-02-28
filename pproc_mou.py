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
import matplotlib 
import kneed 
# plot settings
plt.rc('font', family='serif', size=11)
color_dic = {'fac1':'tan','opt':'royalblue'}

# clean and (re)-build directory tree
shutil.rmtree('pproc')
os.mkdir('pproc')
for d in ['haqpump','qriv','heads']:
    os.mkdir(os.path.join('pproc',d))


# TODO for KNEE POINT REALIZATION 
# -> save pst with knee point parameters
# -> run Marthe model 
# -> plot constraints on head for knee point
# -> plot discharge rates for knee point 
# -> compute number of days alerte / alerte renforcée / crise 


# plot pareto over all generations
# plot annotated pareto of last generation
# insert obj of fac 1 in pareto

# analyze dv for knee point
mm = MartheModel('Lizonne.rma', spatial_index = True)

# load historical reference (fac=1)
sim_dir = os.path.split(os.getcwd())[-1].split('master_')[1] 
fac1_rei = pyemu.pst_utils.read_resfile(os.path.join('..', sim_dir,'mou_lizonne_fac1.base.rei'))
fac1_pump = fac1_rei.loc['tot_pump','modelled']
fac1_deficit = fac1_rei.loc['deficit_tot','modelled']

# ---------------------------------------------
# pproc mou  : pareto
# ---------------------------------------------

pst = pyemu.Pst('mou_lizonne.pst')

# pareto archive summary
pasum_df = pd.read_csv('mou_lizonne.pareto.archive.summary.csv')
feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
#feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1, axis=1),:]

ngen = feas_front_df.generation.unique().shape[0]
#ngen=11

cmap = matplotlib.colormaps.get_cmap('gist_heat').reversed()

fig,ax = plt.subplots(1,1,figsize=(6,6))
objs = pst.pestpp_options["mou_objectives"].split(',')
for gen in range(ngen):
    df = feas_front_df.loc[feas_front_df.generation==gen,:]
    sc = ax.scatter(df.loc[:,objs[0]],df.loc[:,objs[1]],c=df.loc[:,'generation'],vmin=0,vmax=ngen,
               cmap=cmap,marker="o", 
               label=f'gen. {gen}')

cbar = fig.colorbar(sc)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Generations', rotation=270)


# had to do 2 seperate loops for knee points to appear on the foreground 
for gen in range(ngen):
    df = feas_front_df.loc[feas_front_df.generation==gen,:]
    # sort pareto 
    sdf = df.sort_values(objs[0])
    # identify knee/elbow pareto optimum 
    kn = kneed.KneeLocator(
        sdf.loc[:,objs[0]],
        sdf.loc[:,objs[1]],
        curve='concave',
        direction='increasing',
        interp_method='interp1d',
    )
    if  kn.knee is None :
        print(f'Could not identify knee point for gen {gen}')
        continue
    dmds = ax.scatter(kn.knee,kn.knee_y,color='blue',
               marker="D",s=80,label='knee point')
    ax.annotate(gen,(kn.knee,kn.knee_y),ha='center',va='center',
                color='white',fontsize=8,)

# fac1 configuration for comparative purpose 
pfac1 = ax.scatter(fac1_deficit,fac1_pump,marker="+", edgecolor='black',color='black',s=60)

ax.legend([dmds,pfac1],['knee points','Factor=1'])
ax.set_xlabel('Total deficit [m$^3$]')
ax.set_ylabel('Total pumping [m$^3$]')
ax.set_xlim(feas_front_df.loc[:,objs[0]].min(),
                      feas_front_df.loc[:,objs[0]].max())
ax.set_ylim(feas_front_df.loc[:,objs[1]].min(),
                      feas_front_df.loc[:,objs[1]].max())


fig.savefig(os.path.join('pproc','allgens_pareto.pdf'), dpi=300)

# ----------------------------------------------------
# pproc mou  : last generation pareto with fac1 config
# ----------------------------------------------------

gen = ngen-1 

fig,ax = plt.subplots(1,1,figsize=(5,5))
objs = pst.pestpp_options["mou_objectives"].split(',')
df = feas_front_df.loc[feas_front_df.generation==gen,:]
ax.scatter(df.loc[:,objs[0]],df.loc[:,objs[1]],
           marker="o",s=60,color='darkred',label=f'Pareto for gen. {gen}',alpha=0.7)

# extract and plot member ids
'''
mids = df.member.str.extract(r'member=(\d*)_', expand=False)
for mid, x, y in zip(mids,df.loc[:,objs[0]], df.loc[:,objs[1]]):
    ax.annotate(mid,(x,y),ha='left',va='center',
                xytext=(3,0),textcoords='offset points',
                fontsize=6)
'''

# knee point 
sdf = df.sort_values(objs[0])
lastgenkn = kneed.KneeLocator(
        sdf.loc[:,objs[0]],
        sdf.loc[:,objs[1]],
        curve='concave',
        direction='increasing',
        interp_method='interp1d',
    )

#kn.plot_knee_normalized()

# knee
ax.scatter(lastgenkn.knee,lastgenkn.knee_y,marker='D',color='royalblue',s=60, label='Knee point')

# fac1 configuration for comparative purpose 
ax.scatter(fac1_deficit,fac1_pump,marker="+", edgecolor='black',color='black',s=60,label='Factor=1')

ax.legend(loc='lower left')

ax.set_xlabel('Total deficit [m$^3$]')
ax.set_ylabel('Total pumping [m$^3$]')
ax.set_xlim(feas_front_df.loc[:,objs[0]].min(),
                      feas_front_df.loc[:,objs[0]].max())
ax.set_ylim(feas_front_df.loc[:,objs[1]].min(),
                      feas_front_df.loc[:,objs[1]].max())

fig.savefig(os.path.join('pproc',f'lastgen_pareto.pdf'),dpi=300)

# ---------------------------------------------
# pproc mou : dv over generations (convergence check)
# ---------------------------------------------

dic_props = dict(color = 'k', lw =  0.5)
meanprops = dict(marker = 'o', mfc = 'none', ms = 3, lw = 0.5, mec = 'k')
flierprops = dict(marker = 'o', mfc = 'none', ms = 3, lw = 0.1, mec = 'k')

def plot_multviolins(ens_list, ax, showpoints=False):
    # boxplot
    bplot = ax.boxplot(ens_list, widths = 0.4, showbox = True, showcaps = False, 
                showmeans = True, showfliers = False, boxprops = dic_props, 
                whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
                meanprops = meanprops, flierprops = flierprops)
    # violin plot (white)
    vp = ax.violinplot(ens_list, widths = 0.8, points = 100, showmeans = False, 
                       showmedians = False, showextrema = False)
    for pc in vp['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('white')
        pc.set_alpha(1)
    # violin plot (color)
    vp = ax.violinplot(ens_list, widths = 0.8, points = 100, showmeans = False,
                       showmedians = False, showextrema = False)
    for i, pc in enumerate(vp['bodies']):
        pc.set_facecolor('blue')
        pc.set_edgecolor('blue')
        pc.set_alpha(0.4)
    if showpoints != False: # Points
        for i, tick in enumerate(xticks_list):
            y = ens_list[i]
            x = np.random.normal(tick, 0.04, size = len(y))
            ax.plot(x, y, 'k.', alpha = 0.1)
    ax.set_xticklabels([])
    ax.set_axisbelow(True)
    ax.yaxis.grid()

# function to read decision variable population
def read_dv_pop(csvfile):
    dv_pop = pd.read_csv(csvfile)
    dv_pop.set_index('real_name',drop=True,inplace=True)
    # drop adjustable model parameters, keep dvars only
    dvpars = dv_pop.columns[dv_pop.columns.str.contains('fac')]
    dv_pop = dv_pop.loc[:,dvpars]
    fac_type = [ cname.split('pump')[0] for cname in dv_pop.columns]
    fac_id = np.array([ cname.split('_')[1] for cname in dv_pop.columns]).astype(int)
    fac_istep = np.array([ cname.split('_')[2] for cname in dv_pop.columns]).astype(int)
    fac_date = mm.mldates[fac_istep]
    fac_year = fac_date.year
    fac_month = fac_date.month
    dv_pop.columns = pd.MultiIndex.from_frame(
            pd.DataFrame({'type':fac_type,'id':fac_id,
                          'year':fac_year,'month':fac_month
                          })
            )
    return(dv_pop.T)

# plot rivpump factors 
fig,axs=plt.subplots(3,1,figsize=(8,8))
for l,ax in zip([1,3,5],axs):
    ens_list = []
    for gen in range(ngen):
        dv_pop = read_dv_pop(f'mou_lizonne.{gen}.dv_pop.csv')
        ptype = 'aq'
        ens_list.append(dv_pop.loc[(ptype,l,slice(None))].values.ravel())
    plot_multviolins(ens_list, ax)
    ax.set_title(f'Pumping factors for layer {l}')
    ax.set_ylabel('Factor value [-]')

ax.set_xticklabels(range(ngen))
ax.set_xlabel('Generations')
fig.tight_layout()
fig.savefig(os.path.join('pproc','evol_aqpump.pdf'),dpi=300)


# plot river factors 
fig,ax=plt.subplots(1,1,figsize=(8,4))
ens_list = []
for gen in range(ngen):
    dv_pop = read_dv_pop(f'mou_lizonne.{gen}.dv_pop.csv')
    ptype = 'riv'
    ens_list.append(dv_pop.loc[(ptype,slice(None),slice(None))].values.ravel())

plot_multviolins(ens_list, ax)
ax.set_xticklabels(range(ngen))
ax.set_xlabel('Generations')
ax.set_ylabel('Factor value [-]')
ax.set_title(f'Pumping factors for river reaches')
fig.tight_layout()
fig.savefig(os.path.join('pproc','evol_rivpump.pdf'),dpi=300)

# ---------------------------------------------
# dv analysis for last generation (gen = ngen)
# ---------------------------------------------

# plot aqpump factors 
dv_pop = read_dv_pop(f'mou_lizonne.{gen}.dv_pop.csv')
opt_year = dv_pop.index.levels[2].max()
months = dv_pop.index.levels[3]
layers = dv_pop.xs('aq').index.get_level_values(0).unique()
fig, axs=plt.subplots(len(months),1,figsize=(4,8))
for m,ax in zip(months,axs):
    ens_list = [ dv_pop.loc[('aq',l,opt_year,m)] for l in layers ]
    plot_multviolins(ens_list, ax)
    ax.set_title(f'Pumping factors for month {m}')
    ax.set_ylabel('Factor value [-]')

ax.set_xticklabels(layers)
ax.set_xlabel('Layer')
fig.tight_layout()
fig.savefig(os.path.join('pproc',f'lastgen_aqpumpfac.pdf'),dpi=300)

# plot rivpump factors 
dv_pop = read_dv_pop(f'mou_lizonne.{gen}.dv_pop.csv')
opt_year = dv_pop.index.levels[2].max()
months = dv_pop.index.levels[3]
reach_ids = dv_pop.xs('riv').index.get_level_values(0).unique()
fig, axs=plt.subplots(len(months),1,figsize=(8,8))
for m,ax in zip(months,axs):
    ens_list = [ dv_pop.loc[('riv',r,opt_year,m)] for r in reach_ids ]
    plot_multviolins(ens_list, ax)
    ax.set_title(f'Pumping factors for month {m}')
    ax.set_ylabel('Factor value [-]')

ax.set_xticklabels(reach_ids)
ax.set_xlabel('Reach id')
fig.tight_layout()
fig.savefig(os.path.join('pproc',f'lastgen_rivpumpfac.pdf'),dpi=300)



# ---------------------------------------------
# dv of pareto aqpump and rivpump (gen = ngen)
# ---------------------------------------------

# violin plots 
dv_pop = read_dv_pop(f'mou_lizonne.archive.dv_pop.csv')
opt_year = dv_pop.index.levels[2].max()
dv_pop = dv_pop.loc[slice(None),slice(None),opt_year,slice(None),:]

mmonths = {6:'June',7:'July',8:'August',9:'September'}

# plot aqpump and rivpump factors 
fig,axs=plt.subplots(1,4,figsize=(8,4),sharey=True)
for m,ax in zip([6,7,8,9],axs):
    ens_list = []
    ptype = 'aq'
    for l in [1,3,5]:
        ens_list.append(dv_pop.loc[('aq',l,m)].values.ravel())
    ens_list.append(dv_pop.loc[('riv',slice(None),m)].values.ravel())
    ax = plot_multviolins(ens_list, ax)
    ax.set_xticklabels(['l2','l4','l6','riv'])
    ax.set_title(mmonths[m])

axs[0].set_ylabel('Factor value [-]')
fig.tight_layout()

fig.savefig(os.path.join('pproc','evol_aqpump.pdf'),dpi=300)
#fig.savefig(os.path.join('evol_aqpump.pdf'),dpi=300)


# --- plot decorated rivpump factors map

# load shapefile of simulated river network
simriv_shp = os.path.join('gis','sim_riv.shp')
simriv_gdf = gpd.read_file(simriv_shp)
simriv_gdf.set_index(simriv_gdf.val.astype(int),inplace=True)

# load basin outline 
basin_shp = os.path.join('..','data','SIG','BassinLizonne.shp')
basin_gdf = gpd.read_file(basin_shp)

# load histo file and convert to gpd
histo_df = marthe_utils.read_histo_file(mm.mlname+'.histo')
histo_gdf = gpd.GeoDataFrame(histo_df,
                       geometry = gpd.points_from_xy(histo_df['x'], 
                                                     histo_df['y']),
                       crs = 2154) # EPSG:RGF93

gstations_gdf = histo_gdf.loc[histo_gdf['type']=='Débit_Rivi']

gstations_gdf['label2']= gstations_gdf.label.replace(['P8284010', 'P8215010', 'P7250001', 'P7270001'],
                                       ['GS1','GS4','GS3','GS2'])

# load factors 
rivfacvals = dv_pop.loc[('riv',slice(None),slice(None))]
mrivfacvals = rivfacvals.mean(axis=1).unstack()
simriv_gdf = simriv_gdf.merge(mrivfacvals,left_index=True, right_index=True)

# plot maps
fig,axs = plt.subplots(1,4,figsize=(10,3.5))
for m,ax in zip(mmonths.keys(),axs):
    # basin outline
    ax = basin_gdf.boundary.plot(ax=ax,color='black', linewidth=0.5)
    # pumping factors
    ax = simriv_gdf.plot(ax=ax,column=m,vmin=0,vmax=10)
    # gaging stations
    ax = gstations_gdf.plot(marker='^',c='red',
        edgecolor='black', markersize=22, label='Gaging stations', ax=ax)
    gstations_gdf.apply(lambda x: ax.annotate(
        text=x['label2'], xy=x.geometry.coords[0], 
        xytext=(3, -5), textcoords="offset points",
        fontsize=6), axis=1)
    ax.set_title(mmonths[m])
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

fig.tight_layout()

# append colorbar
fig.subplots_adjust(bottom=0.20)
p0 = axs[0].get_position().get_points().flatten()
p1 = axs[1].get_position().get_points().flatten()
p2 = axs[2].get_position().get_points().flatten()
p3 = axs[3].get_position().get_points().flatten()
cbar_ax = fig.add_axes([p1[0], 0.2, p2[2]-p1[0], 0.03]) # (left, bottom, width, height)
fig.colorbar(axm.get_children()[1],cax=cbar_ax, orientation='horizontal', label='factor')

fig.savefig(os.path.join('pproc',f'pareto_rivpumpfac_map.png'),dpi=300)

# ---------------------------------------------
# objective values at knee point  
# ---------------------------------------------

# read dv and obs pop archives
dv_pop = pd.read_csv('mou_lizonne.archive.dv_pop.csv',index_col=0)

obs_pop = pd.read_csv('mou_lizonne.archive.obs_pop.csv')
obs_pop.set_index(obs_pop.real_name.str.extract(r'member=(\d*)_', expand=False).astype(int),inplace=True)

# realization number of knee point 
realkn = obs_pop.index[(obs_pop.deficit_tot==lastgenkn.knee) & (obs_pop.tot_pump==lastgenkn.knee_y)][0].astype(int)

# --- compare objectives 

# member realkn
kn_pump = obs_pop.loc[realkn,'tot_pump']
kn_deficit = obs_pop.loc[realkn,'deficit_tot']

pump_df = pd.DataFrame({
    'fac1':fac1_rei.loc[fac1_rei.index.str.startswith('tot_'),'modelled'],
    'kn':obs_pop.loc[realkn,obs_pop.columns.str.startswith('tot_')],
    })


deficit_df = pd.DataFrame({
    'fac1':fac1_rei.loc[fac1_rei.index.str.startswith('deficit'),'modelled'],
    'kn':obs_pop.loc[realkn,obs_pop.columns.str.startswith('deficit')],
    })


# set more explicit names
pump_df.index = pump_df.index.map({'tot_aqpump':'aq', 'tot_pump':'tot', 'tot_rivpump':'riv'})
pump_df.sort_index(inplace=True)
deficit_df.index = pd.Series(deficit_df.index).apply(lambda x: x.split('_')[1])

# plot 
color_dic = {'fac1':'tan','kn':'seagreen','m4975':'royalblue'}
fig,axs=plt.subplots(1,2,figsize=(9,5))
axl = pump_df.plot.bar(ax=axs[0],legend=True,
                       color=color_dic)
axl.set_title('Pumping')
axl.set_xlabel('')
axl.set_ylabel('m$^3$')
axl.set_xticks(axl.get_xticks(), axl.get_xticklabels(), rotation=45, ha='right')

axr = deficit_df.plot.bar(ax=axs[1],legend=False,
                          color=color_dic)
axr.set_title('River deficit')
axr.set_xlabel('')
axr.set_xticks(axr.get_xticks(), axr.get_xticklabels(), rotation=45, ha='right')
axr.set_ylabel('m$^3$')
fig.tight_layout()

fig.savefig(os.path.join('pproc','kneepoint_objvals.pdf'),dpi=300)

# ---------------------------------------------
# plot pumping factors (dv) for real realkn
# ---------------------------------------------

# read dv and obs pop archives
dv_pop = pd.read_csv('mou_lizonne.archive.dv_pop.csv')
dv_pop.set_index(dv_pop.real_name.str.extract(r'member=(\d*)_', expand=False).astype(int),inplace=True)
dv_pop = dv_pop.loc[:,dv_pop.columns.str.contains('fac') ]
# select realization corresponding to kneepoint 
real = realkn

# --- plot aqpump factors  
par = dv_pop.iloc[:,1:].loc[real]
fac_type = [ pname.split('pump')[0] for pname in par.index]
fac_id = np.array([ pname.split('_')[1] for pname in par.index]).astype(int)
fac_istep = np.array([ pname.split('_')[2] for pname in par.index]).astype(int)
fac_date = mm.mldates[fac_istep]
fac_year = fac_date.year
fac_month = fac_date.month
opt_year = max(fac_year)
par.index = pd.MultiIndex.from_frame(
        pd.DataFrame({'type':fac_type,'id':fac_id,
                      'year':fac_year,'month':fac_month
                      })
        )
aqfacvals = par.loc[('aq',slice(None),opt_year)]
fig,ax = plt.subplots(1,1,figsize=(5,5))
aqfacvals.unstack().T.plot(ax=ax,kind='bar',legend=True)
ax.set_ylabel('Factor value [-]')
ax.axhline(y=1,color='grey', ls='--')
fig.tight_layout()
fig.savefig(os.path.join('pproc',f'kneepoint_aqpumpfac_bars.pdf'),dpi=300)

# --- plot aqpump factors  
par = dv_pop.iloc[:,1:].loc[real]
fac_type = [ pname.split('pump')[0] for pname in par.index]
fac_id = np.array([ pname.split('_')[1] for pname in par.index]).astype(int)
fac_istep = np.array([ pname.split('_')[2] for pname in par.index]).astype(int)
fac_date = mm.mldates[fac_istep]
fac_year = fac_date.year
opt_year = max(fac_year)
fac_month = fac_date.month
par.index = pd.MultiIndex.from_frame(
        pd.DataFrame({'type':fac_type,'id':fac_id,
                      'year':fac_year,'month':fac_month
                      })
        )
rivfacvals = par.loc[('riv',slice(None),opt_year)]
vmin,vmax= rivfacvals.min(),rivfacvals.max()
# aggregate over years 
fig,ax = plt.subplots(1,1,figsize=(4,4))
pd.DataFrame(rivfacvals).boxplot(ax=ax,by='month',whis=(0,100),showmeans=True)
ax.set_ylabel('Average factor value [-]')
fig.tight_layout()
fig.savefig(os.path.join('pproc',f'kneepoint_rivpumpfac_bars.pdf'),dpi=300)


# --- plot rivpump factors map

# load shapefile of simulated river network
simriv_shp = os.path.join('gis','sim_riv.shp')
simriv_gdf = gpd.read_file(simriv_shp)
simriv_gdf.set_index(simriv_gdf.val.astype(int),inplace=True)

#def barplot_rivfac(dv_pop,real,simriv_gdf):
par = dv_pop.iloc[:,1:].loc[real]
fac_type = [ pname.split('pump')[0] for pname in par.index]
fac_id = np.array([ pname.split('_')[1] for pname in par.index]).astype(int)
fac_istep = np.array([ pname.split('_')[2] for pname in par.index]).astype(int)
fac_date = mm.mldates[fac_istep]
fac_year = fac_date.year
opt_year = max(fac_year)
fac_month = fac_date.month
par.index = pd.MultiIndex.from_frame(
    pd.DataFrame({'type':fac_type,'id':fac_id,
                  'date':fac_date.strftime('%Y-%m'),
                  'year':fac_year
                  })
    )
rivfacvals = par.loc[('riv',slice(None),slice(None),opt_year)]
simriv_gdf = simriv_gdf.merge(rivfacvals.unstack(),left_index=True, right_index=True)
plot_dates = rivfacvals.index.get_level_values('date').unique()
fig,axs = plt.subplots(1,4,figsize=(10,4))
for d,ax in zip(plot_dates,axs.ravel()):
    simriv_gdf.plot(ax=ax,column=d,vmin=0,vmax=2)
    ax.set_title(d)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

fig.subplots_adjust(bottom=0.2)
p0 = axs[0].get_position().get_points().flatten()
p1 = axs[1].get_position().get_points().flatten()
p2 = axs[2].get_position().get_points().flatten()
p3 = axs[3].get_position().get_points().flatten()
# add ax rect : (left, bottom, width, height)
cbar_ax = fig.add_axes([p1[0], 0.2, p2[2]-p1[0], 0.05])
fig.colorbar(ax.get_children()[0],cax=cbar_ax, orientation='horizontal', label='factor')

fig.savefig(os.path.join('pproc',f'kneepoint_rivpumpfac_map.pdf'),dpi=300)

#====================================================
# plot aquifer constraints
#====================================================

# ---  uncertainties on gw levels in pumped cells derived from par stack at 1st iteration

# load hist file to get aqpump data
histo_df = marthe_utils.read_histo_file('Lizonne.histo')
aqpumpcells = histo_df.loc[histo_df.index.str.startswith('aqpump')]
aqpumpcells.index= aqpumpcells.label.apply(lambda x: x.split('_')[1]).astype(int)
aqpumpcells.index.name='node'

# load obs stack
obs_pop = pd.read_csv('mou_lizonne.0.obs_stack.csv')
h_pop = obs_pop.loc[:,obs_pop.columns.str.startswith('h_')].copy()
h_pop.index.name='real'
# re-index with nodes and time steps 
nodes, tsteps=h_pop.columns.str.extract(r'h_(\d*)_n(\d*)').T.values
h_pop.columns=pd.MultiIndex.from_frame(pd.DataFrame(
    {'tstep':tsteps.astype(int),'node':nodes.astype(int)}))

# multiindex df with time series
h = h_pop.T.unstack()

# std of min h over obs stack 
hmin = h.loc[520:].min() # here, the min is on the time series
hstd = hmin.groupby(level='node').std()
hstd = pd.DataFrame({'std':hstd}).merge(aqpumpcells[['x','y','layer']],left_index=True,right_index=True)

axs = hstd.hist(column='std',by='layer')
axs[0,0].figure.savefig(os.path.join('pproc','uncert_hist.pdf'),dpi=300)

# convert 
hstd_gdf = gpd.GeoDataFrame(hstd,
                       geometry = gpd.points_from_xy(hstd.x, 
                                                     hstd.y),
                       crs = 2154)

# load basin outline 
basin_shp = os.path.join('..','data','SIG','BassinLizonne.shp')
basin_gdf = gpd.read_file(basin_shp)

# plot map
fig,axs=plt.subplots(1,3,figsize=(8,5))
vmin,vmax= 0,round(hstd['std'].max(),-1) # round
for i,l in enumerate([2,4,6]):
    ax=hstd_gdf.loc[hstd_gdf.layer==l].plot(marker='o',
                                     column='std',
                                       vmin=vmin,
                                       vmax=vmax,
                                     edgecolor='black',
                                     ax=axs[i])
    #ax.legend(title=f'Layer {l}',frameon=False,fontsize=7)
    ax = basin_gdf.boundary.plot(ax=ax,color='black')    
    ax.set_title(f'Layer {l}')

fig.tight_layout()
cax = fig.add_axes([0.25,0.2,0.4,0.03])
fig.colorbar(axs[0].get_children()[0],cax,orientation='horizontal', label='std of min sim. head [m]')

fig.savefig(os.path.join('pproc','uncert_maps.pdf'),dpi=300)


# --- plot series of gw levels (heads) at pumped aquifer cells (where constraints apply) 

# load constraints on h
aqpumplim_df = pd.read_csv('aqpump_lim.csv',index_col=0)
hmin = aqpumplim_df['hmin'].copy() # here, the min is the constraint value
hmin.index = [int(x.split('_')[1]) for x in hmin.index]
hmin.index.name = 'node'

def plot_constraints_tseries(obs_pop_file, gen):
    obs_pop = pd.read_csv(obs_pop_file) 
    h_pop = obs_pop.loc[:,obs_pop.columns.str.startswith('h_')].copy()
    h_pop.index.name='real'

    # re-index with nodes and time steps 
    nodes, tsteps=h_pop.columns.str.extract(r'h_(\d*)_n(\d*)').T.values
    h_pop.columns=pd.MultiIndex.from_frame(pd.DataFrame(
        {'tstep':tsteps.astype(int),'node':nodes.astype(int)}))

    # multiindex df with time series
    h = h_pop.T.unstack()

    # distance to hmin (neg if constraint is not satisfied)
    dh = h.sub(hmin,level='node')

    # identify nodes were constraint is not satisfied from tstep=520 (after initialization)
    n_nonfeas_reals = (dh.loc[dh.index>520].min(axis=0) < 0).groupby('node').sum()
    infeasible_cells = n_nonfeas_reals.loc[n_nonfeas_reals>0].index

    print(f'WARNING: {infeasible_cells} infeasible cells were found at generation {gen}')

    # plot on multi-page pdf 
    from matplotlib.backends.backend_pdf import PdfPages
    filename = os.path.join('pproc',f'constraints_{gen}.pdf')
    figsize=(8, 10.5)
    nr, nc = 4, 2

    # list of pumped aquifer nodes (cell id)
    nodes = h.columns.get_level_values(1).unique()

    figs = []
    ax_count = 0
    print(f'Generating figures for generation {gen}...')
    for n in nodes:
        # new figure (pdf page)
        if ax_count % (nr * nc) == 0:
            ax_count = 0
            fig, ax_mat = plt.subplots(nr, nc,figsize=figsize)
            axs=ax_mat.ravel()
        # new ax (plot)
        ax = axs.ravel()[ax_count]
        ax = h.loc[:,(slice(None),n)].plot(ax=ax,legend=False)
        l=aqpumpcells.loc[n,'layer']
        ax.set_title(f'cell {n} - layer {l}')
        ax.axhline(hmin.loc[n],c='grey',ls='--')
        ax.axvline(520,c='grey',ls='--')
        ax.set_xlabel('')
        ax_count += 1
        # save fig
        if ax_count == nr*nc-1 :
            fig.tight_layout()
            figs.append(fig)

    print('Writing pdf...')
    with PdfPages(filename) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)


obs_pop_file = 'mou_lizonne.0.obs_pop.chance.csv'
plot_constraints_tseries(obs_pop_file,0)
obs_pop_file = f'mou_lizonne.{ngen-1}.chance.obs_pop.csv'
plot_constraints_tseries(obs_pop_file,ngen-1)

# run pproc_opt.py
#====================================================
#=> plot aquifer constraints for single real
#====================================================

# run pproc_opt.py

