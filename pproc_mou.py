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
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib 
import kneed 
# plot settings
plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54
mdcol_width = 14/2.54
dbcol_width = 19/2.54
color_dic = {'fac1':'tan','opt':'royalblue'}

# clean and (re)-build directory tree
shutil.rmtree('pproc')
os.mkdir('pproc')
for d in ['haqpump','qriv','heads']:
    os.mkdir(os.path.join('pproc',d))


# analyze dv for knee point
mm = MartheModel('Lizonne.rma', spatial_index = True)

# load historical reference (fac=1)
sim_dir = os.path.split(os.getcwd())[-1].split('master_')[1] 
fac1_rei = pyemu.pst_utils.read_resfile(os.path.join('..', sim_dir,'mou_lizonne_fac1.base.rei'))
fac1_pump = fac1_rei.loc['tot_pump','modelled']
fac1_deficit = fac1_rei.loc['deficit_tot','modelled']


# generation selected for analysis
opt_gen = 10

# ---------------------------------------------
# pproc mou  : pareto
# ---------------------------------------------

pst = pyemu.Pst('mou_lizonne.pst')

# summary of pareto dominant solutions for each generation
pasum_df = pd.read_csv('mou_lizonne.pareto.archive.summary.csv')
feas_front_df = pasum_df.loc[pasum_df.apply(lambda x: x.nsga2_front==1 and x.is_feasible==1,axis=1),:]
feas_front_members = feas_front_df.loc[feas_front_df.generation==opt_gen,'member'].values
ngen = feas_front_df.generation.unique().shape[0]

cmap = matplotlib.colormaps.get_cmap('gist_heat').reversed()

fig,ax = plt.subplots(1,1,figsize=(sgcol_width,sgcol_width))
objs = pst.pestpp_options["mou_objectives"].split(',')
for gen in range(ngen):
    df = feas_front_df.loc[feas_front_df.generation==gen,:]
    sc = ax.scatter(df.loc[:,objs[0]]*1e-6,df.loc[:,objs[1]]*1e-6,c=df.loc[:,'generation'],vmin=0,vmax=ngen,
               cmap=cmap,marker="o", 
               label=f'gen. {gen}')

cbar = fig.colorbar(sc, label='Generations',orientation='horizontal',
             cax=ax.inset_axes((0.4, 0.20, 0.5, 0.05)))

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
    dmds = ax.scatter(kn.knee*1e-6,kn.knee_y*1e-6,color='blue',
               marker="D",s=80,label='knee point')
    ax.annotate(gen,(kn.knee*1e-6,kn.knee_y*1e-6),ha='center',va='center',
                color='white',fontsize=8,)

# fac1 configuration for comparative purpose 
pfac1 = ax.scatter(fac1_deficit*1e-6,fac1_pump*1e-6,marker="+", edgecolor='black',color='black',s=60)

ax.legend([dmds,pfac1],['knee points','Factor=1'])
ax.set_xlabel('Total deficit [Mm$^3$]')
ax.set_ylabel('Total pumping [Mm$^3$]')
'''
ax.set_xlim(feas_front_df.loc[:,objs[0]].min(),
                      feas_front_df.loc[:,objs[0]].max())
ax.set_ylim(feas_front_df.loc[:,objs[1]].min(),
                      feas_front_df.loc[:,objs[1]].max())
ax.set_ylim([0,32])
ax.set_xlim([0,9])
'''
fig.tight_layout()
fig.savefig(os.path.join('pproc','allgens_pareto.pdf'), dpi=300)

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
    return(ax)

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
    return(dv_pop)

# plot rivpump factors 
fig,axs=plt.subplots(3,1,figsize=(8,8))
for l,ax in zip([1,3,5],axs):
    ens_list = []
    for gen in range(ngen):
        dv_pop = read_dv_pop(f'mou_lizonne.{gen}.dv_pop.csv')
        ptype = 'aq'
        ens_list.append(dv_pop.T.loc[(ptype,l,slice(None))].values.ravel())
    plot_multviolins(ens_list, ax)
    ax.set_title(f'Pumping factors for layer {l}')
    ax.set_ylabel('Factor value [-]')

ax.set_xticklabels(range(ngen))
ax.set_xlabel('Generations')
fig.tight_layout()
fig.savefig(os.path.join('pproc',f'aqpumpfac_gen{gen}.pdf'),dpi=300)


# plot river factors 
fig,ax=plt.subplots(1,1,figsize=(8,4))
ens_list = []
for gen in range(ngen):
    dv_pop = read_dv_pop(f'mou_lizonne.{gen}.dv_pop.csv')
    ptype = 'riv'
    ens_list.append(dv_pop.T.loc[(ptype,slice(None),slice(None))].values.ravel())

plot_multviolins(ens_list, ax)
ax.set_xticklabels(range(ngen))
ax.set_xlabel('Generations')
ax.set_ylabel('Factor value [-]')
ax.set_title(f'Pumping factors for river reaches')
fig.tight_layout()
fig.savefig(os.path.join('pproc',f'rivpumpfac_gen{gen}.pdf'),dpi=300)

# ---------------------------------------------
# absolute pumping values
# ---------------------------------------------

# --- load original aquifer pumping 
org_parfile = os.path.join('pest','par','aqpump.dat.org')
keys=['boundname','layer','istep']
parkmi, parvals = pest_utils.parse_mlp_parfile(org_parfile, keys=keys, value_col=1, btrans='none')
parvals.index=parkmi
aqpump_org = parvals.groupby(level=('layer','istep')).sum()
dates = mm.mldates[aqpump_org.index.get_level_values('istep')]
aqpump_org.index=pd.MultiIndex.from_frame(pd.DataFrame({
            'type':'aq',
            'id':aqpump_org.index.get_level_values('layer'),
            'year':dates.year,
            'month':dates.month
            }))

# --- load original river pumping 

org_parfile = os.path.join('pest','par','rivpump.dat.org')
keys=['prefix','aff_r','istep']
parkmi, parvals = pest_utils.parse_mlp_parfile(org_parfile, keys=keys, value_col=1, btrans='none')
parvals.index=parkmi
rivpump_org = parvals.groupby(level=('aff_r','istep')).sum()
dates = mm.mldates[rivpump_org.index.get_level_values('istep')]
rivpump_org.index=pd.MultiIndex.from_frame(pd.DataFrame({
            'type':'riv',
            'id':rivpump_org.index.get_level_values('aff_r'),
            'year':dates.year,
            'month':dates.month
            }))

# concat and convert to m3/month
monthly_m3sec_to_m3 = -1*(365.25/12)*86400 # (days in a month) * (seconds in a day)
pump_org = pd.concat([aqpump_org,rivpump_org])*monthly_m3sec_to_m3

# ---------------------------------------------
# dv analysis for opt_gen
# ---------------------------------------------



# plot aqpump factors 

# load dv and obs of pareto members at gen=opt_gen
dv_pop = read_dv_pop(f'mou_lizonne.{opt_gen}.dv_pop.csv')
dv_pop = dv_pop.loc[feas_front_members] # subset to pareto members
obs_pop = pd.read_csv(f'mou_lizonne.{opt_gen}.chance.obs_pop.csv',index_col=0)
obs_pop = obs_pop.loc[feas_front_members] # subset to pareto members

facvals = dv_pop.T
pump_opt = facvals.mul(pump_org.loc[facvals.index],axis=0)
mmonths = {6:'June',7:'July',8:'August',9:'September'}


# plot
fig,axs=plt.subplots(2,4,figsize=(dbcol_width,0.6*dbcol_width))
for m,ax in zip([6,7,8,9],axs[0,:]):
    year = dv_pop.columns.get_level_values('year').max()
    # plot riv factors (all reaches)
    ax.bar(0,pump_org.loc[('riv',slice(None),year,m)].sum(),color='grey',alpha=0.5)
    # reduce value of river pumping for plotting 
    pump_opt_vals =[min(p,4e6) for p in  pump_opt.loc[('riv',slice(None),year,m)].sum(axis=0).values]
    ax.plot([0]*len(pump_opt_vals),pump_opt_vals,ls='',marker='+', c='k')
    # plot aq factors 
    for i,l in enumerate([1,3,5]):
        ax.bar(i+1,pump_org.loc[('aq',l,year,m)],color='grey',alpha=0.5,label='Original')
        pump_opt_vals = pump_opt.loc[('aq',l,year,m),:].values
        ax.plot([i+1]*len(pump_opt_vals),pump_opt_vals,ls='',marker='+', c='k',label= 'Optimized')
    # rename ticks 
    ax.xaxis.set_ticks(range(i+2))
    if m>6:
        ax.set_yticklabels([])
    ax.set_xticklabels(['RIV','COST','TURO','CENO'],fontsize=8)
    ax.set_title(mmonths[m])
    ax.set_ylim(0,4e6)


ytlbls = axs[0,0].get_yticklabels()
ytlbls[-1].set_text('>4')
axs[0,0].set_yticklabels(ytlbls)
handles, labels = axs[0,0].get_legend_handles_labels()
axs[0,0].legend(handles=[handles[0],handles[-1]],labels=[labels[0],labels[-1]],loc='upper right')
axs[0,0].set_ylabel('Pumping [Mm$^3$/month]')

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
rivfacvals = dv_pop.T.loc[('riv',slice(None),slice(None))]
mrivfacvals = rivfacvals.mean(axis=1).unstack()
mrivfacvals.index = mrivfacvals.index.droplevel('year')
simriv_gdf = simriv_gdf.merge(mrivfacvals,left_index=True, right_index=True)

# plot maps
for m,ax in zip(mmonths.keys(),axs[1,:]):
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
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('X (eastward)')
    ax.add_artist(ScaleBar(1))

axs[1,0].set_ylabel('Y (northward)')
fig.tight_layout()

# append colorbar
fig.subplots_adjust(bottom=0.20)
p0 = axs[1,0].get_position().get_points().flatten()
p1 = axs[1,1].get_position().get_points().flatten()
p2 = axs[1,2].get_position().get_points().flatten()
p3 = axs[1,3].get_position().get_points().flatten()
cbar_ax = fig.add_axes([p1[0], 0.1, p2[2]-p1[0], 0.02]) # (left, bottom, width, height)
fig.colorbar(ax.get_children()[1],cax=cbar_ax, orientation='horizontal', label='River Pumping Factors [-]')


# set the spacing between subplots
fig.subplots_adjust(left=0.08)
fig.subplots_adjust(wspace=0.12)
fig.subplots_adjust(hspace=0.30)
fig.subplots_adjust(right=0.975)
fig.subplots_adjust(bottom=0.20)
fig.subplots_adjust(top=0.925)

# 
line = matplotlib.lines.Line2D((0.025,0.975),(0.55,0.55),color='k',linewidth=1,transform=fig.transFigure)
fig.lines.append(line)

fig.text(0.008,0.95,'(a)',fontsize=14,transform=fig.transFigure)
fig.text(0.008,0.49,'(b)',fontsize=14,transform=fig.transFigure)

fig.savefig(os.path.join('pproc',f'pareto_rivpumpfac.pdf'),dpi=300)


# ---------------------------------------------
# objective values at knee point  
# ---------------------------------------------

# load dv and obs of pareto members at gen=opt_gen
dv_pop = read_dv_pop(f'mou_lizonne.{opt_gen}.dv_pop.csv')
dv_pop = dv_pop.loc[feas_front_members] # subset to pareto members
obs_pop = pd.read_csv(f'mou_lizonne.{opt_gen}.chance.obs_pop.csv',index_col=0)

# knee point of last gen
df = feas_front_df.loc[feas_front_df.generation==opt_gen,:]
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

# realization number of knee point 
realkn = obs_pop.index[(obs_pop.deficit_tot==kn.knee) & (obs_pop.tot_pump==kn.knee_y)][0]

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

# select realization corresponding to kneepoint 
real = realkn

# read dv and obs pop archives
dv_pop = pd.read_csv('mou_lizonne.archive.dv_pop.csv',index_col=0)
par = dv_pop.loc[realkn,dv_pop.columns.str.contains('fac') ].T

# --- plot aqpump factors  
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

simriv_gdf = simriv_gdf.merge(rivfacvals.unstack(),left_index=True, right_index=True)
months = rivfacvals.index.get_level_values('month').unique()
fig,axs = plt.subplots(1,4,figsize=(10,4))
for m,ax in zip(months,axs.ravel()):
    simriv_gdf.plot(ax=ax,column=m,vmin=0,vmax=2)
    ax.set_title(m)
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


# load hist file to get aqpump data
histo_df = marthe_utils.read_histo_file('Lizonne.histo')
aqpumpcells = histo_df.loc[histo_df.index.str.startswith('aqpump')]
aqpumpcells.index= aqpumpcells.label.apply(lambda x: x.split('_')[1]).astype(int)
aqpumpcells.index.name='node'

# load obs stack
obs_pop = pd.read_csv('mou_lizonne.0.obs_stack.csv',index_col='real_name')


# get heads and re-index with nodes and time steps 
h_pop = obs_pop.loc[:,obs_pop.columns.str.startswith('h_')].copy()
nodes, tsteps=h_pop.columns.str.extract(r'h_(\d*)_n(\d*)').T.values
h_pop.columns=pd.MultiIndex.from_frame(pd.DataFrame(
    {'tstep':tsteps.astype(int),'node':nodes.astype(int)}))
h = h_pop.T.unstack() # multiindex df with time series

# std of min h over obs stack 
hmin = h.loc[520:].min() # min of time series after warm up 
hstd = hmin.groupby(level='node').std()
hstd = pd.DataFrame({'std':hstd}).merge(aqpumpcells[['x','y','layer']],left_index=True,right_index=True)


# convert to geopandas 
hstd_gdf = gpd.GeoDataFrame(hstd,
                       geometry = gpd.points_from_xy(hstd.x,
                                                     hstd.y),
                       crs = 2154)

# load basin outline 
basin_shp = os.path.join('..','data','SIG','BassinLizonne.shp')
basin_gdf = gpd.read_file(basin_shp)

# plot map
layer_dic = {2:'COST',4:'TURO',6:'CENO'}
#layer_dic = {1:'COST',3:'TURO',5:'CENO'}


# plot
fig,axs=plt.subplots(2,2,figsize=(0.80*dbcol_width,0.60*dbcol_width))

# uncertainties on total deficit
deficit_tot = obs_pop['deficit_tot'].div(1e6) # Mm3
stack_size = deficit_tot.shape[0]
ecdf = pd.Series(np.arange(1,stack_size+1)/stack_size,index=deficit_tot.sort_values().values)
ax = axs[0,0]
ax = deficit_tot.hist(ax=ax,color='grey',grid=False)
ax.set_ylabel('Frequency (counts)')
ax.set_xlabel('Mm$^3$')
ax.set_title('(a) Total River Deficit')
ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
twax = ax.twinx()
twax.plot(ecdf)
twax.axvline(deficit_tot.loc['ffpt_center'],ls='--',color='green',label='Center')
twax.plot(ecdf,color='red',label='CDF')
twax.set_ylabel('Cumulative density')
twax.set_ylim(0,1)
twax.legend(loc='upper right',bbox_to_anchor=(0.99,0.7))
ax.set_box_aspect(0.85)

# uncertainties map for heads 
vmin,vmax= 0,15 #to improve readability
for l,ax,label in zip([2,4,6],axs.ravel()[1:],['b','c','d']):
    ax=hstd_gdf.loc[hstd_gdf.layer==l].plot(marker='o',
                                     column='std',
                                       vmin=vmin,
                                       vmax=vmax,
                                       edgecolor='k',
                                     ax=ax)
    ax = basin_gdf.boundary.plot(ax=ax,color='black',linewidth=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('X (eastward)')
    ax.set_ylabel('Y (northward)')
    ax.set_title(f'({label}) {layer_dic[l]} Aquifer')
    ax.add_artist(ScaleBar(1))

#fig.tight_layout()
cax = fig.add_axes([0.85,0.2,0.02,0.6]) # left, bot, width, height
cmap = fig.colorbar(axs[1,1].get_children()[0],cax,orientation='vertical', label='Standard deviation of minimum head [m]')

cax.set_yticks(np.arange(vmin,vmax+1,5))
ytlbls = cax.get_yticklabels()
ytlbls[-1].set_text('>15')
cax.set_yticklabels(ytlbls)


# set the spacing between subplots
fig.subplots_adjust(left=0.025)
fig.subplots_adjust(wspace=0.)
fig.subplots_adjust(hspace=0.50)
fig.subplots_adjust(right=0.90)
fig.subplots_adjust(bottom=0.10)
fig.subplots_adjust(top=0.9)

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

