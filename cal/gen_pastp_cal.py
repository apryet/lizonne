import os, sys
import pandas as pd
import numpy as np
import datetime 

from pymarthe import MartheModel
from pymarthe.moptim import MartheOptim
from pymarthe.mfield import MartheField
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
import matplotlib.pyplot as plt
import pyemu
import geopandas as gpd

# -----------------------------------------------------------
# -- Read model
mm = MartheModel('Lizonne.rma',spatial_index=True)
#mm = MartheModel('Lizonne.rma',spatial_index='Lizonne_si')
# -----------------------------------------------------------

# export to shp
# NOTE : bug, handle EPSG code 
for l in range(mm.nlay):
    mm.prop['permh'].to_shapefile(
            os.path.join('gis',f'permh_{l+1:02d}.shp'),
            l,
            epsg=2154)
'''
# NOTE : bug, inconsistent map 
mm.get_outcrop().plot(cmap='tab20', masked_values=[-9999])
plt.show()
'''

# load basin outline 
basin_shp = os.path.join('..','data','SIG','BassinLizonne.shp')
basin_gdf = gpd.read_file(basin_shp)

# -----------------------------------------------------------
# load pumping data  
# -----------------------------------------------------------

# Read pumping data from excel file 
pump_data_xls = os.path.join('..','data','prelevements','Prelevements_V3-Lizonne.xlsx')

# selection of years (column names in excel file) for both groundwater and surface withdrawals
ycols = [y for y in range(2003,2020)]

# ===========================================================
# -------- processing groundwater withdrawals ---------------
# ===========================================================

# wells (aquifer)
aqpump_df = pd.read_excel(pump_data_xls,'Eaux-souterraines')

# remove wells without coordinates or layer id
aqpump_df.dropna(how='any',subset=['XLIIe','YLIIe','Code entite'],inplace=True)

# rename Type for use
aqpump_df['use'] = aqpump_df.Usage

# -----------------------------------------------------------
# convert to geopandas 
# -----------------------------------------------------------

aq_gdf = gpd.GeoDataFrame(aqpump_df,
                       geometry = gpd.points_from_xy(aqpump_df['XLIIe'], 
                                                     aqpump_df['YLIIe']),
                       crs = 27572) # EPSG:L2e

# project L2E to RGF93 and get corresponding coordinates 
aq_gdf = aq_gdf.to_crs(epsg=2154) # EPSG:RGF93
aq_gdf['x'] = aq_gdf.geometry.x
aq_gdf['y'] = aq_gdf.geometry.y

# -- assign model layer from 'code entité'

# layer ids are 0-based, layer=-1 when out of the model domain
layer_dic = {
        'Campanien': -1, # OC: le Campanien étant absent je ne prendrai pas en compte ces prélèvements
        'Coniacien - Santonien': 1,
        'Coniacien + Turonien': 3, # OC [...] OK pour l’affecter au Turonien qui est a priori le plus productif
        'Turonien + Cénomanien': 3, #idem
        'Coniacien + Turonien + Cénomanien': 3, # idem
        'Turonien': 3,
        'Cénomanien': 5,
        'Jurassique': -1 # OC : ne pas reporter les prélèvements jurassique sur d’autres nappes
        }

# map layer number (0-based) from ids 
aq_gdf['layer'] = aq_gdf['Code entite'].map(layer_dic)

# drop rows with layer = -1 
aq_gdf.drop(aq_gdf[aq_gdf.layer==-1].index,inplace=True)

# get model node from point layer and coordinates 
aq_gdf['node'] = mm.get_node(aq_gdf['x'],aq_gdf['y'],aq_gdf['layer'])

# --------------------------------------------
# -- plot map of groundwater wells
# --------------------------------------------
# layers ids are 0-based
geol_marker_dic = {
        0 : '.',
        1 : 'o', #Coniacien-Santonien (COST)
        2 : '.',
        3 :'v', #Turonien (TURO)
        4 :'.',
        5 :'*' #Cénomanien (CENO)
        }

geol_color_dic = {
        0 : 'lightgrey', #Eponte au sommet du Coniacien-Santonien (EPCS)
        1 :'khaki', #Coniacien-Santonien (COST)
        2 : 'lightgrey', #Eponte au sommet du Turonien (EPTU)
        3 : 'greenyellow',#Turonien (TURO)
        4 : 'lightgrey', #Eponte au sommet du Cenomanien (EPCE)
        5 : 'darkgreen'#Cénomanien (CENO)
        }

geol_id_dic = {
        0 :'EPCS',
        1 : 'COST',
        2 : 'EPTU',
        3 : 'TURO',
        4 : 'EPCE',
        5 : 'CENO'
        }



fig, ax = plt.subplots(figsize=(5,4))

# plot layers 
for l in reversed(range(mm.nlay)):
    permh_shp = os.path.join('gis',f'permh_{l+1:02d}.shp')
    permh_gdf = gpd.read_file(permh_shp)
    permh_gdf.plot(color=geol_color_dic[l],ax=ax,label=geol_id_dic[l],alpha=0.9)

ax.legend(title='Couche')

# plot basin
ax = basin_gdf.boundary.plot(ax=ax,color='black')


# plot with distinct marker by layer
for l in sorted(aq_gdf.layer.unique()):
    aq_gdf.loc[aq_gdf.layer==l].plot(marker=geol_marker_dic[l],
                                     color=geol_color_dic[l],
                                     edgecolor='black',
                                     linewidths=0.5,
                                     markersize=14,
                                     label=geol_id_dic[l],
                                     ax=ax)

ax.legend(title='Aquifer',frameon=False,fontsize=7)
fig.savefig(os.path.join('figs','qaq.png'),dpi=300)
ax_map=ax

# --------------------------------------------
# -- well/cell match
# --------------------------------------------

'''
df = aq_gdf.loc[:,['NUM-LIZONNE','layer','node','use']]
dfs=df.sort_values(['node','layer'])
dfs.reset_index(inplace=True)
dfs.to_excel(os.path.join('..','data','prelevements','corr_forage_cellule_sort.xlsx'))

df.set_index(['layer','node'],inplace=True)
df.sort_index(inplace=True)
df.to_excel(os.path.join('..','data','prelevements','corr_forage_cellule_idx.xlsx'))
''' 
# --------------------------------------------
# -- temporal disaggregation by months
# --------------------------------------------

# qdf with single-indexed rows (date) and multi-indexed columns (node, layer, use)
qdf = pd.DataFrame(aq_gdf[['use','layer','node']+ycols]).set_index(['use','layer','node']).T
qdf.index= pd.to_datetime(qdf.index,format='%Y')

# generate monthly time sequence at month start (MS)
start, end = pd.to_datetime(['2003-01-01','2019-12-31'])
months = pd.date_range(start, end, freq='MS') 

# update index, fill monthly records with annual values
dti_uni = qdf.index.union(months).drop_duplicates()
dfm = qdf.reindex(dti_uni).interpolate('pad',limit_direction='forward')

# generate df with monthly ratios of annual pumping 
dfw = dfm.copy()
dfw.iloc[:,:]=1./12 # values for AEP and IND uses
dfw.iloc[:,dfw.columns.get_level_values('use')=='IRRIGATION']=0. 
dfw.loc[dfw.index.month==6,dfw.columns.get_level_values('use')=='IRRIGATION']=0.05 
dfw.loc[dfw.index.month==7,dfw.columns.get_level_values('use')=='IRRIGATION']=0.45
dfw.loc[dfw.index.month==8,dfw.columns.get_level_values('use')=='IRRIGATION']=0.45
dfw.loc[dfw.index.month==9,dfw.columns.get_level_values('use')=='IRRIGATION']=0.05

# term-by-term product to get the dataframe with monthly values (m3) 
dff = dfm*dfw

# --------------------------------------------
# -- spatial aggregation by model cell & layer
# --------------------------------------------

# ---- groupby by node + layer and sum 
dfg = dff.groupby(['node','layer'],axis=1).sum()


# --------------------------------------------
# -- write annual files from average of monthly values
# --------------------------------------------
ndays_in_a_month = 365.25/12
dfa = dfg.mean().to_frame('V')*-1./(ndays_in_a_month*86400) # m3/month to m3/s

# get indices for given node (back to 1-based !)
dfa['C'] = mm.query_grid(node=dfa.index.get_level_values('node').values, target='j').values+1
dfa['L'] = mm.query_grid(node=dfa.index.get_level_values('node').values, target='i').values+1
dfa['P'] = dfa.index.get_level_values('layer').values+1  
# output filename
qfilename = os.path.join('prelevements_sout_mensuels',f'Q_inter_annuel.txt')
# write to csv
dfa[['V','C','L','P']].to_csv(qfilename, sep='\t', header=False, index=False, float_format='%.6f')

# --------------------------------------------
# -- write monthly files
# --------------------------------------------

for m in months:
    # transpose to get current month df 
    # with the given month as column name and layer, node as row index
    dft = dfg.loc[dfg.index==m].T
    # drop na, get value, colonne,ligne,plan
    dft.dropna(subset=m,inplace=True)
    # number of days this month
    ndays = (m.replace(month = m.month % 12 +1, day = 1)-datetime.timedelta(days=1)).day
    # convert m3/month to m3/s 
    dft['V'] = dft[m]*-1./(ndays*86400)
    # get indices for given node (back to 1-based !)
    dft['C'] = mm.query_grid(node=dft.index.get_level_values('node').values, target='j').values+1
    dft['L'] = mm.query_grid(node=dft.index.get_level_values('node').values, target='i').values+1
    dft['P'] = dft.index.get_level_values('layer').values+1  # back to 1-based !
    # output filename
    qfilename = os.path.join('prelevements_sout_mensuels',f'Q_{m.year}_{m.month:02d}.txt')
    # write to csv
    dft[['V','C','L','P']].to_csv(qfilename, sep='\t', header=False, index=False, float_format='%.6f')


# --------------------------------------------
# -- plot monthly records of gw withdrawls
# --------------------------------------------

fig, ax = plt.subplots(figsize=(8,3))
import matplotlib.ticker as ticker
ax = dff.groupby(['use'],axis=1).sum()[['AEP','IRRIGATION']].plot.bar(
    stacked=True,
    color=['darkblue','darkorange','darkgreen'],
    ax=ax)
# Make most of the ticklabels empty so the labels don't get too crowded
ticklabels = ['']*len(dff.index)
# Every 12th ticklabel includes the year
ticklabels[::12] = [item.strftime('%m-%Y') for item in dff.index[::12]]
ax.xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
ax.set_ylabel('Prélèvements sout. [m$^3$/mois]')
plt.gcf().autofmt_xdate()
fig.savefig(os.path.join('figs','qsout_records.png'),dpi=300)

# copy gw withdrawal
dff_gw = dff

# ===========================================================
# ------ processing surface water withdrawals ---------------
# ===========================================================

# river withdrawals
riv_df = pd.read_excel(pump_data_xls,'Eaux-de-surface')

# remove pumping stations without coordinates 
riv_df.dropna(how='any', subset=['X_L93New','Y_L93New'], inplace=True)

# rename Usage for use
riv_df['use'] = riv_df.Usage

# drop rows with Infiltration
riv_df.drop(riv_df[riv_df['Milieur récepteur']=='Infiltration'].index,inplace=True)

# -----------------------------------------------------------
# convert to geopandas 
# -----------------------------------------------------------

riv_gdf = gpd.GeoDataFrame(riv_df,
                       geometry = gpd.points_from_xy(riv_df['X_L93New'], 
                                                     riv_df['Y_L93New']),
                       crs = 2154) # EPSG:RGF93

riv_gdf['x'] = riv_gdf.geometry.x
riv_gdf['y'] = riv_gdf.geometry.y

# -----------------------------------------------------------
# find nearest neighbor to simulated river network
# -----------------------------------------------------------

# get MartheField of simulated river network.
# use_imask set to False
riv_mf = MartheField('aff_r',mm.mlname+'.aff_r',mm,use_imask=False)

# save to shape and reload as gdf
simriv_shp = os.path.join('gis','sim_riv.shp')
riv_mf.to_shapefile(simriv_shp,epsg=2154)
simriv_gdf = gpd.read_file(simriv_shp)

# load cleaned shapefile with Qgis (without northern reach)
simriv_gdf = gpd.read_file(os.path.join('gis','sim_riv_clean.shp'))

# load full river network

fullriv_gdf = gpd.read_file(os.path.join('..','data','SIG','Cours-Eau_BDTopage.shp'))

# --------------------------------------------
# -- plot map of river pumping stations 
# --------------------------------------------

marker_dic = {'STEP' : 'o', 'IRRIGATION' :'v', 'INDUSTRIEL' :'*', 'AEP':'+'}
color_dic = {'STEP' :'grey', 'IRRIGATION' :'darkgreen',  'INDUSTRIEL' :'darkorange','AEP':'blue'}

# load shapefile 
basin_shp = os.path.join('..','data','SIG','BassinLizonne.shp')
basin_gdf = gpd.read_file(basin_shp)

# plot basin
fig, ax = plt.subplots(figsize=(5,4))
ax = basin_gdf.boundary.plot(ax=ax,color='black', linewidth=0.5)

# plot actual river network
ax=fullriv_gdf.plot(ax=ax,color='royalblue',alpha=1,linewidth=0.25)

# plot simulated river network
ax=simriv_gdf.plot(ax=ax,color='darkblue',alpha=0.5)


# plot with distinct marker by layer
for u in sorted(riv_gdf.Usage.unique()):
    riv_gdf.loc[riv_gdf.Usage==u].plot(marker=marker_dic[u],
                                     color=color_dic[u],
                                     markersize=14,
                                     label=u,
                                     ax=ax)

ax.legend(frameon=False,fontsize=7)

fig.savefig(os.path.join('figs','qriv.png'),dpi=300)

# --------------------------------------------
# -- get nearest simulated river cell  
# --------------------------------------------

# spatial join based on neared neighbor
riv_gdf = riv_gdf.sjoin_nearest(simriv_gdf, rsuffix='simriv', distance_col='dist')

# get node from river cell coordinates 
riv_gdf['node'] = [mm.get_node(x,y,layer=0)[0] for x,y in zip(riv_gdf.x_simriv,riv_gdf.y_simriv)]

# --------------------------------------------
# -- temporal disaggregation of river withdrawals by months
# --------------------------------------------

# set usage for 'pompe Riviere déconnectée' (remplissage retenue)
# necessary for montlhy disaggregation 
riv_gdf.loc[riv_gdf.Ressource == 'pompe Riviere déconnectée','use'] = 'REMP_RETENUE'

# qdf with single-indexed rows (date) and multi-indexed columns (node, layer, use)
qdf = riv_gdf.set_index(['use','node'])[ycols].T # gotcha!
qdf.index= pd.to_datetime(qdf.index,format='%Y')

# generate monthly time sequence at month start (MS)
start, end = pd.to_datetime(['2003-01-01','2019-12-31'])
months = pd.date_range(start, end, freq='MS') 

# update index, fill monthly records with annual values
dti_uni = qdf.index.union(months).drop_duplicates()
dfm = qdf.reindex(dti_uni).interpolate('pad',limit_direction='forward')

# generate df with monthly ratios of annual pumping 
dfw = dfm.copy()

# uniform disaggregation for Industriel and STEP uses
dfw.iloc[:,:]=1./12 

# negative weights for "STEP" (flowing to the river network)
dfw.iloc[:,dfw.columns.get_level_values('use')=='STEP'] = \
        dfw.loc[:,dfw.columns.get_level_values('use')=='STEP']*-1

# month-dependent disaggregation for agricultural withdrawals
dfw.iloc[:,dfw.columns.get_level_values('use')=='IRRIGATION']=0. 
dfw.loc[dfw.index.month==6,dfw.columns.get_level_values('use')=='IRRIGATION']=0.05 
dfw.loc[dfw.index.month==7,dfw.columns.get_level_values('use')=='IRRIGATION']=0.45
dfw.loc[dfw.index.month==8,dfw.columns.get_level_values('use')=='IRRIGATION']=0.45
dfw.loc[dfw.index.month==9,dfw.columns.get_level_values('use')=='IRRIGATION']=0.05

# month-dependent disaggregation for agricultural withdrawals
dfw.iloc[:,dfw.columns.get_level_values('use')=='REMP_RETENUE']=0. 
dfw.loc[dfw.index.month==1,dfw.columns.get_level_values('use')=='REMP_RETENUE']=0.33
dfw.loc[dfw.index.month==2,dfw.columns.get_level_values('use')=='REMP_RETENUE']=0.33
dfw.loc[dfw.index.month==3,dfw.columns.get_level_values('use')=='REMP_RETENUE']=0.33

# term-by-term product to get the final dataframe with monthly values (m3) 
dff = dfm*dfw

# --------------------------------------------
# -- spatial aggregation of river withdrawals by model cell
# --------------------------------------------

# ---- groupby by node + layer and perform spatial aggregation
dfg = dff.groupby(['node'],axis=1).sum()

# --------------------------------------------
# -- plot monthly records of surface withdrawls
# --------------------------------------------

fig, ax = plt.subplots(figsize=(8,3))
ax = dff.groupby(['use'],axis=1).sum()[['INDUSTRIEL','STEP','IRRIGATION','REMP_RETENUE','AEP']].plot.bar(
        stacked=True,
        color = ['darkorange','darkgrey','darkgreen','green','blue'],
        ax=ax)
ticklabels = ['']*len(dff.index)
ticklabels[::12] = [item.strftime('%m-%Y') for item in dff.index[::12]]
ax.xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
plt.gcf().autofmt_xdate()
ax.set_ylabel('Prélèvements surf. [m$^3$/mois]')
fig.savefig(os.path.join('figs','qsurf_records.png'),dpi=300)

dff_surf=dff

# --------------------------------------------
# -- write monthly files of river withdrawals
# --------------------------------------------

for m in months:
    # transpose to get current month df 
    # with the given month as column name and layer, node as row index
    dft = dfg.loc[dfg.index==m].T
    # drop na
    dft.dropna(subset=m,inplace=True)
    # number of days this month
    ndays = (m.replace(month = m.month % 12 +1, day = 1)-datetime.timedelta(days=1)).day
    # convert m3/month to m3/s, positive for river inflow, negative for withdrawals
    dft['V'] = dft[m]*-1./(ndays*86400)
    # get coordinates of node centers
    dft['XCC'] = mm.query_grid(node=dft.index.get_level_values('node').values, target='xcc').values
    dft['YCC'] = mm.query_grid(node=dft.index.get_level_values('node').values, target='ycc').values
    # output filename
    qfilename = os.path.join('prelevements_superficiels_mensuels',f'Q_{m.year}_{m.month:02d}.txt')
    # write to csv
    dft[['XCC','YCC','V']].to_csv(qfilename, sep='\t', header=False, index=False, float_format='%.6f')

# ==========================================================
# -- write .pastp file (thanks @ebuscarlet and @jpvergnes)
# ===========================================================

# simulation starts on day sim_start at 00:00:00
sim_start = datetime.datetime(2010, 7, 31)
# simulation ends on day sim_end at 23:59:59
sim_end = datetime.datetime(2019, 7, 31)

# days with resolution of the gw equation
n_gwsol_pmonth = 15 # number of solutions per month
days_gwsol = np.linspace(1,28,n_gwsol_pmonth).astype(int)


# sequence of time step starting dates 
sim_days = pd.date_range(sim_start, sim_end, freq='D')

pastp_file = 'Lizonne.pastp'

initlines = """  /DEBIT_RIVI/EDITION      I= 1;L= 0;F= 0;B= 0;
  /QECH_RIV_NAPP/EDITION   I= 1;
  /RUISSEL/EDITION         I= 1;Z= 0;C= 0;
  /DEBIT_DEBORD/EDITION    I= 1;
  /CHARGE/EDITION          I= 1;
  /DEBIT_RESID/EDITION     I= 1;
  /%SATURAT/EDITION        I= 1;
  /RECHARGE/ZONE_CLIM      Z=      *V=       0.4;
"""
i=0
with open(pastp_file, 'w', encoding='ISO-8859-1') as f:
    f.write("Modèle Lizonne 2021 - EAUX-SCARS\n")
    f.write(" /CALCUL_HDYNAM/ACTION    I= 0; File= Calcul_Hyd_hebdo.txt\n")
    f.write(' #<V7.8a># --- Fin du texte libre --- ;'
    'Ne pas modifier/retirer cette ligne\n')
    f.write(' *** Début de la simulation    à la date :'
            f' {sim_start.day:02d}/{sim_start.month:02d}/{sim_start.year:4d} ; ***\n')
    f.write(initlines)
    f.write('  /*****/***** Fin de ce pas\n')
    # iterate over sequence of time step ends
    for pas, date in enumerate(sim_days[1:]):
        f.write(f' *** Le pas :{pas+1}: se termine à la date '
                f': {date.day:02d}/{date.month:02d}/{date.year:4d} ; ***\n')
        # beginning of each month, set pumping data
        if date.day==1:
            qsoutfilename = os.path.join('prelevements_sout_mensuels',
                                     f'Q_{date.year}_{date.month:02d}.txt')
            f.write(f'  /DEBIT/LIST_MAIL         N: {qsoutfilename}\n')
            qrivfilename = os.path.join('prelevements_superficiels_mensuels',
                                     f'Q_{date.year}_{date.month:02d}.txt')
            f.write(f'  /Q_EXTER_RIVI/LIST_MAIL  N: {qrivfilename} '
                    '<Init_a_Zero> <X_Y_V> <Somm_Mail> <Keep_9999>\n')
            f.write('  /CALCUL_HDYNAM/ACTIO     I= 2;\n')
        # solve groundwater  hydrodynamics n_gwsol_pmonth times per months
        if date.day in days_gwsol[1:]:
            f.write('  /CALCUL_HDYNAM/ACTIO     I= 2;\n')
        # beginning of each hydrological year, set climate data and save head
        if (date.month==8 and date.day==1):
            f.write('  /RECHARGE/ZONE_CLIM      Z=      *V=         0;\n')
            f.write(f'  /METEO_PLUIE/FICHIER     N: SAFRAN_Lizonne/Plu+Neige_Lizonne_{date.year}_{date.year+1}\n')
            f.write(f'  /METEO_ETP/FICHIER       N: SAFRAN_Lizonne/ETP_Lizonne_{date.year}_{date.year+1}\n')
            f.write('  /FLUX_PLUV/FICH_METE     N=      *\n')
            f.write('  /FLUX_ETP/FICH_METE      N=      *\n')
            f.write(f'  /CHARGE/EDITION          I= 1;\n')
        f.write('  /*****/***** Fin de ce pas\n')
    f.write(' ***        :     : Fin de la simulation :            ; ***')

  

# -----------------------------------------------------------
#  write surface and groundwater withdrawals over cal period
# -----------------------------------------------------------

dff_gw = dff_gw.loc[dff_gw.index>sim_start]
dff_surf = dff_surf.loc[dff_surf.index>sim_start]

fig, axs = plt.subplots(2,1, sharex=True,figsize=(6,3))
ax1,ax2=axs

ax1 = dff_gw.groupby(['use'],axis=1).sum()[['AEP','IRRIGATION']].plot.bar(
    stacked=True,
    color=['darkblue','darkorange','darkgreen'],
    ax=ax1)
# Make most of the ticklabels empty so the labels don't get too crowded
ticklabels = ['']*len(dff.index)
# Every 12th ticklabel includes the year
ticklabels[::12] = [item.strftime('%m-%Y') for item in dff.index[::12]]
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
ax1.set_ylabel('[m$^3$/month]')
ax1.set_title('Groundwater withdrawals')
plt.gcf().autofmt_xdate()


ax2 = dff_surf.groupby(['use'],axis=1).sum()[['INDUSTRIEL','STEP','IRRIGATION','AEP']].plot.bar(
        stacked=True,
        color = ['darkorange','darkgrey','darkgreen','blue'],
        ax=ax2)
ticklabels = ['']*len(dff.index)
ticklabels[::12] = [item.strftime('%m-%Y') for item in dff.index[::12]]
ax2.xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
plt.gcf().autofmt_xdate()
ax2.set_ylabel('[m$^3$/month]')
ax2.set_title('Surface withdrawals')
fig.tight_layout()

fig.savefig(os.path.join('figs','q_surf_gw_records.png'),dpi=300)


# -----------------------------------------------------------
#  simplify soil map 
# -----------------------------------------------------------

# original soil zones map
zonep_org = MartheField('zonep','Lizonne.zonep_org',mm,use_imask=False)
zs2d_org = zonep_org.as_3darray()[0,:,:]

# convert to shapefile and reload 
zonep_org_shp = os.path.join('figs','zonep_org.shp')
zonep_org.to_shapefile(zonep_org_shp,layer=0)
zonep_org_gdf = gpd.read_file(zonep_org_shp)

# plot original soil zones + basin
fig, ax = plt.subplots(figsize=(5,4))
ax=zonep_org_gdf.plot(ax=ax, column='val', legend=True, categorical=True)
ax = basin_gdf.boundary.plot(ax=ax,color='black', linewidth=0.5)
fig.savefig(os.path.join('figs','zonep_org.png'), dpi=300)

# discard undesired soil zones 
keep_zone_values = [2,8]
zonep2d_simp = zs2d_org.copy()
zonep2d_simp[~np.isin(zonep2d_simp,keep_zone_values)]=np.nan

# fill gaps with nearest neighbor
from scipy import ndimage as nd
ind = nd.distance_transform_edt(np.isnan(zonep2d_simp), return_distances=False, return_indices=True)
zonep2d_simp = zonep2d_simp[tuple(ind)]

# write MartheField with simplified soil zones
zonep3d_simp =  zonep_org.as_3darray()
zonep3d_simp[0,:,:] = zonep2d_simp
zonep_simp = MartheField('zonep',zonep3d_simp,mm,use_imask=False)
zonep_simp.write_data()

# convert to shapefile and reload 
zonep_simp_shp = os.path.join('figs','zonep_simp.shp')
zonep_simp.to_shapefile(zonep_simp_shp,layer=0)
zonep_simp_gdf = gpd.read_file(zonep_simp_shp)

# plot simplified soil zones + basin
fig, ax = plt.subplots(figsize=(5,4))
ax=fullriv_gdf.plot(ax=ax,color='black',alpha=1,linewidth=0.25)
ax=zonep_simp_gdf.plot(ax=ax, column='val', legend=True, categorical=True)
ax = basin_gdf.boundary.plot(ax=ax,color='black', linewidth=0.5)
ax.set_xlim(ax_map.get_xlim())
ax.set_ylim(ax_map.get_ylim())
fig.savefig(os.path.join('figs','zonep_simp.png'), dpi=300)


# -----------------------------------------------------------
# plot observation data
# -----------------------------------------------------------
# load histo file and convert layer ids to 0-based 
histo_df = marthe_utils.read_histo_file(mm.mlname+'.histo')
histo_df['layer'] = histo_df.layer - 1

# convert to geopandas df
histo_gdf = gpd.GeoDataFrame(histo_df,
                       geometry = gpd.points_from_xy(histo_df['x'], 
                                                     histo_df['y']),
                       crs = 2154) # EPSG:RGF93

# get color from layer index (layers are 1-based in histo file!)
histo_gdf['color']=histo_gdf.layer.sub(1).replace(geol_color_dic)

# river discharge in blue
histo_gdf.loc[histo_gdf['type']=='Débit_Rivi','color']='red'

# plot basin
fig, ax = plt.subplots(figsize=(5,4))
ax = basin_gdf.boundary.plot(ax=ax,color='black', linewidth=0.5)

# plot simulated river network
ax=simriv_gdf.plot(ax=ax,color='darkblue',alpha=0.5)



histo_gdf['label2']= histo_gdf.label.replace(['P8284010', 'P8215010', 'P7250001', 'P7270001'],
                                       ['Lizonne','Belle','Sauvanie','Pude'])



# plot river discharge stations 
histo_gdf.loc[histo_gdf['type']=='Débit_Rivi'].plot(
        marker='^',
        c='red',
        edgecolor='black',
        markersize=22,
        label='River discharge station',
        ax=ax)
'''
histo_gdf.loc[histo_gdf['type']=='Débit_Rivi'].apply(
        lambda x: ax.annotate(
            text=x.label2,
            xy=(x.geometry.x,x.geometry.y),
            xytext=(3,-2),textcoords='offset points',
            ha='left',
            va='center',
            fontsize=10),
        axis=1)
'''
ss_obs_cal = ['07574X0014', '07345X0017', '07338X0017', '07338X0016', '07346X0083', '07345X0023', '07333X0027', '07346X0017']
histo_gdf = histo_gdf.loc[ss_obs_cal]

# plot head observations with distinct color by layer
for l,group in histo_gdf.loc[histo_gdf['type']=='Charge'].groupby('layer'):
    group.plot(marker='o',
               c=geol_color_dic[l],
               edgecolor='black',
               markersize=22,
               label=geol_id_dic[l],
               ax=ax)

ax.legend(title='',frameon=False,fontsize=7)


fig.savefig(os.path.join('figs','obs.png'),dpi=300)

# -----------------------------------------------------------
# plot observation data
# -----------------------------------------------------------
