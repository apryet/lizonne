import os, sys, shutil
import pandas as pd
import numpy as np
import datetime 

from pymarthe import MartheModel
from pymarthe.moptim import MartheOptim
from pymarthe.mfield import MartheField
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
import geopandas as gpd

# -- Read model
mm = MartheModel('Lizonne.rma',spatial_index=True)

# -----------------------------------------------------------
# load pumping data  
# -----------------------------------------------------------

# Read pumping data from excel file 
pump_data_xls = os.path.join('..','data','prelevements','Prelevements_V3-Lizonne.xlsx')

# selection of years (column names in excel file) for both groundwater and surface withdrawals
ycols = [y for y in range(2003,2021)]

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
dfg = dfg.groupby(dfg.index.month).mean()
dfg_sout = dfg.copy()

# --------------------------------------------------------------------------
# -- write monthly aquifer pumping files with inter-annual monthly means
# --------------------------------------------------------------------------
for m in dfg.index:
    # transpose to get current month df 
    # with the given month as column name and layer, node as row index
    dft = dfg.loc[dfg.index==m].T
    # drop na, get value, colonne,ligne,plan
    dft.dropna(subset=m,inplace=True)
    # number of days by month
    ndays = 365./12
    # convert m3/month to m3/s 
    dft['V'] = dft[m]*-1./(ndays*86400)
    # get indices for given node (back to 1-based !)
    dft['C'] = mm.query_grid(node=dft.index.get_level_values('node').values, target='j').values+1
    dft['L'] = mm.query_grid(node=dft.index.get_level_values('node').values, target='i').values+1
    dft['P'] = dft.index.get_level_values('layer').values+1  # back to 1-based !
    # output filename
    qfilename = os.path.join('prelevements_sout_mensuels',f'qsout_ia_{m:02d}.txt')
    # write to csv
    dft[['V','C','L','P']].to_csv(qfilename, sep='\t', header=False, index=False, float_format='%.6f')


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
dfg = dfg.groupby(dfg.index.month).mean()
dfg_surf = dfg.copy()

# --------------------------------------------
# -- write monthly files of river withdrawals
# --------------------------------------------

for m in dfg.index:
    # transpose to get current month df 
    # with the given month as column name and layer, node as row index
    dft = dfg.loc[dfg.index==m].T
    # drop na
    dft.dropna(subset=m,inplace=True)
    # number of days this month
    ndays = 365/12
    # convert m3/month to m3/s, positive for river inflow, negative for withdrawals
    dft['V'] = dft[m]*-1./(ndays*86400)
    # get coordinates of node centers
    dft['XCC'] = mm.query_grid(node=dft.index.get_level_values('node').values, target='xcc').values
    dft['YCC'] = mm.query_grid(node=dft.index.get_level_values('node').values, target='ycc').values
    # output filename
    qfilename = os.path.join('prelevements_superficiels_mensuels',f'qsurf_ia_{m:02d}.txt')
    # write to csv
    dft[['XCC','YCC','V']].to_csv(qfilename, sep='\t', header=False, index=False, float_format='%.6f')


# --------------------------------------------
# -- replicate template dir for each gcm/rcm
# --------------------------------------------
print('Settting up full simulation dirs')
cm_df = pd.read_excel('../data/DRIAS_Lizonne/clim_models.xlsx',index_col=0)

for cm_id in cm_df.index:
    # define sim dir and  remove if exists 
    sim_dir = os.path.join('..','sims',f'sim_{cm_id:02d}')
    print(sim_dir)
    if os.path.exists(sim_dir):
           shutil.rmtree(sim_dir)
    # copy from template 
    shutil.copytree(os.path.join('..','tpl_sim'),sim_dir)
    

