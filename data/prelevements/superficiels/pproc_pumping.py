"""
Preprocessing surface withdraws
"""


# ---- Import usefull packages 
import os, sys
import numpy as np
import pandas as pd
import re
dev_ws = os.path.normpath(r"E:\EauxSCAR\pymarthe_dev")
sys.path.append(dev_ws)
from pymarthe import MartheModel



def m3y_to_m3s(v):
    """
    """
    if v == 0:
        return 0
    elif v > 0:
        return -v / (365*86400)
    elif v < 0:
        return v / (365*86400)


def spatial_aggregation(mm, x, y, layer, value, agg = 'sum'):
    """
    """
    # ---- Sample node ids for each xylayer-points

    nodes = mm.get_node(x,y,layer)

    # nodes = [] 
    # for ix,iy,ilay in zip(x,y,layer):
    #     inodes = list(mm.spatial_index.intersection((ix,iy)))
    #     inode = [n for n in inodes if (ilay-1)*mm.ncpl < n < ilay * mm.ncpl][0]
    #     nodes.append(inode)

    # ---- Store in DataFrame
    df = pd.DataFrame( { 'node': nodes,
                         'x': x,
                         'y': y,
                         'layer' : layer,
                         'value' : value  } )

    # ---- Groupby by node and perform aggregation
    agg_df = df.groupby('node', as_index=False).agg({'value': agg})
    agg_df.rename({'value':'agg_value'}, axis = 1, inplace=True)

    # ---- Merge agg/initial DataFrames 
    res = pd.merge(df,agg_df).drop_duplicates('node').drop('value',axis=1)

    # ---- Return aggregate values in DataFrame
    return res



# ---- Read Lizonne model
rma_path = os.path.join('..','..','..','lizonne_v1','Lizonne.rma')
si = rma_path.replace('.rma', '_si')
mm = MartheModel(rma_path, si)




# ---- Read superficial withdraw
filename = 'prelevements_sup_simplifies.csv'
raw_df = pd.read_csv(filename)

# ---- Get xy-coordinates and date columns
x, y = raw_df['X_L93New'], raw_df['Y_L93New']
date_cols = raw_df.columns.str.extract(r"(\d{4})", expand=False).dropna()

# ---- Write external texte file with aggregated transformed value for each simulated year
dfs = []
for date in date_cols:
    f = os.path.join('agg_pumps',f'prelevements_sup_{date}.txt')
    agg_df = spatial_aggregation(mm, x, y, layer=0, value=raw_df[date], agg='sum')
    agg_df['value'] = agg_df['agg_value'].transform(m3y_to_m3s)
    agg_df[['x','y','value']].to_csv(f, header=False, index = False, sep='\t')
    df_ss = agg_df[['node','x','y']]
    df_ss[f'val_{date}'] = agg_df['value']
    dfs.append(df_ss)



# ---- Export pumping points as shapefile
import geopandas as gpd

# After aggregation
df = pd.concat(dfs, axis=1)
df = df.loc[:,~df.columns.duplicated()]
gdf = gpd.GeoDataFrame(df.to_records(index=False),
                       geometry = gpd.points_from_xy(df['x'], df['y']),
                       crs = 'epsg:2154')
gdf.to_file('prelevements_sup_aggregated.shp')












