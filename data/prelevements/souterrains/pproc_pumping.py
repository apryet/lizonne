"""
Preprocessing groundwater withdraw
"""


# ---- Import usefull packages 
import os, sys
import numpy as np
import pandas as pd
import re
dev_ws = os.path.normpath(r"E:\EauxSCAR\pymarthe_dev")
sys.path.append(dev_ws)
from pymarthe import MartheModel


# ---- Extract str date from file name (regex)
def extract_date_from_filename(f):
    """
    """
    re_date = r'prelevements_(\d{4})_(\d{2})_(\d{2}).txt'
    str_date = '-'.join(re.findall(re_date, f)[0])
    return str_date


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



# ---- Read all pumping .txt file and write aggregate pumping 
agg_dfs = []
dfs = []
for f in os.listdir('all_pumps'):
    if f.endswith('.txt'):
        df = pd.read_csv(os.path.join('all_pumps',f),
                         delim_whitespace=True,
                         names=['x','y','layer','q'])
        dfs.append(df)
        agg_df = spatial_aggregation(mm, x=df.x, y=df.y, layer=df.layer.sub(1), value=df.q, agg='sum')
        agg_df = agg_df.loc[[mm.all_active(n) for n in agg_df.node]]
        agg_df['layer'] = agg_df.layer.add(1)
        agg_df[['x','y','layer','agg_value']].to_csv(os.path.join('agg_pumps',f),
                                                     header=False, index = False,
                                                     sep='\t')
        agg_dfs.append(agg_df)



# ---- Export pumping points as shapefile
import geopandas as gpd

# Before aggregation
df = pd.concat(dfs, ignore_index=True)
gdf = gpd.GeoDataFrame(df.to_records(index=False),
                       geometry = gpd.points_from_xy(df['x'], df['y']),
                       crs = 'epsg:2154')
gdf.to_file('before_agg_prelevements.shp')

# After aggregation
agg_df = pd.concat(agg_dfs, ignore_index=True)
agg_gdf = gpd.GeoDataFrame(agg_df.to_records(index=False),
                       geometry = gpd.points_from_xy(agg_df['x'], agg_df['y']),
                       crs = 'epsg:2154')

agg_gdf.to_file('prelevements.shp')











