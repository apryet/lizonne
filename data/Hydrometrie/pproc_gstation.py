"""
Preprocessing gaging station 
"""


# ---- Import usefull packages 
import os, sys
import numpy as np
import pandas as pd
import re
dev_ws = os.path.normpath(r"E:\EauxSCAR\pymarthe_dev")
sys.path.append(dev_ws)
import shutil
from pymarthe.utils import marthe_utils

# ---- Build a 'obs' folder
if os.path.exists('obs'):
    shutil.rmtree('obs')

os.makedirs('obs')


# ---- Set time window
start = '2003-07-31'
end = '2019-07-31'


# -- Iterate over (excel) files in current folder
for f in os.listdir():
    if f.endswith('.xlsx'):
        # -- extract locnme
        name = f.split('-')[0]
        # -- Read excel sheet
        df = pd.read_excel(f, sheet_name=0, usecols=list(range(7)))
        # -- Add DatetimeIndex and subset in time window
        df_ss = df.set_index(pd.to_datetime(df['Date (TU)'])).loc[start:end]
        # -- Build output filename (locnme)
        out = os.path.join('obs', f'{name}.dat')
        marthe_utils.write_obsfile(date=df_ss.index, value=df_ss.iloc[:,1], obsfile=out)


# ---- Export gaging station network as shapefile
import geopandas as gpd

data = [[485597,6470094,0,'P8284010','Lizonne'],
        [501993,6486335,0,'P8215010', 'Belle'],
        [488540,6479509,0,'P7250001','Pude'],
        [488706,6470782,0,'P7270001','Sauvanie']]

df = pd.DataFrame(data, columns = ['x','y','layer','id','name'])

gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x,df.y), crs='epsg:2154')

gdf.to_file('gaging_stations.shp')