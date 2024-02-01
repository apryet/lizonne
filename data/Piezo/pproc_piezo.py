"""
Preprocessing observed heads
"""


# ---- Import usefull packages 
import os, sys
import numpy as np
import pandas as pd
import shutil

dev_ws = os.path.normpath(r"E:\EauxSCAR\pymarthe_dev")
sys.path.append(dev_ws)
from pymarthe.utils import marthe_utils

# ---- Read base data in excel file
df = pd.read_excel('Piezo_Eaux-SCARS-Lizonne.xlsx', sheet_name='Mesures')
df.set_index('Date', inplace=True)

# ---- Build a 'obs' folder
if os.path.exists('obs'):
    shutil.rmtree('obs')

os.makedirs('obs')

# ---- Set time window
start = '2003-07-31'
end = '2019-07-31'

# ---- Write observation file for each set of observation
for name in df['Ouvrage'].unique():
    # -- Subset by name and sort by index (date)
    df_ss = df.loc[df['Ouvrage'] == name].sort_index()
    # -- Subset just data in time window
    df_ss = df_ss.loc[start:end]
    # -- Build output filename (locnme)
    out = os.path.join('obs', f'{name}.dat')
    marthe_utils.write_obsfile(date=df_ss.index, value=df_ss.Mesure, obsfile=out)
    # -- Verify that is readable
    try:
        test = marthe_utils.read_obsfile(out)
    except:
        print('ERROR : Could not re-read observation file just created.')


# ---- Export observation well network as shapefile
import geopandas as gpd
df = pd.read_excel('Piezo_Eaux-SCARS-Lizonne.xlsx', sheet_name='Piezo-Lizonne')
gdf = gpd.GeoDataFrame(df.iloc[:,:6].to_records(index=False),
                       geometry = gpd.points_from_xy(df['X_L2e'], df['Y_L2e']),
                       crs = 'epsg:27572')

gdf.to_file('piezo_L2E.shp')
# has to be projected in 2154 - RG93
