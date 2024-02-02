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
mm = MartheModel('Lizonne.rma')

# ===========================================================
# -- write .pastp file (thanks @ebuscarlet and @jpvergnes)
# ===========================================================

# simulation starts on day sim_start at 00:00:00
sim_start = datetime.datetime(1952, 7, 31)
# simulation ends on day sim_end at 23:59:59
sim_end = datetime.datetime(2099, 7, 31)

# sequence of time step starting dates 
sim_days = pd.date_range(sim_start, sim_end, freq='D')


cm_id = int(os.getcwd().split('_')[-1])
cm_df = pd.read_excel(os.path.join('..','..','data','DRIAS_Lizonne','clim_models.xlsx'),index_col=0)
gcm = cm_df.loc[cm_id,'GCM']
rcm = cm_df.loc[cm_id,'RCM']

initlines = """  /DEBIT_RIVI/EDITION      I= 1;L= 0;F= 0;B= 0;
  /QECH_RIV_NAPP/EDITION   I= 1;
  /RUISSEL/EDITION         I= 1;Z= 0;C= 0;
  /DEBIT_DEBORD/EDITION    I= 1;
  /CHARGE/EDITION          I= 1;
  /DEBIT_RESID/EDITION     I= 1;
  /%SATURAT/EDITION        I= 1;
  /RECHARGE/ZONE_CLIM      Z=      *V=       0.3;
"""


pastp_file = 'Lizonne.pastp'
#i=0

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
        # beginning of each month, set pumping data and solve for head
        if date.day==1:
            qsoutfilename = os.path.join('prelevements_sout_mensuels',
                                     f'qsout_ia_{date.month:02d}.txt')
            f.write(f'  /DEBIT/LIST_MAIL         N: {qsoutfilename}\n')
            qrivfilename = os.path.join('prelevements_superficiels_mensuels',
                                     f'qsurf_ia_{date.month:02d}.txt')
            f.write(f'  /Q_EXTER_RIVI/LIST_MAIL  N: {qrivfilename} '
                    '<Init_a_Zero> <X_Y_V> <Somm_Mail> <Keep_9999>\n')
            f.write('  /CALCUL_HDYNAM/ACTIO     I= 2;\n')
        # beginning of each hydrological year, set climate data, solve for head, and save head
        if (date.month==8 and date.day==1):
            f.write('  /RECHARGE/ZONE_CLIM      Z=      *V=         0;\n')
            if date < pd.to_datetime('2005-08-01'):
                f.write(f'  /METEO_PLUIE/FICHIER   N: ../../data/DRIAS_Lizonne/{gcm}/{rcm}/historical/prtotAdjust_{gcm}_{rcm}_MF-ADAMONT-SAFRAN-1980-2011_historical_day_{date.year}_{date.year+1}\n')
                f.write(f'  /METEO_ETP/FICHIER     N: ../../data/DRIAS_Lizonne/{gcm}/{rcm}/historical/evspsblpotAdjust_Hg0175_{gcm}_{rcm}_MF-ADAMONT-SAFRAN-1980-2011_historical_day_{date.year}_{date.year+1}\n')
            else : 
                f.write(f'  /METEO_PLUIE/FICHIER     N: ../../data/DRIAS_Lizonne/{gcm}/{rcm}/rcp85/prtotAdjust_{gcm}_{rcm}_MF-ADAMONT-SAFRAN-1980-2011_rcp85_day_{date.year}_{date.year+1}\n')
                f.write(f'  /METEO_ETP/FICHIER       N: ../../data/DRIAS_Lizonne/{gcm}/{rcm}/rcp85/evspsblpotAdjust_Hg0175_{gcm}_{rcm}_MF-ADAMONT-SAFRAN-1980-2011_rcp85_day_{date.year}_{date.year+1}\n')
            f.write('  /FLUX_PLUV/FICH_METE     N=      *\n')
            f.write('  /FLUX_ETP/FICH_METE      N=      *\n')
            f.write('  /CALCUL_HDYNAM/ACTIO     I= 2;\n')
            f.write(f'  /CHARGE/EDITION          I= 1;\n')
        f.write('  /*****/***** Fin de ce pas\n')
    f.write(' ***        :     : Fin de la simulation :            ; ***')


# run model (so as to update historiq.prn and chasim.out)
mm.run_model()
shutil.copy('chasim.out','chasim.out.full')


