import os, shutil
import numpy as np
import pandas as pd 
import pymarthe
from pymarthe.utils import marthe_utils,pest_utils

def parfacmul(facfile, org_parfile, parfile, keys):
    # load original pumping values 
    parkmi, parvals = pest_utils.parse_mlp_parfile(org_parfile,
                                                  keys=keys,
                                                  value_col=1,btrans='none')
    parvals.index = parkmi
    # load pumping factors 
    fackmi, facvals = pest_utils.parse_mlp_parfile(facfile,
                                                 keys=keys,
                                                 value_col=1,btrans='none')
    facvals.index= fackmi.droplevel(0)
    # apply factors
    # newparvals = facvals*parvals (highly dangerous, resort parvals !)
    # stackoverflow options :
    # https://stackoverflow.com/questions/41493177/pandas-multiply-dataframes-with-multiindex-and-overlapping-index-levels
    # https://stackoverflow.com/questions/61338300/divide-pandas-dataframe-with-multiindex-by-another-dataframe-with-smaller-multii
    # but I prefer this other approach :
    newparvals = parvals*facvals.loc[parvals.index.droplevel(0)].values

    # write new parameter value file 
    parnmes = ['__'.join(list(map(str, items))) for items in parkmi]
    newparam_df = pd.DataFrame({'parnme':parnmes, 'value':newparvals.values})
    pest_utils.write_mlp_parfile(parfile, newparam_df, trans='none', value_col='value')


def write_pastp(pastp_file, sim_start, sim_end, gcm, rcm):
    datestrfmt = "%d/%m/%Y"
    # sequence of time step starting dates 
    sim_days = pd.date_range(sim_start, sim_end, freq='D')
    init_date = sim_days[0]-pd.Timedelta(1, "d")
    # days with resolution of the gw equation
    n_gwsol_pmonth = 15 # number of solutions per month
    days_gwsol = np.linspace(1,28,n_gwsol_pmonth).astype(int)
    with open(pastp_file, 'w', encoding='ISO-8859-1') as f:
        f.write("Modèle Lizonne 2021 - EAUX-SCARS\n")
        f.write(" /CALCUL_HDYNAM/ACTION    I= 0; File= Calcul_Hyd_hebdo.txt\n")
        f.write(' #<V7.8a># --- Fin du texte libre --- ;'
        'Ne pas modifier/retirer cette ligne\n')
        f.write(' *** Début de la simulation à la date : ' + init_date.strftime(datestrfmt) + ' ; ***\n')
        f.write('  /*****/***** Fin de ce pas\n')
        # iterate over sequence of time step ends
        for pas, date in enumerate(sim_days):
            f.write(f' *** Le pas :{pas+1}: se termine à la date : ' + date.strftime(datestrfmt) + ' ; ***\n')
            # beginning of each month, set pumping data
            if date.day==1:
                qsoutfilename = os.path.join('prelevements_sout_mensuels',
                                         f'qsout_ia_{date.month:02d}.txt')
                f.write(f'  /DEBIT/LIST_MAIL         N: {qsoutfilename}\n')
                qrivfilename = os.path.join('prelevements_superficiels_mensuels',
                                         f'qsurf_ia_{date.month:02d}.txt')
                f.write(f'  /Q_EXTER_RIVI/LIST_MAIL  N: {qrivfilename} '
                        '<Init_a_Zero> <X_Y_V> <Somm_Mail> <Keep_9999>\n')
                f.write('  /CALCUL_HDYNAM/ACTIO     I= 2;\n')
            # solve groundwater  hydrodynamics n_gwsol_pmonth times per months
            if date.day in days_gwsol[1:]:
                f.write('  /CALCUL_HDYNAM/ACTIO     I= 2;\n')
            # beginning of each hydrological year, set climate data
            if (date.month==8 and date.day==1):
                f.write('  /RECHARGE/ZONE_CLIM      Z=      *V=         0;\n')
                if date < pd.to_datetime('2005-08-01'):
                    f.write(f'  /METEO_PLUIE/FICHIER   N: DRIAS_Lizonne/{gcm}/{rcm}/historical/prtotAdjust_{gcm}_{rcm}_MF-ADAMONT-SAFRAN-1980-2011_historical_day_{date.year}_{date.year+1}\n')
                    f.write(f'  /METEO_ETP/FICHIER     N: DRIAS_Lizonne/{gcm}/{rcm}/historical/evspsblpotAdjust_Hg0175_{gcm}_{rcm}_MF-ADAMONT-SAFRAN-1980-2011_historical_day_{date.year}_{date.year+1}\n')
                else : 
                    f.write(f'  /METEO_PLUIE/FICHIER     N: DRIAS_Lizonne/{gcm}/{rcm}/rcp85/prtotAdjust_{gcm}_{rcm}_MF-ADAMONT-SAFRAN-1980-2011_rcp85_day_{date.year}_{date.year+1}\n')
                    f.write(f'  /METEO_ETP/FICHIER       N: DRIAS_Lizonne/{gcm}/{rcm}/rcp85/evspsblpotAdjust_Hg0175_{gcm}_{rcm}_MF-ADAMONT-SAFRAN-1980-2011_rcp85_day_{date.year}_{date.year+1}\n')
                f.write('  /FLUX_PLUV/FICH_METE     N=      *\n')
                f.write('  /FLUX_ETP/FICH_METE      N=      *\n')
            f.write('  /*****/***** Fin de ce pas\n')
        f.write(' ***        :     : Fin de la simulation :            ; ***')



def opt_pastp():
    # --- write pastp 
    cm_id = int(os.getcwd().split('_')[-1])
    year = int(os.getcwd().split('_')[-2])
    cm_df = pd.read_excel('../data/DRIAS_Lizonne/clim_models.xlsx',index_col=0)
    gcm = cm_df.loc[cm_id,'GCM']
    rcm = cm_df.loc[cm_id,'RCM']
    # start of warm up period 
    sim_start = pd.to_datetime(str(year-2)+'-08-01')
    # end of simulation period
    sim_end = pd.to_datetime(str(year)+'-09-30')
    # copy files 
    if os.path.exists('DRIAS_Lizonne'):
        shutil.rmtree('DRIAS_Lizonne')
    shutil.copytree(os.path.join('..','data','DRIAS_Lizonne',f'{gcm}',f'{rcm}'),os.path.join('DRIAS_Lizonne',f'{gcm}',f'{rcm}'))
    # generate pastpfile 
    pastp_file = 'Lizonne.pastp'
    write_pastp(pastp_file, sim_start, sim_end, gcm, rcm)


def opt_parfacmul():
    # --- processing aquifer pumping 
    parfacmul(facfile = os.path.join('pest','par','aqpumpfac.dat'), 
              org_parfile = os.path.join('pest','par','aqpump.dat.org'),
              parfile = os.path.join('pest','par','aqpump.dat'),
              keys=['boundname','layer','istep']
              )
    # --- processing river pumping 
    parfacmul(facfile = os.path.join('pest','par','rivpumpfac.dat'), 
              org_parfile = os.path.join('pest','par','rivpump.dat.org'),
              parfile = os.path.join('pest','par','rivpump.dat'),
              keys=['prefix','aff_r','istep']
              )


def opt_postproc(mm=None,out_dir=os.path.join('pest','sim')):
    station_locs = ['P8284010', 'P8215010', 'P7250001', 'P7270001']
    monthly_m3sec_to_m3 = (365.25/12)*86400 # (days in a month) * (seconds in a day)
    daily_m3sec_to_m3day = 86400 # (seconds in a day)
    # load simulated values 
    prn_df = marthe_utils.read_prn()
    # load aquifer pumping data from MartheModel 
    if mm is None:
        mm = pymarthe.MartheModel('Lizonne.rma', spatial_index = 'Lizonne_si')
    # get simulation period 
    sim_start,sim_end = mm.get_time_window()
    # optimization starts at the beginning of the last year of the simulation period
    opt_start = pd.to_datetime(str(sim_end.year))

    # --- Get total pumping volumes over simulation period (m3)
    # processing aquifer pumping
    mm.load_prop('aqpump')
    aqpump_df = mm.prop['aqpump'].data
    ap_dates = mm.mldates[aqpump_df.istep]
    aqpump_ss = aqpump_df.loc[ap_dates >= opt_start] # subset to optimization period
    aqpump_value = (aqpump_ss.value*(-1*monthly_m3sec_to_m3)).sum() # monthly pumping rates in m3/s to m3
    pest_utils.write_simfile([sim_end], [aqpump_value], os.path.join(out_dir,'tot_aqpump.dat'))

    # processing river pumping
    mm.load_prop('rivpump')
    rivpump_df = mm.prop['rivpump'].data
    rp_dates = mm.mldates[rivpump_df.istep]
    rivpump_ss = rivpump_df.loc[rp_dates >= opt_start] # subset to optimization period
    rivpump_value = (rivpump_ss.value*(-1*monthly_m3sec_to_m3)).sum() # monthly pumping rates in m3/s to m3
    pest_utils.write_simfile([sim_end], [rivpump_value], os.path.join(out_dir,'tot_rivpump.dat'))

    # total pumping 
    tot_value = aqpump_value + rivpump_value
    pest_utils.write_simfile([sim_end], [tot_value], os.path.join(out_dir,'tot_pump.dat'))

    # --- Get river deficit at gaging stations

    # load DOE data 
    doe_df = pd.read_csv('doe.dat',index_col=0)
    # get simulated records 
    qsim_df = prn_df.loc[:,(slice(None), station_locs)].loc[:,'Débit_Rivi']
    qsim_df = qsim_df.loc[qsim_df.index >= opt_start] # subset to optimization period

    # compute discharge deficit (m3 from m3/s at the daily frequency)
    deficit = (doe_df['alerte'] - qsim_df).mul(daily_m3sec_to_m3day*qsim_df.le(doe_df['alerte'])).sum()
    #write total deficit
    pest_utils.write_simfile([sim_end], [deficit.sum()], os.path.join(out_dir,f'deficit_tot.dat'))
    for loc in station_locs:
            value = deficit.loc[loc]  # m3
            pest_utils.write_simfile([sim_end], [value], os.path.join(out_dir,f'deficit_{loc}.dat'))



