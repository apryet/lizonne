import os, glob
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt 
from statsmodels.distributions.empirical_distribution import ECDF

plt.rc('font', family='serif', size=11)


srcdir = os.path.join('DRIAS_Lizonne','CNRM-CERFACS-CNRM-CM5','CNRM-ALADIN63')


def read_var(srcdir, varname):
    hist_files = [f for f in glob.iglob(os.path.join(srcdir,'historical',f'{varname}*'))]
    rcp_files = [f for f in glob.iglob(os.path.join(srcdir,'rcp85',f'{varname}*'))]
    files = hist_files + rcp_files
    df_list=[]
    for f in files : 
        df=pd.read_csv(f,delim_whitespace=True,
                       header=None,index_col=-1,skiprows=1)
        df.index=pd.to_datetime(df.index,format='%d/%m/%Y')
        df.index.name='date'
        df = df.mean(axis=1)
        df_list.append(df)
    return(pd.concat(df_list))


etp = read_var(srcdir,'evspsblpotAdjust')
ptot = read_var(srcdir,'prtotAdjust')
# read histoclim
prn = pd.read_csv(os.path.join('histoclim.prn'),delimiter=r'\s*\t\s*',
                  encoding='latin-1',engine='python',skiprows=2,skipfooter=34)
prn.index= pd.to_datetime(prn.iloc[:,0],format='%d/%m/%Y')
prn.index.name = 'date'
rech = prn['Infiltration']
rech.name= 'rech'

# all clim data 
clim = pd.DataFrame({'ptot':ptot,'etp':etp})
clim = clim.merge(rech,left_index=True, right_index=True)

# annual aggregation (remove incomplete 1st and last years)
climy = clim.groupby(pd.Grouper(freq='Y')).sum().iloc[1:-1,:]

# --- define "current" and "future" periods, and subset
cper_start, cper_end  = pd.to_datetime('1990-01-01'), pd.to_datetime('2019-12-31')
#fper_start , fper_end =  pd.to_datetime('2060-01-01'), pd.to_datetime('2089-12-31')
fper_start , fper_end =  pd.to_datetime('2040-01-01'), pd.to_datetime('2069-12-31')

# subset to fit trend line 
rcp_start =  pd.to_datetime('2005-08-01')
fit_start = rcp_start 
fit_start =sclimy = climy.loc[climy.index>fit_start]
t= (sclimy.index.values - sclimy.index.values.min())/np.timedelta64(1,'D')

# ---plot P, PET, rech
fig,axs=plt.subplots(3,1,figsize=(8,5))
# plot ptot 
axs[0].bar(climy.index,climy['ptot'],color='royalblue',width=180)
# fit trend line 
pfit = np.polyfit(t,sclimy.ptot.values, 1)
axs[0].plot(sclimy.index, np.polyval(pfit,t),color='royalblue',ls='--')
axs[0].set_ylabel('Ptot [mm/y]')
axs[0].set_xticklabels([])
# plot etp 
axs[1].bar(climy.index,climy['etp'],color='darkorange',width=180)
pfit = np.polyfit(t,sclimy.etp.values, 1)
axs[1].plot(sclimy.index, np.polyval(pfit,t),color='darkorange',ls='--')
axs[1].set_ylabel('ETP [mm/y]')

# plot rech
axs[2].bar(climy.index,climy['rech'],color='tan',width=180)
pfit = np.polyfit(t,sclimy.rech.values, 1)
axs[2].plot(sclimy.index, np.polyval(pfit,t),color='tan',ls='--')
axs[2].set_ylabel('Recharge [mm/y]')

for ax in axs:
    ax.axvline(pd.to_datetime(rcp_start),alpha=0.5,ls=':',color='black')
    ax.axvspan(cper_start,cper_end,color='grey',alpha=0.3)
    ax.axvspan(fper_start,fper_end,color='grey',alpha=0.3)

fig.savefig(os.path.join('figs','climy.pdf'),dpi=300)


cclimy = climy.loc[(climy.index.values>cper_start) & (climy.index.values< cper_end)]
fclimy = climy.loc[(climy.index>fper_start) & (climy.index< fper_end)]

# --- ECDF

cecdf = ECDF(cclimy.rech)
fecdf = ECDF(fclimy.rech)


fig,ax=plt.subplots(1,1,figsize=(4,4))
ax.plot(cecdf.x,cecdf.y,label='current period',color='grey')
ax.plot(fecdf.x,fecdf.y,label='future period',color='darkred')
ax.set_xlabel('Recharge [mm/y]')
ax.set_ylabel('ECDF [-]')
'''
ax.axvline(np.quantile(cclimy.rech,0.1),alpha=0.5,ls=':',color='grey')
ax.axvline(np.quantile(cclimy.rech,0.5),alpha=0.5,ls=':',color='grey')
ax.axvline(np.quantile(cclimy.rech,0.9),alpha=0.5,ls=':',color='grey')
ax.axvline(np.quantile(fclimy.rech,0.1),alpha=0.5,ls=':',color='darkred')
ax.axvline(np.quantile(fclimy.rech,0.5),alpha=0.5,ls=':',color='darkred')
ax.axvline(np.quantile(fclimy.rech,0.9),alpha=0.5,ls=':',color='darkred')
'''
ax.axhline(0.1,alpha=0.5,ls=':',color='black')
ax.axhline(0.5,alpha=0.5,ls=':',color='black')
ax.axhline(0.9,alpha=0.5,ls=':',color='black')

fig.tight_layout()

fig.savefig(os.path.join('figs','ecdf.pdf'),dpi=300)


# histograms


fig,axs=plt.subplots(1,3,figsize=(6,3))

for col,ax in zip(cclimy.columns,axs.ravel()):
    vmin = min(cclimy[col].min(),fclimy[col].min())
    vmax = max(cclimy[col].max(),fclimy[col].max())
    bins = np.linspace(vmin,vmax,10)
    cclimy[col].hist(ax=ax,color='grey',bins=bins,grid=False,alpha=0.8,label='current')
    fclimy[col].hist(ax=ax,color='darkred',bins=bins,grid=False,alpha=0.5, label='future')
    ax.set_xlabel(f'{col} [mm/y]')

axs[2].legend()

fig.tight_layout()

fig.savefig(os.path.join('figs','histos.pdf'),dpi=300)


# --- histogram and box plots for recharge 

rech_df = pd.DataFrame({'cur.':cclimy.rech.values,'fut.':fclimy.rech.values})

fig, axs = plt.subplots(2, figsize=(5, 5), sharex=True, # Common x-axis
                       gridspec_kw={"height_ratios": (.7, .3)})
# histogram
rech_df['cur.'].hist(ax=axs[0], bins=bins, color='grey', grid=False,alpha=0.8,label='current')
rech_df['fut.'].hist(ax=axs[0], bins=bins, color='darkred', grid=False,alpha=0.8,label='future')
axs[0].legend()

# boxplot 
medianprops = dict(color='black')
meanprops = dict(linestyle=None,marker='+', markeredgecolor='black', markerfacecolor='black')

data = [rech_df['fut.'].values,rech_df['cur.'].values]
bplots = axs[1].boxplot(data,vert=False,
                        widths=0.4,
                        whis=(10,90),
                        showmeans=True,
                        medianprops=medianprops,
                        meanprops=meanprops,
                        patch_artist=True
                        )

axs[1].set_yticklabels(['fut','cur.'])

# fill with colors
colors = ['darkred','grey']
for patch, color in zip(bplots['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)

axs[1].set_xlabel('Recharge [mm/y]')

fig.tight_layout()


fig.savefig(os.path.join('figs','rech.pdf'),dpi=300)

# --- identify deciles of interest

qs = [0.1,0.5,0.9]

crech_qs = [np.quantile(cclimy.rech,q) for q in qs ]
frech_qs = [np.quantile(fclimy.rech,q) for q in qs ]

cqs_yrs = [ cclimy.index[np.argmin(np.abs(cclimy.rech - rech))].year for rech in crech_qs]
fqs_yrs = [ fclimy.index[np.argmin(np.abs(fclimy.rech - rech))].year for rech in frech_qs]

qs_labels = ['Q10','Q50','Q90']
cyrs  = pd.DataFrame({'rech':crech_qs,'year':cqs_yrs},index=qs_labels)
fyrs = pd.DataFrame({'rech':frech_qs,'year':fqs_yrs},index=qs_labels)

simyrs = pd.concat([cyrs,fyrs],keys=['current','future'])

simyrs.to_excel('simyrs.xlsx')
simyrs.to_csv('simyrs.csv')

