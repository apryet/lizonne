import os, glob
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt 
from matplotlib.lines import Line2D
plt.rc('font', family='serif', size=11)


# --- define "historical" and "future" periods, and subset
# DRIAS "reference period for historical climate"
cper_start, cper_end  = pd.to_datetime('1976-01-01'), pd.to_datetime('2005-12-31')
# DRIAS "mid term horizon"
#fper_start , fper_end =  pd.to_datetime('2041-01-01'), pd.to_datetime('2070-12-31')
# DRIAS "long term horizon"
fper_start , fper_end =  pd.to_datetime('2069-01-01'), pd.to_datetime('2098-12-31')

def read_var(srcdir, varname):
    hist_files = [f for f in glob.iglob(os.path.join(srcdir,'historical',f'{varname}*'))]
    rcp_files = [f for f in glob.iglob(os.path.join(srcdir,'rcp85',f'{varname}*'))]
    safran_files = [f for f in glob.iglob(os.path.join(srcdir,f'{varname}*'))]
    files = hist_files + rcp_files + safran_files
    df_list=[]
    for f in files : 
        df=pd.read_csv(f,delim_whitespace=True,
                       header=None,index_col=-1,skiprows=1)
        df.index=pd.to_datetime(df.index,format='%d/%m/%Y')
        df.index.name='date'
        df = df.mean(axis=1)
        df_list.append(df)
    return(pd.concat(df_list))

# read safran data
safran_dir = os.path.join('data','SAFRAN_Lizonne')
ptot = read_var(safran_dir,'Plu+Neige')
pet = read_var(safran_dir,'ETP')
safran = pd.DataFrame({'ptot':ptot,'pet':pet,'rech':np.nan })
safrany = safran.groupby(pd.Grouper(freq='Y')).sum().iloc[1:-1,:]

# load list of climate models 
cm_df = pd.read_excel(os.path.join('data','DRIAS_Lizonne','clim_models.xlsx'),index_col=0)

climy_dic={}
for i in cm_df.index:
    # read climate data
    gcm = cm_df.loc[i,'GCM']
    rcm = cm_df.loc[i,'RCM']
    cmdata_dir = os.path.join('data','DRIAS_Lizonne',f'{gcm}',f'{rcm}')
    pet = read_var(cmdata_dir,'evspsblpotAdjust')
    ptot = read_var(cmdata_dir,'prtotAdjust')
    # read histoclim
    prn = pd.read_csv(os.path.join('sims',f'sim_{i:02d}','histoclim.prn'),delimiter=r'\s*\t\s*',
                      encoding='latin-1',engine='python',skiprows=2,skipfooter=34)
    prn.index= pd.to_datetime(prn.iloc[:,0],format='%d/%m/%Y')
    prn.index.name = 'date'
    prn_ss = prn[['Ruissellement','Infiltration']]
    prn_ss.columns= ['runoff','rech']
    # all clim data 
    clim = pd.DataFrame({'ptot':ptot,'pet':pet})
    clim = clim.merge(prn_ss,left_index=True, right_index=True)
    # summer runoff
    clim['srunoff']=clim.runoff.values.copy()
    clim.loc[~clim.index.month.isin([6,7,8,9]),'srunoff']=0
    # annual aggregation (remove incomplete 1st and last years)
    # annual rech of year Y from 1-oct-(Y-1) to 30-sept-(Y)
    start = pd.Timestamp(f'{clim.index.min().year+1}-10-01')
    end = pd.Timestamp(f'{clim.index.max().year-1}-10-01')
    bins = pd.date_range(start,end, freq='12MS')
    climy = clim.groupby(pd.cut(clim.index, bins=bins)).sum()
    climy.index=bins[1:] # year Y
    # add to compil dic 
    climy_dic[i]=climy


'''
# test 

dates = pd.date_range('1980-08-01','1983-07-31')
df = pd.DataFrame({'v':[0]},index=dates)

df.loc['1981-09-01']=1
df.loc['1981-11-01']=1
df.loc['1982-01-05']=5

# not effective
df.groupby(pd.Grouper(freq="12MS", offset="10MS",origin='epoch')).sum()

# effective 
start = pd.Timestamp(f'{df.index.min().year+1}-10-01')
end = pd.Timestamp(f'{df.index.max().year-1}-10-01')
bins = pd.date_range(start,end, freq='12MS')
dfy = df.groupby(pd.cut(df.index, bins=bins)).sum()
dfy.index=bins[1:] # year Y
print(dfy)

'''

# concat to single df 
climy = pd.concat(climy_dic.values(),keys=climy_dic.keys(),axis=1)
climy.to_csv('climy.csv')
# subset to fit trend line
rcp_start= pd.to_datetime('2005-08-01') # rcp start 
fit_start = rcp_start
fit_end = pd.to_datetime('2098-12-31')

# subset to fit period and average over cms (!)
tclimy = climy.loc[(climy.index>fit_start) & (climy.index<fit_end)].groupby(level=1,axis=1).mean()
t= (tclimy.index.values - tclimy.index.values.min())/np.timedelta64(1,'D')

# --- subset over "historical" and "future" periods
cclimy = climy.loc[(climy.index.values >= cper_start) & (climy.index.values <= cper_end)].stack(level=0)
fclimy = climy.loc[(climy.index >=fper_start) & (climy.index <= fper_end)].stack(level=0)

# ------------- plot records ------------
# cm colors 
clist = ['darkgreen','red','purple','royalblue']

fig,axs=plt.subplots(4,1,sharex=True,figsize=(10,7))
# ---- total precip
ptot=climy.xs('ptot',1,1)
ptot.plot(style='+',ax=axs[0],color='black',lw=0.5,ms=3, alpha=0.5,legend=False)
ptot.mean(axis=1).rolling(window=10,center=True).mean().plot(ax=axs[0],color='tomato',ls='-', lw=2,legend=False)
safrany.ptot.plot(ax=axs[0],style='x',color='green')
axs[0].grid(which='both')
axs[0].set_xticklabels([])
axs[0].set_ylabel('Ptot [mm/y]')

# ---- total pet
pet=climy.xs('pet',1,1)
pet.plot(style='+',ax=axs[1],color='black',lw=0.5,ms=3, alpha=0.5,legend=False)
pet.mean(axis=1).rolling(window=10,center=True).mean().plot(ax=axs[1],color='tomato',ls='-', lw=2,legend=False)
safrany.pet.plot(ax=axs[1],style='x',color='green')
axs[1].grid(which='both')
axs[1].set_xticklabels([])
axs[1].set_ylabel('PET [mm/y]')

# ---- runoff
runoff=climy.xs('runoff',1,1)
runoff.plot(style='+',ax=axs[2],color='black',lw=0.5,ms=3, alpha=0.5,legend=False)
runoff.mean(axis=1).rolling(window=10,center=True).mean().plot(ax=axs[2],color='tomato',ls='-', lw=2,legend=False)
axs[2].grid(which='both')
axs[2].set_xticklabels([])
axs[2].set_ylabel('Runoff [mm/y]')

# ---- recharge
rech=climy.xs('rech',1,1)
rech.plot(style='+',ax=axs[3],color='black',lw=0.5,ms=3, alpha=0.5,legend=False)
rech.mean(axis=1).rolling(window=10,center=True).mean().plot(ax=axs[3],color='tomato',ls='-', lw=2,legend=False)
axs[3].grid(which='both')
axs[3].set_ylabel('Recharge [mm/y]')
axs[3].set_xlabel('')

for ax in axs:
    ax.axvline(pd.to_datetime(rcp_start),alpha=0.5,lw=1.5,ls=':',color='black')
    ax.axvspan(cper_start,cper_end,color='grey',alpha=0.3)
    ax.axvspan(fper_start,fper_end,color='grey',alpha=0.3)

lls =  [ Line2D([0], [0], label=f'Climate models', marker='+', linestyle='', color='black',alpha=0.5)]
lls += [Line2D([0], [0], label='Multi-model 10-year moving average',linestyle='-', color='tomato')]
lls += [Line2D([0], [0], label='SAFRAN', linestyle='',marker='x', color='green')]
fig.legend(handles=lls,loc='upper center',ncols=6,facecolor='white', framealpha=1)
fig.tight_layout()
fig.savefig(os.path.join('figs','long_term_records.pdf'),dpi=300)



# ---------- plot histograms and box plot for all scenarios ------------

fig, axs = plt.subplots(2,4, figsize=(9, 4), sharex='col', # Common x-axis
                       gridspec_kw={"height_ratios": (.7, .3)})

for i, col in enumerate(['ptot','pet','runoff','rech']):
    vmin = min(cclimy[col].min(),fclimy[col].min())
    vmax = max(cclimy[col].max(),fclimy[col].max())
    bins = np.linspace(vmin,vmax,10)
    # histograms 
    cclimy[col].hist(ax=axs[0,i],color='grey',bins=bins,grid=False,alpha=0.8,label='historical')
    fclimy[col].hist(ax=axs[0,i],color='darkred',bins=bins,grid=False,alpha=0.5, label='future')
    # boxplot 
    medianprops = dict(color='black')
    meanprops = dict(linestyle=None,marker='+', markeredgecolor='black', markerfacecolor='black')
    data = [fclimy[col].values, cclimy[col].values]
    flierprops = dict(marker='.', markersize=2,alpha=0.5)
    bplots = axs[1,i].boxplot(data,vert=False,
                            widths=0.4,
                            whis=(5,95),
                            showmeans=True,
                            medianprops=medianprops,
                            meanprops=meanprops,
                            patch_artist=True,
                            flierprops=flierprops
                            )
    for patch, color in zip(bplots['boxes'], ['darkred','grey']):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    axs[1,i].set_yticklabels(['',''])
    axs[1,i].set_xlabel(f'{col} [mm/y]')

axs[0,0].legend(loc='upper left')
axs[1,0].set_yticklabels(['fut','hist.'])
fig.tight_layout()

fig.savefig(os.path.join('figs','cf_allvars.pdf'),dpi=300)

# --- identify deciles of interest

qs = [0.05,0.5,0.95]

crech_qs = [np.quantile(cclimy.rech,q) for q in qs ]
frech_qs = [np.quantile(fclimy.rech,q) for q in qs ]

cqs_idx = [ cclimy.index[np.argmin(np.abs(cclimy.rech - rech))] for rech in crech_qs]
fqs_idx = [ fclimy.index[np.argmin(np.abs(fclimy.rech - rech))] for rech in frech_qs]

qs_labels = ['Q5','Q50','Q95']
cyrs  = pd.DataFrame({'rech' : crech_qs,
                      'cm'   : [ idx[1] for idx in cqs_idx],
                      'year' : [ idx[0].year for idx in cqs_idx]
                      }, index=qs_labels)
fyrs  = pd.DataFrame({'rech' : frech_qs,
                      'cm'   : [ idx[1] for idx in fqs_idx],
                      'year' : [ idx[0].year for idx in fqs_idx]
                      }, index=qs_labels)

simyrs = pd.concat([cyrs,fyrs],keys=['historical','future'])

simyrs.to_csv('simyrs.csv')


#---------------------------------------------------------------------
# ---------- plot scatter runoff against rech ------------
climy = pd.read_csv('climy.csv',header=[0,1],parse_dates=True)
fig, ax = plt.subplots(1,2, figsize=(8, 4)) # Common x-axis
axs[0].scatter(climy.loc[:,(slice(None),'rech')],climy.loc[:,(slice(None),'runoff')],
           marker='+',c='darkred')
axs.set_xlabel('Recharge [mm/y]')
axs[0].set_ylabel('Runoff [mm/y]')

axs[.scatter(climy.loc[:,(slice(None),'rech')],climy.loc[:,(slice(None),'srunoff')],
           marker='+',c='darkred')
axs[1].set_xlabel('Recharge [mm/y]')
axs[1].set_ylabel('Summer Runoff [mm/y]')
axs[1].set_title('Summer Runoff vs Annual Recharge')

fig.tight_layout()




