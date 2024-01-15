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
    # annual aggregation (remove incomplete 1st and last years)
    climy = clim.groupby(pd.Grouper(freq='Y')).sum().iloc[1:-1,:]
    # add to compil dic 
    climy_dic[i]=climy

# concat to single df 
climy = pd.concat(climy_dic.values(),keys=climy_dic.keys(),axis=1)

# subset to fit trend line
rcp_start= pd.to_datetime('2005-08-01') # rcp start 
fit_start = rcp_start
fit_end = pd.to_datetime('2098-12-31')

# subset to fit period and average over 4 cms (!)
tclimy = climy.loc[(climy.index>fit_start) & (climy.index<fit_end)].groupby(level=1,axis=1).mean()
t= (tclimy.index.values - tclimy.index.values.min())/np.timedelta64(1,'D')

# --- subset over "historical" and "future" periods
cclimy = climy.loc[(climy.index.values >= cper_start) & (climy.index.values <= cper_end)].stack(level=0)
fclimy = climy.loc[(climy.index >=fper_start) & (climy.index <= fper_end)].stack(level=0)

# ------------- plot records ------------
# cm colors 
clist = ['darkgreen','red','purple','royalblue']

fig,axs=plt.subplots(4,1,sharex=True,figsize=(12,7))
# ---- total precip
ptot=climy.xs('ptot',1,1)
ptot.plot(style='+',ax=axs[0],color=clist, alpha=0.5,legend=False)
ptot.rolling(window=10,center=True).mean().plot(ax=axs[0],color=clist,alpha=0.75,legend=False)
safrany.ptot.plot(ax=axs[0],style='x',color='orange')
safrany.ptot.rolling(window=10,center=True).mean().plot(ax=axs[0],color='orange',alpha=0.75,legend=False)
# fit trend line 
pfit = np.polyfit(t,tclimy.ptot.values, 1)
pd.DataFrame({'reg':np.polyval(pfit,t)},index=tclimy.index).plot(ax=axs[0],color='black',ls='--',legend=False)
axs[0].set_xticklabels([])
axs[0].set_ylabel('Ptot [mm/y]')

# ---- total pet
pet=climy.xs('pet',1,1)
pet.plot(style='+',ax=axs[1],color=clist, alpha=0.5,legend=False)
pet.rolling(window=10,center=True).mean().plot(ax=axs[1],color=clist,alpha=0.75,legend=False)
safrany.pet.plot(ax=axs[1],style='x',color='orange')
safrany.pet.rolling(window=10,center=True).mean().plot(ax=axs[1],color='orange',alpha=0.75,legend=False)
# fit trend line 
pfit = np.polyfit(t,tclimy.pet.values, 1)
pd.DataFrame({'reg':np.polyval(pfit,t)},index=tclimy.index).plot(ax=axs[1],color='black',ls='--',legend=False)
axs[1].set_xticklabels([])
axs[1].set_ylabel('ETP [mm/y]')

# ---- runoff
runoff=climy.xs('runoff',1,1)
runoff.plot(style='+',ax=axs[2],color=clist, alpha=0.5,legend=False)
runoff.rolling(window=10,center=True).mean().plot(ax=axs[2],color=clist,alpha=0.75,legend=False)
# fit trend line 
pfit = np.polyfit(t,tclimy.runoff.values, 1)
pd.DataFrame({'reg':np.polyval(pfit,t)},index=tclimy.index).plot(ax=axs[2],color='black',lw=1.5,ls='--',legend=False)
axs[2].set_ylabel('Runoff [mm/y]')

# ---- recharge
rech=climy.xs('rech',1,1)
rech.plot(style='+',ax=axs[3],color=clist, alpha=0.5,legend=False)
rech.rolling(window=10,center=True).mean().plot(ax=axs[3],color=clist,alpha=0.75,legend=False)
# fit trend line 
pfit = np.polyfit(t,tclimy.rech.values, 1)
pd.DataFrame({'reg':np.polyval(pfit,t)},index=tclimy.index).plot(ax=axs[3],color='black',lw=1.5,ls='--',legend=False)
axs[3].set_ylabel('Recharge [mm/y]')

for ax in axs:
    ax.axvline(pd.to_datetime(rcp_start),alpha=0.5,ls=':',color='black')
    ax.axvspan(cper_start,cper_end,color='grey',alpha=0.3)
    ax.axvspan(fper_start,fper_end,color='grey',alpha=0.3)

lls =  [ Line2D([0], [0], label=f'DRIAS CM.{i+1}', marker='+', color=c) for i,c in enumerate(clist)]
lls += [Line2D([0], [0], label='SAFRAN', linestyle=None,marker='x', color='orange')]
lls += [Line2D([0], [0], label='trend',linestyle='--', color='black')]
fig.legend(handles=lls,loc='upper center',ncols=6)

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
    bplots = axs[1,i].boxplot(data,vert=False,
                            widths=0.4,
                            whis=(0.05,95),
                            showmeans=True,
                            medianprops=medianprops,
                            meanprops=meanprops,
                            patch_artist=True
                            )
    for patch, color in zip(bplots['boxes'], ['darkred','grey']):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    axs[1,i].set_yticklabels(['',''])
    axs[1,i].set_xlabel(f'{col} [mm/y]')

axs[0,0].legend(loc='upper left')
axs[1,0].set_yticklabels(['fut','cur.'])
fig.tight_layout()

fig.savefig(os.path.join('figs','cf_allvars.pdf'),dpi=300)


# --- plot histogram and box plots of recharge for each scenario ------------

fig, axs = plt.subplots(2,4, figsize=(9, 4), sharex=True, # Common x-axis
                       gridspec_kw={"height_ratios": (.7, .3)})
col='rech'

for i,cm_id in enumerate(cm_df.index):
    # subset to periods and single climate models (cm_id)
    scclimy = climy.loc[(climy.index.values>cper_start) & (climy.index.values< cper_end),(cm_id,slice(None))].stack(level=0)
    sfclimy = climy.loc[(climy.index>fper_start) & (climy.index< fper_end),(cm_id,slice(None))].stack(level=0)
    # histogram
    scclimy[col].hist(ax=axs[0,i], bins=bins, color='grey', grid=False,alpha=0.8,label='historical')
    sfclimy[col].hist(ax=axs[0,i], bins=bins, color='darkred', grid=False,alpha=0.8,label='future')
    axs[0,i].set_title(f'CM {cm_id}')
    # boxplot 
    medianprops = dict(color='black')
    meanprops = dict(linestyle=None,marker='+', markeredgecolor='black', markerfacecolor='black')
    data = [sfclimy[col].values, scclimy[col].values]
    bplots = axs[1,i].boxplot(data,vert=False,
                            widths=0.4,
                            whis=(0.05,95),
                            showmeans=True,
                            medianprops=medianprops,
                            meanprops=meanprops,
                            patch_artist=True
                            )
    axs[1,i].set_yticklabels(['fut','cur.'])
    # fill with colors
    colors = ['darkred','grey']
    for patch, color in zip(bplots['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    axs[1,i].set_xlabel('Recharge [mm/y]')

axs[0,0].legend(loc='upper left')
fig.tight_layout()

fig.savefig(os.path.join('figs',f'cf_rech_cms.pdf'),dpi=300)

# --- identify deciles of interest

qs = [0.05,0.5,0.95]

crech_qs = [np.quantile(cclimy.rech,q) for q in qs ]
frech_qs = [np.quantile(fclimy.rech,q) for q in qs ]

cqs_idx = [ cclimy.index[np.argmin(np.abs(cclimy.rech - rech))] for rech in crech_qs]
fqs_idx = [ fclimy.index[np.argmin(np.abs(fclimy.rech - rech))] for rech in frech_qs]

qs_labels = ['Q10','Q50','Q90']
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

