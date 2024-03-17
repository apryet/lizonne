import os, sys, glob
import shutil
import pandas as pd
import numpy as np
import pyemu
from pymarthe import MartheModel
from pymarthe.utils import marthe_utils, shp_utils, pest_utils, pp_utils
from pymarthe.mfield import MartheField, MartheFieldSeries
from pymarthe.helpers.postprocessing import PestPostProcessing
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
# plot settings
plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54
mdcol_width = 14/2.54
dbcol_width = 19/2.54

# ---------------------------------------------
# load pst and observation data from .config  
# ---------------------------------------------

pstfile = 'cal_lizonne.pst'
pst = pyemu.Pst(pstfile)
case = pst.filename.split('.')[0]

hdic, pdics, odics = pest_utils.read_config(pst.filename.replace('.pst','.config'))

ogdates = { odic['locnme']:pd.to_datetime(odic['dates_out'].split('|')) for odic in odics}

# iteration id for posterior distribution 
pt_id = pst.control_data.noptmax
pt_id = 3

# ---------------------------------------------
# load prior and last iteration observation ensembles  
# ---------------------------------------------

# prior ensemble of simulated values 
pr_oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(f"{case}.0.obs.csv"))
pr_oe._df.index=pr_oe.index.astype(str) #trick since index mixed types in oe.index (should find a fix)

# posterior (last iteration) ensemble of simulated values 
pt_oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(f"{case}.{pt_id}.obs.csv"))
pt_oe._df.index=pt_oe.index.astype(str) #trick since index mixed types in oe.index (should find a fix)

# observation with noise 
obswns = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(f"{case}.obs+noise.csv"))

# ---------------------------------------------
# load prior and last iteration parameter ensembles  
# ---------------------------------------------

# load parameter ensembles 
pe_dic = { it:pyemu.ParameterEnsemble.from_csv(
    pst=pst,filename=os.path.join(f"{case}.{it}.par.csv")
    ) for it in range(pt_id+1)}


# load prior parameter ensemble 
pr_pe = pe_dic[0]

# load posterior parameter ensemble 
pt_pe = pe_dic[pt_id]

# ---------------------------------------------
# save pe-oe for restart 
# ---------------------------------------------
'''
nreals_restart = 116

best_reals = pr_oe.phi_vector.sort_values().index[0:nreals_restart]

pr_oe_best = pr_oe.loc[best_reals].copy()
pr_oe_best.to_csv(filename=f'{case}.restart.obs.csv')

pr_pe_best = pr_pe.loc[best_reals.astype(str)].copy()
pr_pe_best.to_csv(filename=f'{case}.restart.par.csv')

obswns_best = obswns.loc[best_reals].copy()
obswns_best.to_csv(filename=f'{case}.restart.obs+noise.csv')


pst.pestpp_options['ies_restart_parameter_ensemble'] = f'{case}.restart.par.csv'
pst.pestpp_options['ies_observation_ensemble'] = f'{case}.restart.obs+noise.csv'
pst.pestpp_options['ies_restart_observation_ensemble'] = f'{case}.restart.obs.csv'
pst.pestpp_options['ies_bad_phi'] = 1e10 
pst.control_data.noptmax=2
pst.write(f'{case}.pst',version=2)
'''
# ---------------------------------------------
#  filtering with total phi values 
# ---------------------------------------------

# filter out bad phi
def filter_ens(pe,oe,qt,th):
    pv = oe.phi_vector
    th = min(pv.quantile(qt),th)
    keep_reals = pv.loc[pv<th].index
    return(pe.loc[keep_reals,:], oe.loc[keep_reals,:])

# remove extreme values from prior ensembles 
qt=1.0
th=1e8
fpr_pe, fpr_oe = filter_ens(pr_pe,pr_oe,qt=qt,th=th)

'''
# extract and save filtered posterior ensembles 
qt=0.5
th=1e8
fpt_pe, fpt_oe =filter_ens(pt_pe,pt_oe,qt=qt,th=th)
fpt_pe.to_csv(os.path.join(f"{case}.fpt.par.csv"))
'''
# ---------------------------------------------
# more granular filtering by phi contributions    
# ---------------------------------------------

def get_phi_df(oe):
    cols = oe._df.columns
    pst = oe.pst
    ogroups = pst.observation_data.loc[cols].groupby("obgnme").groups
    res = pd.DataFrame(data={'name':cols,
                            'group':pst.observation_data.loc[cols,'obgnme'].values,
                            'modelled':np.nan,
                            'residual':np.nan
                                })
    res.index = res.name
    obs = pst.observation_data.loc[cols, ['obsval','weight']]
    phi_rows = []
    for idx in oe._df.index.values:
        res.loc[cols,'modelled'] = oe._df.loc[idx, cols]
        res.loc[cols,'residual'] = res.loc[cols,'modelled'] - obs.loc[cols,'obsval']
        contribs =  pyemu.Pst.get_phi_components(ogroups,res,obs,None,None)
        phi_rows.append(contribs)
    return(pd.DataFrame(data=phi_rows, index=oe.index, columns=ogroups))


phi_df = get_phi_df(pt_oe)
phi_df = phi_df.loc[:,~phi_df.columns.str.endswith('mf')]

# take best phi values by quantiles 
phi_accept_by_ogroup = phi_df <  phi_df.quantile(0.75)
# observation groups of absolute head values 
hlocs = phi_df.columns[phi_df.columns.str.contains('x') & ~phi_df.columns.str.endswith('mf')]
# filter 
phi_accept = phi_accept_by_ogroup[hlocs].apply(lambda x : np.logical_and.reduce(x), axis=1)
keep_reals = phi_df.index[phi_accept]

# append base realization
keep_reals = pd.Index.union(keep_reals,["base"])

fpt_pe = pt_pe.loc[keep_reals,:]
fpt_oe = pt_oe.loc[keep_reals,:]

fpt_pe.to_csv('cal_lizonne.fpt.par.csv')


# ---------------------------------------------
# evolution of phi and hist for last iteration
# ---------------------------------------------
fig, axes = plt.subplots(1, 2, sharey=True, figsize=(10,3.5))
# left
ax = axes[0]
phi = pd.read_csv(os.path.join(f"{case}.phi.actual.csv"),index_col=0)
phi.index = phi.total_runs
phi.iloc[:,6:].apply(np.log10).plot(legend=False,lw=0.5,color='k', ax=ax)
ax.set_title(r'Actual ')
ax.set_ylabel(r'log ')
# right
ax = axes[-1]
phi = pd.read_csv(os.path.join(f"{case}.phi.meas.csv"),index_col=0)
phi.index = phi.total_runs
phi.iloc[:,6:].apply(np.log10).plot(legend=False,lw=0.2,color='r', ax=ax)
ax.set_title(r'Measured+Noise ')
fig.tight_layout()
fig.savefig(os.path.join('pproc','phi_evol.pdf'),dpi=300)

# ---------------------------------------------
# compare phi distrib for prior and last iter
# ---------------------------------------------

pr_logphi = pr_oe.phi_vector.apply(np.log10)
pt_logphi = pt_oe.phi_vector.apply(np.log10)

bins=np.histogram(np.hstack((pr_logphi)), bins=40)[1] #get the bin edges

bins=np.histogram(np.hstack((pr_logphi,pt_logphi)), bins=40)[1] #get the bin edges

fig,ax = plt.subplots(1,1,figsize=(3.5,3.5))
pr_logphi.hist(bins=bins,ax=ax,fc="0.5",ec="none",alpha=0.5,density=False,label='prior')
pt_logphi.hist(bins=bins,ax=ax,fc="b",ec="none",alpha=0.5,density=False,label='posterior')
ax.legend()
_ = ax.set_xlabel('log10($\Phi$)')
fig.tight_layout()

fig.savefig(os.path.join('pproc','phi.pdf'),dpi=300)

# ---------------------------------------------
#  time series   
# ---------------------------------------------

# get ensemble of time series for given observation group 
def get_og_ts(oe,onames,odates, trans):
    ts = oe._df.loc[:,onames].T.apply(trans)
    ts.index = odates 
    return(ts)

def plot_tseries_ensembles(pr_oe, pt_oe, obswns, ognmes, ogdates, trans=None, ylabel='',legend=True ):
    # get the observation data from the control file and select 
    obs = pst.observation_data.copy()
    fig,axes = plt.subplots(len(ognmes),1,sharex=True,figsize=(dbcol_width,0.8*dbcol_width))
    if trans==None:
        trans=[lambda x : x]*len(ognmes)
    if not type(trans)==list:
        trans = [trans]*len(ognmes)
    # for each observation group (i.e. timeseries)
    for ax,og,t in zip(axes,ognmes,trans):
        # get values 
        oobs = obs.loc[obs.obgnme==og.lower(),:].copy()
        onames = oobs.obsnme.values
        odates = ogdates[og]
        # plot prior
        if pr_oe is not None :
            ts = get_og_ts(pr_oe,onames,odates, trans=t)
            ts.plot(ax=ax,color='grey',lw=0.5,alpha=0.10,legend=False)
        # plot posterior
        if pt_oe is not None :
            ts = get_og_ts(pt_oe,onames,odates,trans=t)
            ts.plot(ax=ax,color='red',lw=0.5,alpha=0.50,legend=False)
            ts['base'].plot(ax=ax,color='green',alpha=1,lw=2,legend=False)
        # plot measured+noise 
        if obswns is not None :
            ts = get_og_ts(obswns,onames,odates,trans=t)
            ts.plot(ax=ax,color='blue',lw=0.5,alpha=0.05,legend=False)
        # plot obs
        ax.plot(odates, oobs.obsval.apply(t).values,color="black",alpha=1,lw=1)
        ax.set_title(og,loc="left")
        ax.set_ylabel(ylabel)
        lpr = Line2D([0], [0], label='Sim. prior', color='grey')
        lpt = Line2D([0], [0], label='Sim. posterior', color='red')
        lbase = Line2D([0], [0], label='Sim. base', color='green')
        lobs = Line2D([0], [0], label='Observed', color='black')
        lobsn = Line2D([0], [0], label='Obs.+noise', color='blue')
        if legend:
            ax.legend(handles=[lpr,lpt,lbase,lobs,lobsn],loc='upper left',ncols=5)
            plot_legend=False
    fig.tight_layout()
    return fig

# gaging stations
gstations = ['P7250001','P7270001','P8215010','P8284010']

fig = plot_tseries_ensembles(pr_oe, pt_oe, obswns , gstations, ogdates,trans=lambda x : 10**x, ylabel='River discharge [m$^3$/s]')
fig.savefig(os.path.join('pproc','pr_pt_qsimobs.png'),dpi=300)

fig = plot_tseries_ensembles(None, pt_oe, obswns , gstations, ogdates,trans=lambda x : 10**x, ylabel='River discharge [m$^3$/s]')
fig.savefig(os.path.join('pproc','pt_qsimobs.png'),dpi=300)

fig = plot_tseries_ensembles(None, fpt_oe, obswns , gstations, ogdates,trans=lambda x : 10**x, ylabel='River discharge [m$^3$/s]')
fig.savefig(os.path.join('pproc','fpt_qsimobs.png'),dpi=300)

# selection of observation wells 
obswells = ['07333X0027','07345X0023','07346X0017','07346X0083','07574X0014']

fig = plot_tseries_ensembles(pr_oe, pt_oe, obswns , obswells, ogdates,trans=lambda x : x, ylabel='Water level [m NGF]')
fig.savefig(os.path.join('pproc','pr_pt_hsimobs.png'),dpi=300)

fig = plot_tseries_ensembles(None, pt_oe, obswns , obswells, ogdates,trans=lambda x : x, ylabel='Water level [m NGF]')
fig.savefig(os.path.join('pproc','pt_hsimobs.png'),dpi=300)

fig = plot_tseries_ensembles(None, fpt_oe, obswns , obswells, ogdates,trans=lambda x : x, ylabel='Water level [m NGF]')
fig.savefig(os.path.join('pproc','fpt_hsimobs.png'),dpi=300)

#   ---------  selection of the two (article)

obslocs = pd.read_excel(os.path.join('..','data','SIG','obsloc_labels.xlsx'),index_col='id')
locs = ['P8284010','P7250001','07333X0027','07346X0017']
labels = obslocs.loc[locs].label
units = ['m$^3$/s']*2 + ['m a.s.l.']*2
trans = [lambda x : 10**x,lambda x : 10**x,lambda x : x,lambda x : x]

fig = plot_tseries_ensembles(None, fpt_oe, obswns , locs, ogdates,trans=trans, ylabel='',legend=False)

cal_start = pd.to_datetime('2012-08-01')
cal_end =pd.to_datetime('2019-07-31')

fig.axes[0].set_xlim(cal_start,cal_end)
fig.axes[0].set_xticks([])
fig.axes[0].set_xticks([],minor=True)
fig.axes[0].xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonth=1,bymonthday=1))
fig.axes[0].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m'))
fig.axes[0].xaxis.set_minor_locator(matplotlib.dates.MonthLocator())
fig.axes[-1].tick_params(axis='x', which='major', pad=0)

for ax,unit,label in zip(fig.axes,units,labels):
    ax.set_ylabel(unit)
    ax.set_title(f' {label}',loc="left",y=0.80)


fig.axes[0].set_ylim(0,45)
fig.axes[1].set_ylim(0,5)
fig.axes[2].set_ylim(94,135)
#fig.axes[3].set_ylim(100,175)

lpt = Line2D([0], [0], label='Posterior', color='red')
lbase = Line2D([0], [0], label='Center', color='green')
lobs = Line2D([0], [0], label='Obs.', color='black')
lobsn = Line2D([0], [0], label='Obs.+noise', color='blue')
fig.axes[0].legend(handles=[lpt,lbase,lobs,lobsn],loc='upper right',ncols=4,fontsize=11)
fig.axes[0].margins(x=0)
fig.align_ylabels()
fig.tight_layout()
fig.subplots_adjust(hspace=0.10)

fig.savefig(os.path.join('pproc','fpt_q_and_h_simobs.png'),dpi=150)
fig.savefig(os.path.join('pproc','fpt_q_and_h_simobs.pdf'),dpi=150)

# ---------------------------------------------
# final ensemble selection      
# ---------------------------------------------

# see Shafii et al 2015 : http://dx.doi.org/10.1016/j.jhydrol.2015.01.051
# compute reliability
def get_reliability(pt_oe, obs):
    pt_oe = pt_oe
    lqt = pt_oe._df.quantile(qt_lb)
    uqt = pt_oe._df.quantile(qt_ub)
    obs_in_bds = obs.between(lqt,uqt)
    # re= (number of obs in pred. interval)/(total number of obs).
    r = obs_in_bds.sum()/obs_in_bds.shape[0]
    return(r)

# compute (normalized) sharpness 
def get_sharpness(pt_oe,pr_oe=None):
    pt_wdth = pt_oe._df.quantile(qt_ub)-pt_oe._df.quantile(qt_lb)
    # When pr_oe is not provided, return absolute sharpness value
    # (not relevant for heterogeneous observational dataset).
    if pr_oe is None:
        s = pt_wdth.mean()
    # When pr_oe is provided, return normalized sharpness indicator
    # s=1 for single value
    # s=0 for prediction interval equal to posterior interval
    # see 
    else :
        pr_wdth = pr_oe._df.quantile(qt_ub)-pr_oe._df.quantile(qt_lb)
        wdth_ratio = 1 - pt_wdth/pr_wdth
        s = wdth_ratio.mean()
    return(s)

qt_lb=0.05
qt_ub=0.95

# -----------------------------------------------------
# --- reliability and sharpness analysis   
# -----------------------------------------------------

# --- plot reliability and sharpness through iterations  
fig,ax=plt.subplots(1,1,figsize=(6,4),sharey=True)
nit = phi.shape[0]

# all obs 
rvals,svals,dvals=[],[],[]
for i in range(nit):
    # load ensemble of simulated observations at IES iteration i
    ioe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(f"{case}.{i}.obs.csv"))
    # compute sharpness and reliability considering all subset
    r=get_reliability(ioe,obswns.loc['base'])
    s= get_sharpness(ioe,pr_oe)
    d = np.sqrt((1-r)**2+(1-s)**2)
    rvals.append(r)
    svals.append(s)
    dvals.append(d)

rsd_iter = rvals, svals, dvals

ax.plot(rvals,'--+', color='royalblue',label='Reliability')
ax.plot(svals,'--+', color='darkgreen',label='Sharpness')
ax.plot(dvals,'-+', color='black', label='Distance to optimum')
ax.axvline(3,color='grey')
ax.set_xlabel('IES iterations')
ax.legend('lower right')
fig.savefig(os.path.join('pproc','rsd_through_iterations.pdf'),dpi=300)

# --- plot reliability and sharpness through final filtered ensemble size   
fig,ax=plt.subplots(1,1,figsize=(6,4),sharey=True)

# posterior observation ensemble sorted by decreasing phi
spt_oe = fpt_oe.loc[fpt_oe.phi_vector.sort_values(ascending=False).index]

# compute distance to optimum (r=1,s=1) to get optimum ffpt ensemble size

# all obs. 
rvals,svals,dvals=[],[],[]
for n in range(spt_oe.shape[0]):
    r=get_reliability(spt_oe.iloc[:n+1],obswns.loc['base'])
    s= get_sharpness(spt_oe.iloc[:n+1],pr_oe)
    d = np.sqrt((1-r)**2+(1-s)**2)
    rvals.append(r)
    svals.append(s)
    dvals.append(d)

ax.plot(rvals,'--+', color='royalblue',label='Reliability')
ax.plot(svals,'--+', color='darkgreen',label='Sharpness')
ax.plot(dvals,'-+', color='black', label='Distance to optimum')
ax.set_xlabel('Number of realizations in the final filtered posterior ensemble')
ax.legend('lower right')
fig.tight_layout()
fig.savefig(os.path.join('pproc','rsd_through_ensemble_size.pdf'),dpi=300)

# optimum size for the final filtered posterior ensemble
nreals_ffpt_ens=np.argmin(dvals)
print(f'Optimum number of best realizations is: {nreals_ffpt_ens}')


# final ensemble selection 
ffpt_oe = spt_oe.iloc[:nreals_ffpt_ens] # set to nreals_ffpt_ens after examination of rsd.pdf
ffpt_ids = ffpt_oe.index
ffpt_pe = pt_pe.loc[ffpt_ids]

# identify new base (center) obs realization from normalized distance to mean, and update index 
center_real = (ffpt_oe - ffpt_oe.mean()).div(ffpt_oe.mean()).apply(np.linalg.norm, axis=1).idxmin()
ffpt_pe._df.index = ffpt_pe.index.str.replace(center_real,'ffpt_center')
ffpt_oe._df.index = ffpt_oe.index.str.replace(center_real,'ffpt_center')

# write final parameter ensemble 
ffpt_pe.to_csv('cal_lizonne.ffpt.par.csv')

# plot final filtered ensemble 
fig = plot_tseries_ensembles(None, ffpt_oe, obswns , gstations, ogdates,trans=lambda x : 10**x, ylabel='River discharge [m$^3$/s]')
fig.savefig(os.path.join('pproc','ffpt_qsimobs.png'),dpi=300)

fig = plot_tseries_ensembles(None, ffpt_oe, obswns , obswells, ogdates,trans=lambda x : x, ylabel='Water level [m NGF]')
fig.savefig(os.path.join('pproc','ffpt_hsimobs.png'),dpi=300)


# ---------------------------------------------
# evolution of phi and hist for last iteration + rsd 
# ---------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(dbcol_width,0.4*dbcol_width),width_ratios=[6.5,4])
# phi evol
phi_color = 'grey'
ax = axes[0]
phi = pd.read_csv(os.path.join(f"{case}.phi.actual.csv"),index_col=0)
phi.iloc[:,6:].apply(np.log10).plot(legend=False,lw=0.5,color=phi_color,alpha=0.8, ax=ax)

ax.set_ylabel('log10($\Phi$)')
ax.set_xlabel('IES iterations')
ax.axvline(0,color='grey',ls='--',lw=3,alpha=0.6)
ax.axvline(3,color='darkred',ls='--',lw=3,alpha=0.6)
ax.text(0.5,7.4,'a)',fontsize=20)

phi_evol_handle = Line2D([0], [0], label='$\Phi$ (232 realizations)', color=phi_color, marker=None, linestyle= '-')
ax.legend(handles=[phi_evol_handle], alignment='left', loc='lower left',fontsize=9,  bbox_to_anchor=(0.08, 0))

# rsd 
twax= ax.twinx()
rsd_color = 'blue'
rvals, svals, dvals = rsd_iter
rl = twax.plot(rvals,'--+', color='royalblue',label='Reliability')
sl = twax.plot(svals,'--*', color='darkgreen',label='Sharpness')
dl = twax.plot(dvals,'-', color='black', label='Distance/opt.')
twax.set_ylabel('Performance ratios [-]')
# legend
lns = rl + sl + dl
labs = [l.get_label() for l in lns]
twax.legend(lns, labs,loc='upper right', bbox_to_anchor=(1.01, 1.01),alignment='left',fontsize=9)

ax.set_ylabel('log10($\Phi$)')

# histogram
ax = axes[1]

bins=np.histogram(np.hstack((pr_logphi,pt_logphi)), bins=40)[1] #get the bin edges
pt_logphi.hist(bins=bins,ax=ax,fc="darkred",ec="none",alpha=0.6,density=False,label='Posterior (It. 3)')
pr_logphi.hist(bins=bins,ax=ax,fc="0.5",ec="none",alpha=0.6,density=False,label='Prior (It. 0)')
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.legend(loc='upper right')
ax.set_ylabel('Frequency')
_ = ax.set_xlabel('log10($\Phi$)')

ax.text(5.05,42,'b)',fontsize=20)
fig.subplots_adjust(wspace=0.30)

fig.tight_layout()

fig.savefig(os.path.join('pproc','phi_rsd.pdf'),dpi=300)


# ---------------------------------------------
# fit statistics   
# ---------------------------------------------

# get all metrics 
m = pyemu.utils.metrics.calc_metric_ensemble(ffpt_oe._df,pst)

# drop fluctuations
m = m.loc[:,~m.columns.str.endswith('mf')]


# multiindex
locs = [cname.split('_')[-1] for cname in m.columns]
metrics  = ['_'.join(cname.split('_')[:-1]) for cname in m.columns]
m.columns = pd.MultiIndex.from_frame(
        pd.DataFrame(
            {'metric':metrics,
            'loc':locs}
            )
        )


# label dics 
label_dic = {id.lower():label for id,label in zip(obslocs.index,obslocs.label)}
m.columns = m.columns.set_levels(m.columns.levels[1].map(label_dic,na_action='ignore'),level=1)

allgstations = obslocs.loc[obslocs.index.str.startswith('P'),'label'].sort_values()
allobswells = obslocs.loc[obslocs.index.str.contains('X'),'label'].sort_values()

fig, axs = plt.subplots(1,2,figsize=(dbcol_width,0.33*dbcol_width),width_ratios=[len(allgstations),len(allobswells)])
ax1 = m.loc[:,('KGE',allgstations)].boxplot(showfliers=False,ax=axs[0])
ax1.set_xticklabels(allgstations.values)
ax1.set_ylabel('KGE [-]')
ax1.set_title('Gaging stations (KGE)')
ax2 = m.loc[:,('RMSE',allobswells)].boxplot(showfliers=False, ax=axs[1])
ax2.set_xticklabels(allobswells.values)
ax2.set_title('Observation wells (RMSE)')
ax2.set_ylabel('RMSE [m]')
fig.tight_layout()
fig.savefig(os.path.join('pproc','metrics_boxplots.pdf'),dpi=300)

# ---------------------------------------------
# parameters evolution  
# ---------------------------------------------

dic_props = dict(color = 'k', lw =  0.5)
meanprops = dict(marker = 'o', mfc = 'none', ms = 3, lw = 0.5, mec = 'k')
flierprops = dict(marker = 'o', mfc = 'none', ms = 3, lw = 0.1, mec = 'k')

def plot_multviolins(ens_list, ax, showpoints=False):
    bplot = ax.boxplot(ens_list, widths = 0.4, showbox = True, showcaps = False, 
                showmeans = True, showfliers = False, boxprops = dic_props,
                whiskerprops = dic_props, capprops = dic_props, medianprops = dic_props, 
                meanprops = meanprops, flierprops = flierprops)
    # violin plot (white)
    vp = ax.violinplot(ens_list, widths = 0.8, points = 100, showmeans = False, 
                       showmedians = False, showextrema = False)
    for pc in vp['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('white')
        pc.set_alpha(1)
    # violin plot (color)
    vp = ax.violinplot(ens_list, widths = 0.8, points = 100, showmeans = False,
                       showmedians = False, showextrema = False)
    for i, pc in enumerate(vp['bodies']):
        pc.set_facecolor('blue')
        pc.set_edgecolor('blue')
        pc.set_alpha(0.4)
    if showpoints != False: # Points
        for i, tick in enumerate(xticks_list):
            y = ens_list[i]
            x = np.random.normal(tick, 0.04, size = len(y))
            ax.plot(x, y, 'k.', alpha = 0.1)
    ax.set_xticklabels([])
    ax.set_axisbelow(True)
    ax.yaxis.grid()


par = pst.parameter_data

df_list = []

for i in range(pt_id+1):
    pe_df = pe_dic[i]._df.copy(deep=True)
    pe_df.columns = pd.MultiIndex.from_frame(
    pd.DataFrame({
        'pargp':par.loc[pe_df.columns,'pargp'].values,
        'layer':pe_df.columns.str.extract(r"_l(\d+)_z",expand=False).astype(float)-1,
        'parnme':pe_df.columns.values
        }))
    df_list.append(pe_df)

pes = pd.concat(df_list,keys=range(pt_id+1))

pfs =['emmca','emmli','permh']

for prefix in pfs:
    fig, axs = plt.subplots(6,1,figsize=(10,18))
    for l in range(6):
        parvals_list=[pes.loc[(i,slice(None)),\
                (pes.columns.get_level_values(0).str.startswith(prefix),l,slice(None))
                              ].values.ravel() \
                for i in range(pt_id+1)]
        ax = plot_multviolins(parvals_list, axs[l])
        axs[l].set_title(f'{prefix} for layer {l+1}')
    fig.tight_layout()
    fig.savefig(os.path.join('pproc',f'evol_{prefix}.pdf'),dpi=300)
    plt.close(fig)


pfs =['cap_sol_progr', 'equ_ruis_perc', 't_demi_percol','perm_r_zpc']
fig, axs = plt.subplots(len(pfs),1,figsize=(10,18))
for i,prefix in enumerate(pfs):
    parvals_list=[pes.loc[(i,slice(None)),\
            (pes.columns.get_level_values(0).str.startswith(prefix),slice(None),slice(None))
                          ].values.ravel() \
            for i in range(pt_id+1)]
    ax = plot_multviolins(parvals_list, axs[i])
    axs[i].set_title(f'{prefix}')

fig.tight_layout()
fig.savefig(os.path.join('pproc',f'evol_other_params.pdf'),dpi=300)
plt.close(fig)


# ---------------------------------------------
# parameters pt - pr distributions 
# ---------------------------------------------

# get parameter df and append phi value value 
pr_pe_df = pr_pe._df.copy(deep=True)
pt_pe_df = pt_pe._df.copy(deep=True)

# set multiindexed columns with parmameter groups (for hist grouping)
par = pst.parameter_data
parnmes= pr_pe_df.columns
colmix = pd.MultiIndex.from_frame(
        pd.DataFrame({
            'pargp':par.loc[parnmes,'pargp'],
            'parnme':parnmes
            }))

pr_pe_df.columns = colmix
pt_pe_df.columns = colmix


# plot hydraulic properties (per layer)
pfs =['emmca','emmli','permh','perm_r']
for prefix in pfs:
    pr_df = pr_pe_df.loc[:,(pr_pe_df.columns.get_level_values(0).str.startswith(prefix),slice(None))].melt()
    pt_df = pt_pe_df.loc[:,(pt_pe_df.columns.get_level_values(0).str.startswith(prefix),slice(None))].melt()
    axs = pr_df.hist(by='pargp',fc="0.5",ec="none",alpha=0.5,density=False,label='prior')
    axs = pt_df.hist(ax=axs,by='pargp',fc="b",ec="none",alpha=0.5,density=False,label='posterior')
    try : 
        fig = axs[0].get_figure()
    except : 
        fig = axs.get_figure()
    fig.tight_layout()
    fig.savefig(os.path.join('pproc',f'pr_pt_hist_{prefix}.pdf'),dpi=300)

# plot soil parameters 
pfs =['cap_sol_progr', 'equ_ruis_perc', 't_demi_percol']
fig,axs=plt.subplots(1,3,figsize=(12,6))
for i,prefix in enumerate(pfs):
    pr_df = pr_pe_df.loc[:,(pr_pe_df.columns.get_level_values(0).str.startswith(prefix),slice(None))].melt()
    pt_df = pt_pe_df.loc[:,(pt_pe_df.columns.get_level_values(0).str.startswith(prefix),slice(None))].melt()
    ax = pr_df.hist(ax=axs[i],by='pargp',fc="0.5",ec="none",alpha=0.5,density=False,label='prior')
    ax = pt_df.hist(ax=axs[i],by='pargp',fc="b",ec="none",alpha=0.5,density=False,label='posterior')
    fig.tight_layout()

fig.savefig(os.path.join('pproc','pr_pt_hist_surf.pdf'),dpi=300)


# ---------------------------------------------
# extract realization
# ---------------------------------------------
# best real
real = pt_oe.index[pt_oe.phi_vector.argmin()]
# user defined real
par_set = pt_pe.loc[str(real)]
pst.parameter_data.loc[par_set.index,'parval1']=par_set.values
pst.control_data.noptmax=0
pst.write(f'caleval_r{real}_it{pt_id}_best.pst')

# write pst with posterior base pe
real = 'ffpt_center'
par_set = ffpt_pe.loc[real]
pst.parameter_data.loc[par_set.index,'parval1']=par_set.values
pst.control_data.noptmax=0
pst.write('caleval_ffpt_center.pst')


'''
# run pest with center real from final filtered posterior ensemble (ffpt)
pyemu.os_utils.run('pestpp-ies caleval_ffpt_center.pst')

# copy model files with center filtered posterior real to sim dir
tpl_sim_dir = os.path.join('..','tpl_sim')

# post-proc center real
print('post processing center realization...'
exec(open('../pproc_cal.py').read())

print('Copying model files of center real to {tpl_sim_dir}')
for f in glob.glob("Lizonne.*"): shutil.copy(f, tpl_sim_dir)

# copying the calibration heads (they will be used to define the constraints on gw level)
shutil.copy('chasim.out',os.path.join(tpl_sim_dir,'chasim.out.cal'))

# copying the final filtered posterior parameter ensemble (it will constitute the par stack for MOU)
shutil.copy('cal_lizonne.ffpt.par.csv',os.path.join(tpl_sim_dir,'cal_lizonne.ffpt.par.csv'))

'''
