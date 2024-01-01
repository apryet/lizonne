import os, sys
import pandas as pd
import numpy as np
from pymarthe import MartheField,MartheModel
from matplotlib import pyplot as plt
import scipy

'''
This script smooths the geometry by setting cells as active / inactive in the permh field.
1) removes all cells below thk_thd
2) smooth the geometry given kernel and ncell
3) removes all cells below thk_thd in smoothed geometry (iterative approach)

Geometry is not modified, cells are just activated/deactivated in the hsubs and permh files. 

Current version is only valid for regular grids.

'''
# ---------------------------------
#  -------- settings --------------
# ---------------------------------

# thickness below which cells are discarded
thk_thd_1 = 5. # m 
thk_thd_2 = 2. # m 

# 2d 3x3 filter kernel
kernel = np.ones((3,3))

# number of active cells in the kernel to consider the cell as active 
ncell = 1


# ---------------------------------
# ----- crop model domain ---------
# ---------------------------------

'''
permh = mm.prop['permh']
permh_data = permh.data.copy()
permh_data[(permh.data.x>510000) & (permh.data.layer==0)]=0
permh_data[(permh.data.x>514500) & (permh.data.layer==4)]=0
permh_data[(permh.data.x>514500) & (permh.data.layer==5)]=0
permh.set_data(permh_data)
permh.write_data() # NOTE should reload imask when altering permh !
'''

# ---------------------------------
#  ------ pre-processing ----------
# ---------------------------------

# load model (just for geometry data)
mm = MartheModel('Lizonne.rma')

# load original permh field
permh_mf = MartheField('permh','Lizonne.permh.bckup', mm, use_imask=False)
permh_3d = permh_mf.as_3darray()

# load original topo, husbs
topo = MartheField('topog','Lizonne.topog.bckup', mm,use_imask=False).as_3darray()[0,:,:]
hsubs_3d = MartheField('hsubs','Lizonne.hsubs.bckup', mm,use_imask=False).as_3darray()

# set no data values to nan
permh_3d[permh_3d==0]=np.nan
hsubs_3d[hsubs_3d==9999]=np.nan
topo[topo==9999]=np.nan

# ---------------------------------------------
#  ------ remove cells below thk_thd  ---------
# ---------------------------------------------

# arrays of interface elevation
itfs_3d = np.concatenate((topo[None,:,:],hsubs_3d),axis=0)

# arrays of thickness (accounting for gaps)
mask = np.isnan(itfs_3d)
idx = np.where(~mask,(np.ones(mask.shape[::-1])*np.arange(mask.shape[0])).T,0)
idx = np.maximum.accumulate(idx, axis=0,dtype=int)
itfs_3d[mask] = itfs_3d[ idx[mask], np.nonzero(mask)[1], np.nonzero(mask)[2]]
thk_3d = -1*np.diff(itfs_3d,axis=0)

# save a copy of original fields (for comparative purpose)
permh_3d_org = permh_3d.copy()
hsubs_3d_org = hsubs_3d.copy()
thk_3d_org = thk_3d.copy()

# remove cells below thickness threshold
permh_3d[thk_3d<thk_thd_1] = np.nan

# infer original imask 
imask_3d = np.ones_like(permh_3d,dtype=int)
imask_3d[np.isnan(permh_3d)] = 0

# list of 2d arrays
newmask_list, newpermh_list, newhsubs_list = [], [], []

# initialize array of nan (will be filled)
nans = np.empty(topo.shape)
nans[:]=np.nan

# initialize coordinate arrays for interpolation 
x, y = np.arange(0, topo.shape[1]), np.arange(0, topo.shape[0])
xx, yy = np.meshgrid(x, y)

# ----------------------------------------------
# ----- smooth imask to remove isolated patches
# ----------------------------------------------

for l in range(mm.nlay):
    print(f'Smoothing layer {l+1:02}...')
    # original imask, permh, hsubs (2d)
    imask = imask_3d[l,:,:]
    permh = permh_3d[l,:,:]
    hsubs = hsubs_3d[l,:,:]   
    # filter over imask (3x3 kernel moving average) 
    fmask = scipy.ndimage.convolve(imask.astype(float), kernel*1./kernel.size)
    # cell as active (newmask=1) when >ncell out of (kernel.size)^2 cells in the kernel are active 
    newmask = np.zeros_like(imask,dtype=int) # 1 if cell is active, otherwise 0
    newmask[fmask>ncell/kernel.size]=1
    # -- interpolate permh 
    ma = np.ma.masked_invalid(permh)
    x1, y1, v1 = xx[~ma.mask], yy[~ma.mask], ma[~ma.mask]
    ipermh = scipy.interpolate.griddata((x1, y1), v1, (xx, yy))
    newpermh = nans.copy()
    newpermh[newmask==1] = ipermh[newmask==1]
    # -- interpolate hsubs (for newly activated cells without geometry) 
    ma = np.ma.masked_invalid(hsubs)
    x1, y1, v1 = xx[~ma.mask], yy[~ma.mask], ma[~ma.mask]
    ihsubs = scipy.interpolate.griddata((x1, y1), v1, (xx, yy),method='nearest')
    newhsubs = nans.copy()
    newhsubs[newmask==1] = ihsubs[newmask==1]
    # permh to nan where geometry is nan
    newpermh[np.isnan(newhsubs)]=np.nan
    # set data
    newmask_list.append(newmask)
    newpermh_list.append(newpermh)
    newhsubs_list.append(newhsubs)


# 3d grids from stack
newmask_3d = np.stack(newmask_list)
newpermh_3d = np.stack(newpermh_list)
newhsubs_3d = np.stack(newhsubs_list)

# interpolate topo (for newly activated cells without geometry)
ma = np.ma.masked_invalid(topo)
x1, y1, v1 = xx[~ma.mask], yy[~ma.mask], ma[~ma.mask]
itopo = scipy.interpolate.griddata((x1, y1), v1, (xx, yy),method='nearest')
# cell consider as valid where any of newpermh_3d or former topo is invalid
newmask = np.logical_or( np.any(~np.isnan(newpermh_3d), axis=0), ~np.isnan(topo))
newtopo = nans.copy()
newtopo[newmask==1] = itopo[newmask==1]

# ----------------------------------------------
# ----- fix new geometry (deactivate invalid cells)
# ----------------------------------------------

# initialize counter
nfix = 1

# iteratively fix geometry 
while nfix>0 :
    # arrays of interface elevations
    itfs_3d = np.concatenate((newtopo[None,:,:],newhsubs_3d),axis=0)
    # arrays of thickness (accounting for gaps)
    mask = np.isnan(itfs_3d)
    idx = np.where(~mask,(np.ones(mask.shape[::-1])*np.arange(mask.shape[0])).T,0)
    idx = np.maximum.accumulate(idx, axis=0,dtype=int)
    itfs_3d[mask] = itfs_3d[ idx[mask], np.nonzero(mask)[1], np.nonzero(mask)[2]]
    thk_3d = -1*np.diff(itfs_3d,axis=0)
    # deactivate cells below thickness threshold
    nfix = np.logical_and(thk_3d<thk_thd_2,~np.isnan(newpermh_3d)).sum()
    print(f'Fixed {nfix} cells with invalid thicknesses')
    newpermh_3d[thk_3d<thk_thd_2] = np.nan
    newhsubs_3d[thk_3d<thk_thd_2] = np.nan
    print('Minimum thickness : ' + str(np.nanmin(thk_3d[~np.isnan(newpermh_3d)])))


# ----------------------------------------------
# ----- plotting
# ----------------------------------------------


print('Plots for permh...')
for l in range(mm.nlay):
    fig, ax = plt.subplots(figsize=(5,4))
    ax = ax.imshow(permh_3d_org[l,:,:])
    fig.colorbar(ax)
    fig.savefig(os.path.join('figs',f'l{l+1:02d}_permh_org.png'))

for l in range(mm.nlay):
    fig, ax = plt.subplots(figsize=(5,4))
    ax = ax.imshow(permh_3d[l,:,:])
    fig.colorbar(ax)
    fig.savefig(os.path.join('figs',f'l{l+1:02d}_permh_crop.png'))

for l in range(mm.nlay):
    fig, ax = plt.subplots(figsize=(5,4))
    ax = ax.imshow(newpermh_3d[l,:,:])
    fig.colorbar(ax)
    fig.savefig(os.path.join('figs',f'l{l+1:02d}_permh_new.png'))

print('Plots for hsubs...')
for l in range(mm.nlay):
    fig, ax = plt.subplots(figsize=(5,4))
    ax = ax.imshow(hsubs_3d_org[l,:,:])
    fig.colorbar(ax)
    fig.savefig(os.path.join('figs',f'l{l+1:02d}_hsubs_org.png'))

for l in range(mm.nlay):
    fig, ax = plt.subplots(figsize=(5,4))
    ax = ax.imshow(newhsubs_3d[l,:,:])
    fig.colorbar(ax)
    fig.savefig(os.path.join('figs',f'l{l+1:02d}_hsubs_new.png'))

print('Plots for thicknesses...')

for l in range(mm.nlay):
    fig, ax = plt.subplots(figsize=(5,4))
    ax = ax.imshow(thk_3d_org[l,:,:])
    fig.colorbar(ax)
    fig.savefig(os.path.join('figs',f'l{l+1:02d}_thk_org.png'))

for l in range(mm.nlay):
    fig, ax = plt.subplots(figsize=(5,4))
    ax = ax.imshow(thk_3d[l,:,:])
    fig.colorbar(ax)
    fig.savefig(os.path.join('figs',f'l{l+1:02d}_thk_new.png'))

plt.close('all')

# ----------------------------------------------
# ----- save to files 
# ----------------------------------------------

# nan to no data values
newpermh_3d[np.isnan(newpermh_3d) ]=0
newpermh_3d[np.isnan(newhsubs_3d) ]=0
newhsubs_3d[np.isnan(newhsubs_3d)]=9999
newtopo[np.isnan(newtopo)]=9999

# append 5 layers with 0 to fit MartheField format (quick and dirty)
newtopo_3d = np.concatenate((newtopo[None,:,:],np.zeros_like(newpermh_3d[1:,:,:])),axis=0)

# convert to MartheField
newpermh_mf = MartheField('permh', newpermh_3d, mm, use_imask=False)
newhsubs_mf = MartheField('hsubs', newhsubs_3d, mm, use_imask=False)
newtopo_mf = MartheField('topog', newtopo_3d, mm, use_imask=False)

# write data to file 
newpermh_mf.write_data()
newhsubs_mf.write_data()
newtopo_mf.write_data()

