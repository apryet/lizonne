
# ------ parameterization of soil properties ----

print('Parameterize soil properties ...')
mm.load_prop('soil')
ms = mm.prop['soil']
skmi = pest_utils.get_kmi(ms, keys=['soilprop', 'zone'], istep=0) # Constant in time
mopt.add_param(parname='soil', mobj=ms, kmi=skmi, parlbnd=0, parubnd=1e5)

# --- Parameterization of  perm_r ----
# -> zpc based on aff_r id

print('Parameterize `perm_r` field ...')
# set up izone for perm_r based on reach id  
mm.load_prop('aff_r',use_imask=False) # reach id 
aff_r = mm.prop['aff_r']
iperm_r = MartheField('iperm_r', 1, mm,use_imask=False)
iperm_r.data['value'] = -1 * aff_r.data['value'] # negative values for ZPCs

# initialize perm_r field from scratch 
perm_r = MartheField('perm_r', 0, mm,use_imask=False)
perm_r.data['value'][aff_r.data['value']!=0]=1e-5
mm.prop['perm_r'] = perm_r

# set up parameter
mopt.add_param(parname='perm_r', mobj=mm.prop['perm_r'],
               izone=iperm_r,
               trans='log10', btrans='lambda x: 10**x',
               parlbnd=1e-8, parubnd=1, defaultvalue=1e-5,
               use_imask=False)

# --- parameterization of distributed hydraulic properties (permh, emmca, emmli)

# set pilot point spacings depending on layer type 
spacing_ep_layer = 8*1000 # in meters
spacing_aq_layer = 4*1000 # in meters
buffer = 4*1000 # in meters
ep_mask = mm.layers_infos.name.str.startswith('Eponte')
spacings = ep_mask.map({ True : spacing_ep_layer, False : spacing_aq_layer})

# --- Parameterization of  permh  ----
print('Parameterize `permh` field ...')

# izone value of 1 => pilot points for all layers 
ipermh = MartheField('ipermh', 1, mm, use_imask=True)

# pilot points set up
pp_permh = pp_utils.PilotPoints(ipermh)
for l,s in enumerate(spacings):
    pp_permh.add_spacing_pp(layer=l, zone=1, xspacing=s, yspacing=s, buffer=buffer)


# add parameters 
mopt.add_param(parname='permh', mobj=mm.prop['permh'],
               izone=ipermh, pp_data=pp_permh.to_pp_data(),
               trans='log10', btrans='lambda x: 10**x',
               parlbnd=1e-12, parubnd=1) 

# compute kriging factors
mopt.write_kriging_factors(pp_permh.extract_vgm_range(), parname='permh', cleanup=True, save_cov=True)

# --- Parameterization of  emmca  ---
print('Parameterize `emmca` field ...')

# izone definition: 
iemmca = MartheField('iemmca', 1, mm, use_imask=True)
iemmca.set_data(-1,layer=0)

# pilot points set up
pp_emmca = pp_utils.PilotPoints(iemmca)
for l,s in enumerate(spacings):
    # if pp for this layer
    if iemmca.data[iemmca.data['layer']==l]['value'].max() > 0 :
        # seed pps
        pp_emmca.add_spacing_pp(layer=l, zone=1, xspacing=s, yspacing=s, buffer=buffer)

# add parameters 
mm.load_prop('emmca')
mopt.add_param(parname='emmca', mobj=mm.prop['emmca'],
               izone=iemmca, pp_data=pp_emmca.to_pp_data(),
               trans='log10', btrans='lambda x: 10**x',
               parlbnd=1e-15, parubnd=1) 

# compute kriging factors
mopt.write_kriging_factors(pp_emmca.extract_vgm_range(), parname='emmca', cleanup=False, save_cov=True)

# --- Parameterization of  emmli ---
# --- parameterization of distributed hydraulic properties (emmli)
# zpc for all layers except layer 2:"Coniacien-Santonien" (layer_id = 1)
print('Parameterize `emmli` field ...')

# initialize emmli field from scratch
emmli_value = 0.1 
mm.prop['emmli'] =  MartheField('emmli', emmli_value, mm)

# izone definition: pp for aquifers only, zpc for aquitards
iemmli = MartheField('iemmli', 1, mm, use_imask=True)
_ = [ iemmli.set_data(-1,layer=i) for i,ep in enumerate(ep_mask) if ep]

# pilot points set up
pp_emmli = pp_utils.PilotPoints(iemmli)
#pp_emmli.add_spacing_pp(layer=1, zone=1, xspacing=spacings[1], yspacing=spacings[1])

for l,s in enumerate(spacings):
    # if pp for this layer
    if iemmli.data[iemmli.data['layer']==l]['value'].max() > 0 :
        # seed pps
        pp_emmli.add_spacing_pp(layer=l, zone=1, xspacing=s, yspacing=s, buffer=buffer)

mopt.add_param(parname='emmli', mobj=mm.prop['emmli'],
               izone=iemmli, pp_data=pp_emmli.to_pp_data(),
               trans='log10', btrans='lambda x: 10**x',
               parlbnd=1e-4, parubnd=1, defaultvalue=emmli_value)

# compute kriging factors
mopt.write_kriging_factors(pp_emmli.extract_vgm_range(), parname='emmli', cleanup=False, save_cov=True)

