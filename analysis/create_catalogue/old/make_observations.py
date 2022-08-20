

import numpy as np
import h5py
from flags.pz import eazy
from flare.filters import add_filters


np.random.seed(4484)

z = 10.00

scenario = 'constant'
id = rf'{scenario}-{z:.2f}'
hf = h5py.File(f'data/{scenario}/{z:.2f}.hf', 'r')

hf.visit(print)

cat = h5py.File(f'data/{scenario}/{z:.2f}_obs.hf', 'w')

SNR = 20.

filters = [f'JWST.NIRCAM.{f}' for f in ['F090W','F115W','F150W','F200W','F277W','F356W','F410M','F444W']]

F = add_filters(filters) # make FLARE filter object

noise = {'JWST.NIRCAM.F090W': 6.67, 'JWST.NIRCAM.F115W': 6.08, 'JWST.NIRCAM.F150W': 4.83, 'JWST.NIRCAM.F200W': 3.98, 'JWST.NIRCAM.F277W': 4.70, 'JWST.NIRCAM.F356W': 3.91, 'JWST.NIRCAM.F410M': 7.56, 'JWST.NIRCAM.F444W': 5.95}

reference_filter = 'JWST.NIRCAM.F200W'


N = len(hf['diagnostics/beta'][()])

cat['obs/filters'] = filters
cat['obs/pivwv'] = [F[filter].pivwv() for filter in filters]

# --- redcale flux and add noise

scaling = (SNR*noise[reference_filter])/hf['fnu/'+reference_filter][()]

for filter in filters:
    cat['obs/'+filter+'/flux_err'] = np.ones(N)*noise[filter]
    cat['obs/'+filter+'/flux'] = hf['fnu/'+filter][()]*scaling + noise[filter]*np.random.normal(size=N)

for filter in filters:
    snr = cat['obs/'+filter+'/flux'][()]/cat['obs/'+filter+'/flux_err'][()]
    print(rf'{filter} {np.median(snr):.2f}')





# --- photometric redshifts


run_eazy = True

if run_eazy:

    for template, template_combos in [( 'tweak_fsps_QSF_12_v3', 'a')]:


        # --- initialise EAZY fitter
        pz = eazy.Eazy(id, F)

        pz.params['TEMPLATES_FILE'] = f'templates/{template}.spectra.param'
        pz.params['TEMPLATE_COMBOS'] = template_combos
        pz.params['OBS_SED_FILE'] = 'y'                  # Write out observed SED/object, .obs_sed
        pz.params['TEMP_SED_FILE'] = 'y'                   # Write out best template fit/object, .temp_sed
        pz.params['POFZ_FILE'] = 'y'

        # --- create input catalogue from HDF5 object
        pz.create_input_catalogue_from_HDF5(cat)
        pz.run()

        eazy.append_EAZY_output_to_HDF5(f'EAZY/outputs/{id}', cat, read_pz = True, read_template_norm = True, group_name = 'pz/eazy/'+template)


hf.flush()
cat.flush()
