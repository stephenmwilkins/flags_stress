


import numpy as np
import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh
import flare.filters
import flare.plt as fplt
from interrogator.sed.core import rebin
import h5py
from flags.pz import eazy
import time


uniform = lambda low, high, N: np.random.uniform(low = low, high = high, size = N)

np.random.seed(4484)


template_combos = {}
template_combos['tweak_fsps_QSF_12_v3'] = 'a'
template_combos['Larson22'] = 'a'
template_combos['Wilkins22.bpass-2.2.1'] = 'a'


def generate_galaxies():

    parameter_values = {}
    for parameter in parameters:
        parameter_values[parameter] = np.random.uniform(*parameter_range[parameter], size = N)

    parameter_values['log10duration'][1] = np.log10(cosmo.age(z).to('yr').value)

    print(parameter_values['log10duration'])

    hf = h5py.File(f'data/{scenario}/{z:.2f}.hf', 'w')

    hf[f'parameters/z'] = np.ones(N)*z

    for parameter in parameters:
        hf[f'parameters/{parameter}'] = parameter_values[parameter]

    hf[f'diagnostics/beta'] = np.zeros(N)

    for filter in filters:
        hf[f'fnu/{filter}'] = np.zeros(N)

    for i, (log10duration, log10Z, fesc, log10tau_V) in enumerate(zip(parameter_values['log10duration'], parameter_values['log10Z'],parameter_values['fesc'],parameter_values['log10tau_V'])):
        sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': log10duration, 'log10Z': log10Z, 'log10M*': 1.})
        SED = SPS.get_Lnu(sfzh, {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))
        SED.total.get_fnu(cosmo, z) # calculate observer frame wavelength
        SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
        for filter in filters:
            hf[f'fnu/{filter}'][i] = SED.total.Fnu[filter]
        beta = SED.total.return_beta()
        hf[f'diagnostics/beta'][i] = beta
        print(i, beta, log10duration, log10Z, fesc, log10tau_V)

    return hf

def make_observations(hf = False):

    if not hf:
        hf = h5py.File(f'data/{scenario}/{z:.2f}.hf', 'a')
    N = len(hf['parameters/z'][:])


    hf['obs/filters'] = filters
    hf['obs/pivwv'] = [F[filter].pivwv() for filter in filters]

    # --- rescale flux and add noise

    scaling = (SNR*noise[reference_filter])/hf['fnu/'+reference_filter][()]

    for filter in filters:
        hf['obs/'+filter+'/flux_err'] = np.ones(N)*noise[filter]
        hf['obs/'+filter+'/flux'] = hf['fnu/'+filter][()]*scaling + noise[filter]*np.random.normal(size=N)

    return hf

# --- photometric redshifts


def run_eazy(hf = False):

    print(hf)

    if not hf:
        hf = h5py.File(f'data/{scenario}/{z:.2f}.hf', 'a')
    N = len(hf['parameters/z'][:])


    id = rf'{scenario}-{z:.2f}'

    for template in templates:

        # --- initialise EAZY fitter
        pz = eazy.Eazy(id, F)

        pz.params['TEMPLATES_FILE'] = f'templates/{template}.spectra.param'
        pz.params['TEMPLATE_COMBOS'] = template_combos[template]
        pz.params['OBS_SED_FILE'] = 'y'                  # Write out observed SED/object, .obs_sed
        pz.params['TEMP_SED_FILE'] = 'y'                   # Write out best template fit/object, .temp_sed
        pz.params['POFZ_FILE'] = 'y'

        # --- create input catalogue from HDF5 object
        pz.create_input_catalogue_from_HDF5(hf)
        pz.run()

        eazy.append_EAZY_output_to_HDF5(f'EAZY/outputs/{id}', hf, read_pz = True, read_template_norm = True, group_name = 'pz/eazy/'+template)

    return hf






if __name__ == "__main__":


    hf = False

    parameter_range = {}
    parameter_range['log10duration'] = [6., 9.]
    parameter_range['log10Z'] = [-5., -1.4]
    parameter_range['fesc'] = [0., 1]
    parameter_range['log10tau_V'] = [-2., 1.]

    parameters = list(parameter_range.keys())

    cosmo = flare.default_cosmo()

    scenario = 'constant'
    SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')

    z = np.float(sys.argv[1])
    N = int(sys.argv[2])

    filters = [f'JWST.NIRCAM.{f}' for f in ['F090W','F115W','F150W','F200W','F277W','F356W','F410M','F444W']]
    F = flare.filters.add_filters(filters, new_lam = SPS.lam*(1.+z))

    hf = generate_galaxies()

    noise = {'JWST.NIRCAM.F090W': 6.67, 'JWST.NIRCAM.F115W': 6.08, 'JWST.NIRCAM.F150W': 4.83, 'JWST.NIRCAM.F200W': 3.98, 'JWST.NIRCAM.F277W': 4.70, 'JWST.NIRCAM.F356W': 3.91, 'JWST.NIRCAM.F410M': 7.56, 'JWST.NIRCAM.F444W': 5.95}
    reference_filter = 'JWST.NIRCAM.F200W'
    SNR = 20.

    hf = make_observations(hf = hf)

    templates = ['tweak_fsps_QSF_12_v3','Larson22', 'Wilkins22.bpass-2.2.1']
    # templates = ['Larson22']
    # templates = ['Wilkins22.bpass-2.2.1']
    # templates = ['tweak_fsps_QSF_12_v3']

    hf = run_eazy(hf = hf)

    hf.flush()
