

import sys
import numpy as np
import h5py
from unyt import yr, Myr
import time

from synthesizer.filters import SVOFilterCollection
from synthesizer.grid import SpectralGrid
from synthesizer.binned import SFH, ZH, generate_sfzh, SEDGenerator
from synthesizer.plt import single, single_histxy, mlabel
from flags.pz import eazy








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

    hf = h5py.File(f'out/{scenario}/{z:.2f}.hf', 'w')

    hf[f'parameters/z'] = np.ones(N)*z

    for parameter in parameters:
        hf[f'parameters/{parameter}'] = parameter_values[parameter]

    hf[f'diagnostics/beta'] = np.zeros(N)

    for filter in filters:
        hf[f'fnu/{filter}'] = np.zeros(N)

    for i, (log10duration, log10Z, fesc, fesc_LyA, log10tauV) in enumerate(zip(parameter_values['log10duration'], parameter_values['log10Z'],parameter_values['fesc'], parameter_values['fesc_LyA'], parameter_values['log10tauV'])):


        sfh = SFH_({'duration': 10**log10duration * Myr })
        Zh = ZH_({'log10Z': log10Z})

        sfzh = generate_sfzh(grid.log10ages, grid.metallicities, sfh, Zh, stellar_mass = 1E8)

        galaxy = SEDGenerator(grid, sfzh)
        galaxy.pacman(fesc = fesc, fesc_LyA = fesc_LyA, tauV = 10**log10tauV)

        sed = galaxy.spectra['total'] # choose total SED
        sed.get_fnu(cosmo, z) # generate observed frame spectra

        # --- measure broadband fluxes
        fluxes = sed.get_broadband_fluxes(fc)

        for filter in filters:
            hf[f'fnu/{filter}'][i] = fluxes[filter].value
            print(filter, fluxes[filter].value)

        beta = sed.return_beta()
        hf[f'diagnostics/beta'][i] = beta
        # print(i, beta, log10duration, log10Z, fesc, fesc_LyA, log10tauV)
        # print(fluxes)

    return hf



def make_observations(hf = False):

    if not hf:
        hf = h5py.File(f'data/{scenario}/{z:.2f}.hf', 'a')
    N = len(hf['parameters/z'][:])


    hf['obs/filters'] = filters
    hf['obs/pivwv'] = [fc.filter[filter].pivwv() for filter in filters]

    # --- rescale flux and add noise

    scaling = (SNR*noise[reference_filter])/hf['fnu/'+reference_filter][()]

    for filter in filters:
        hf['obs/'+filter+'/flux_err'] = np.ones(N)*noise[filter]
        hf['obs/'+filter+'/flux'] = hf['fnu/'+filter][()]*scaling + noise[filter]*np.random.normal(size=N)

        snr = hf['obs/'+filter+'/flux'][()]/hf['obs/'+filter+'/flux_err'][()]
        print(filter, np.median(snr))


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
        pz = eazy.Eazy(id, fco)

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

    scenario = 'constant'

    parameter_range = {}

    parameter_range['log10Z'] = [-5., -1.4]
    parameter_range['fesc'] = [0., 1]
    parameter_range['fesc_LyA'] = [0., 1]
    parameter_range['log10tauV'] = [-2., 1.0]

    if scenario == 'constant':
        parameter_range['log10duration'] = [0., 3.]
        SFH_ = SFH.Constant # constant star formation
        ZH_ = ZH.deltaConstant # constant metallicity


    parameters = list(parameter_range.keys())

    from astropy.cosmology import Planck18 as cosmo

    grid_name = 'fsps-v3.2_Chabrier03_cloudy-v17.03_log10Uref-2'
    grid = SpectralGrid(grid_name)

    z = 7.+0.01*(float(sys.argv[1])-1)
    z = 7.
    N = 10

    # --- calculate broadband luminosities
    filters = [f'JWST/NIRCam.{f}' for f in ['F090W', 'F115W','F150W','F200W','F277W','F356W','F410M','F444W']] # define a list of filter codes
    fc = SVOFilterCollection(filters, new_lam = grid.lam * (1.+z)) # used for synthesizer
    fco = SVOFilterCollection(filters) # used for EAZY

    t1 = time.time()

    hf = generate_galaxies()

    noise = {'JWST/NIRCam.F090W': 6.67, 'JWST/NIRCam.F115W': 6.08, 'JWST/NIRCam.F150W': 4.83, 'JWST/NIRCam.F200W': 3.98, 'JWST/NIRCam.F277W': 4.70, 'JWST/NIRCam.F356W': 3.91, 'JWST/NIRCam.F410M': 7.56, 'JWST/NIRCam.F444W': 5.95}
    reference_filter = 'JWST/NIRCam.F200W'
    SNR = 20.

    hf = make_observations(hf = hf)

    templates = ['tweak_fsps_QSF_12_v3','Larson22', 'Wilkins22.bpass-2.2.1']
    # templates = ['Larson22']
    # templates = ['Wilkins22.bpass-2.2.1']
    templates = ['tweak_fsps_QSF_12_v3']

    hf = run_eazy(hf = hf)

    t2 = time.time()
    print((t2-t1)/N)

    hf.flush()
