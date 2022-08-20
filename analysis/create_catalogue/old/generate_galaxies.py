

import numpy as np
import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh
import flare.filters
import flare.plt as fplt
from interrogator.sed.core import rebin
import h5py
import time

uniform = lambda low, high, N: np.random.uniform(low = low, high = high, size = N)


# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')

print(SPS.grid['log10Z'])
print(SPS.grid['log10age'])

# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)




parameter_range = {}
parameter_range['log10_duration'] = [6., 9.]
parameter_range['log10Z'] = [-5., -1.4]
parameter_range['fesc'] = [0., 1]
parameter_range['log10tau_V'] = [-2., 1.]

# parameter_range['log10_duration'] = [7.5, 8.5]
# parameter_range['log10Z'] = [-1.5,-2.5]
# parameter_range['fesc'] = [0., 1]
# parameter_range['log10tau_V'] = [-3., -2.]

parameters = list(parameter_range.keys())

cosmo = flare.default_cosmo()
filters = [f'JWST.NIRCAM.{f}' for f in ['F090W','F115W','F150W','F200W','F277W','F356W','F410M','F444W']]
print(filters)




z = 10.0
N = 1000


# for z in np.arange(5., 10., 0.1)

F = flare.filters.add_filters(filters, new_lam = SPS.lam*(1.+z))

parameter_values = {}
for parameter in parameters:
    parameter_values[parameter] = np.random.uniform(*parameter_range[parameter], size = N)

parameter_values['log10_duration'][1] = np.log10(cosmo.age(z).to('yr').value)


hf = h5py.File(f'data/constant/{z:.2f}.hf', 'w')


hf[f'parameters/z'] = np.ones(N)*z

for parameter in parameters:
    hf[f'parameters/{parameter}'] = parameter_values[parameter]

hf[f'diagnostics/beta'] = np.zeros(N)

for filter in filters:
    hf[f'fnu/{filter}'] = np.zeros(N)


t1 = time.time()

for i, (log10_duration, log10Z, fesc, log10tau_V) in enumerate(zip(parameter_values['log10_duration'], parameter_values['log10Z'],parameter_values['fesc'],parameter_values['log10tau_V'])):

    sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': log10_duration, 'log10Z': log10Z, 'log10M*': 1.})
    SED = SPS.get_Lnu(sfzh, {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))
    SED.total.get_fnu(cosmo, z) # calculate observer frame wavelength
    SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
    for filter in filters:
        hf[f'fnu/{filter}'][i] = SED.total.Fnu[filter]

    beta = SED.total.return_beta()
    hf[f'diagnostics/beta'][i] = beta
    print(beta, log10_duration, log10Z, fesc, log10tau_V)

hf.flush()

t2 = time.time()
print(f'total time take per object: {(t2-t1)/N}')
