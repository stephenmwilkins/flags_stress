



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



SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')

model = 'bpass-2.2.1'

out_dir = f'Wilkins22/{model}'
param_filename = f'Wilkins22.{model}.spectra.param'

print(SPS.grid['log10age'])
print(SPS.grid['log10Z'])

i = 0



# parameters = {}
#
# parameters['sfzh'] = {}
# parameters['sfzh'] = {'model': 'constant', 'log10_duration': 7., 'log10Z': -2.4, 'log10M*': 1.}

param_file = []

def write_template(name, SED):
    global i, param_file
    i += 1
    llam = SED.total.lnu/SED.total.lam**2
    # np.savetxt(f'{out_dir}/{name}.dat', np.array([np.round(SED.total.lam, 3), np.round(llam, 3)]).T)
    np.savetxt(f'{out_dir}/{i}.dat', np.array([np.round(SED.total.lam, 3), np.round(llam, 3)]).T)
    # param_file.append(f'{i} {out_dir}/{name}.dat 1.0 {10**(log10age-9):.2f} 1.0')
    param_file.append(f'{i} templates/{out_dir}/{i}.dat 1.0 {10**(log10age-9):.2f} 1.0')


# --- young with no nebular emission
log10age = 7.

sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': log10age, 'log10Z': -2.4, 'log10M*': 1.})
SED = SPS.get_Lnu(sfzh, {'fesc': 1.0}, dust = None)
write_template(f'{model}-constant-log10age:7-log10Z:-2.4-fesc:1-tauV:0.0', SED)


# --- young with nebular emission
sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': log10age, 'log10Z': -2.4, 'log10M*': 1.})
SED = SPS.get_Lnu(sfzh, {'fesc': 0.0}, dust = None)
write_template(f'{model}-constant-log10age:7-log10Z:-2.4-fesc:0-tauV:0.0', SED)


log10age = 8.

# --- older with nebular emission
sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': log10age, 'log10Z': -2.4, 'log10M*': 1.})
SED = SPS.get_Lnu(sfzh, {'fesc': 0.0}, dust = None)
write_template(f'{model}-constant-log10age:8-log10Z:-2.4-fesc:0-tauV:0.0', SED)


# --- dusty constant templates
for tauV in [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]:
    sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': log10age, 'log10Z': -2.4, 'log10M*': 1.})
    SED = SPS.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': np.log10(tauV)}, dust = ('simple', {'slope':-1}))
    write_template(f'{model}-constant-log10age:8-log10Z:-2.4-fesc:0-tauV:{tauV:.2f}', SED)


# --- instantaneous bursts
for log10age in [8.0, 8.5, 9.0, 9.5, 10.0]:
    sfzh, sfr = interrogator.sed.sfzh.instantaneous(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10age': log10age, 'log10Z': -2.4, 'log10M*': 1.})
    SED = SPS.get_Lnu(sfzh, {'fesc': 0.0}, dust = None)
    write_template(f'{model}-instant-log10age:{log10age:.1f}-log10Z:-2.4-fesc:0-tauV:0.0', SED)

print(param_file)


with open(param_filename, 'w') as f:
    for line in param_file:
        f.write(f'{line}\n')
