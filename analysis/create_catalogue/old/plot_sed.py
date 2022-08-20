

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib import cm
import cmasher as cmr

from flags.pz import eazy

import flare.plt as fplt



template_set = 'tweak_fsps_QSF_12_v3.spectra.param'
templates = eazy.get_templates(template_set)



left  = 0.15
bottom = 0.15
width = 0.8
height = 0.8

lam_range = [3, 4]

fig = plt.figure(figsize = (5., 4))
ax = fig.add_axes((left, bottom, width, height))

scenario, id = 'constant-10.00', 12


obs_sed = ascii.read(f'EAZY/outputs/{scenario}_{id}.obs_sed')
ax.scatter(obs_sed['lambda'].data, obs_sed['flux_cat'].data)


# --- 1: plot the best fitting template
z, template_norm, lam, fnu = eazy.read_temp_sed(f'EAZY/outputs/{scenario}', id)
s = lam<50000.
ax.plot(lam[s], fnu[s])

# print(template_norm)
# print(eazy.read_template_norm(f'EAZY/outputs/{scenario}', id))


# --- 2: plot individual templates with the appropriate weighting
for i, t in enumerate(templates.values()):
    s = t.lam<50000.
    ax.plot(t.lam[s]*(1+z), template_norm[i]*t.fnu[s], lw = 1, alpha = 0.3)

# --- 3: plot the weighted combination of templates (should be the same as 1)
temp_lam, temp_fnu = eazy.get_template_grid(templates)
fnu = np.sum(temp_fnu*template_norm, axis=0)
ax.plot(temp_lam[s]*(1+z), fnu[s], lw = 3, alpha = 0.1)


fig.savefig('figs/sed.pdf')
