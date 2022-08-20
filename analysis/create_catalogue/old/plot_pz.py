

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib import cm
import cmasher as cmr

from flags.pz import eazy

import flare.plt as fplt



left  = 0.15
bottom = 0.15
width = 0.8
height = 0.8

lam_range = [3, 4]

fig = plt.figure(figsize = (5., 4))
ax = fig.add_axes((left, bottom, width, height))

scenario = 'EAZY/outputs/constant-10.00'
id = 12

z, chi2, pz = eazy.read_pz(scenario, id)

ax.plot(z, pz)


fig.savefig('figs/pz_test.pdf')
