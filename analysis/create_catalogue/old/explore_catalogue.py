

import numpy as np
import h5py




z = 5.00

hf = h5py.File(f'data/constant/{z:.2f}.hf', 'r')

beta = hf['diagnostics/beta']

print(np.min(beta), np.median(beta), np.max(beta))


min_beta = np.min(beta)
