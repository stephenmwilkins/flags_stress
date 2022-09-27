

import sys
import numpy as np
import h5py


keys = []

def append(k):
    keys.append(k)




if __name__ == "__main__":

    scenario = 'constant'

    hfo = h5py.File(f'out/{scenario}/all.h5', 'w')

    # --- calculate length of output array
    N = 0
    for i in range(1000):
        z = 7. + (0.01 * i)
        with h5py.File(f'out/{scenario}/{z:.2f}.hf', 'r') as hf:
            N += len(hf['parameters/z'])

    print(f'total number of objects: {N}')


    # --- make a list of keys
    z = 7.00
    with h5py.File(f'out/{scenario}/{z:.2f}.hf', 'r') as hf:
        hf.visit(append)

    for k in keys:
        print(k)
        hfo.create_dataset(k, np.empty(N))


    rn = 0

    for i in range(1000):

        z = 7. + (0.01 * i)

        with h5py.File(f'out/{scenario}/{z:.2f}.hf', 'r') as hf:

            n += len(hf['parameters/z'])

            for k in keys:
                hfo[k][rn:rn+n] = hf[k][()]

            rn += n

    hfo.flush()
