

import sys
import numpy as np
import h5py


datasets = []

def append(k, obj):
    if type(obj)==h5py._hl.dataset.Dataset:
        datasets.append(k)




if __name__ == "__main__":

    scenario = sys.argv[1]

    hfo = h5py.File(f'out/{scenario}/all.h5', 'w')

    # --- calculate length of output array
    N = 0
    for i in range(2000):
        z = 0.01 * (i+1)
        with h5py.File(f'out/{scenario}/{z:.2f}.hf', 'r') as hf:
            N += len(hf['parameters/z'])

    print(f'total number of objects: {N}')


    # --- make a list of keys
    z = 7.00
    with h5py.File(f'out/{scenario}/{z:.2f}.hf', 'r') as hf:
        hf.visititems(append)

    # datasets.remove('obs/filters')
    # datasets.remove('obs/pivwv')


    for k in datasets:
        print(k)
        hfo.create_dataset(k, data = np.empty(N))


    rn = 0

    for i in range(2000):

        z = 0.01 * (i+1)

        with h5py.File(f'out/{scenario}/{z:.2f}.hf', 'r') as hf:

            n = len(hf['parameters/z'])
            print(i, n)

            for k in datasets:
                hfo[k][rn:rn+n] = hf[k][()]

            rn += n

    hfo.flush()
