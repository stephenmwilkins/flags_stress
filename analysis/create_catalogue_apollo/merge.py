

import sys
import numpy as np
import h5py


keys2 = []

def append(k):
    keys2.append(k)




if __name__ == "__main__":

    scenario = 'constant'

    hfo = h5py.File(f'out/{scenario}/all.h5', 'w')

    # --- calculate length of output array
    N = 0
    for i in range(1000):
        z = 7. + (0.01 * i)

        with h5py.File(f'out/{scenario}/{z}.hf', 'r') as hf:
            N += len(hf['z'])

    print(f'total number of objects: {N}')



    with h5py.File(f'out/{scenario}/{z}.hf', 'r') as hf:

        keys = list(hf.keys())
        hf.visit(append)

    print(keys)
    print(keys2)



    #
    #
    # o = {k: np.empty(N) for k in keys}
    #
    # print(o)
    #
    # rn = 0
    #
    # for z in np.arange(5., 16., 1.):
    #
    #     input = h5py.File(f'{t}_{z}.h5', 'r')
    #
    #     n = len(input['z'])
    #     print(rn, n)
    #     for k in keys:
    #         print(len(input[k][()]))
    #
    #         o[k][rn:rn+n] = input[k][()]
    #
    #     rn += n
