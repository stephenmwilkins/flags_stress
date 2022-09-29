


import h5py




if __name__ == "__main__":

    scenario = 'constant'


    # --- calculate length of output array
    N = 1000
    for i in range(1000):
        z = 7. + (0.01 * i)

        try:
            with h5py.File(f'out/{scenario}/{z:.2f}.hf', 'r') as hf:
                n = len(hf['parameters/z'])
                if n != N:
                    print(f'failed N ({n})')
                    print(f'qsub -t {i+1} -jc test run.sh')
        except:
            print(f'qsub -t {i+1} -jc test run.sh')
