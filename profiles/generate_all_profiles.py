"""
This mini script grabs all the cosmo files and generates splitted profiles according to the 
approach of Hsieh for the three-profile splitting, and according to Mullins for generation
of a single profile.  Note that the appropriate effective areas are used for each approach.

Also, the single-profiles are generated according to the averaging approach of Hsieh.
"""
import glob
import os
import timeit
import functools
import multiprocessing

from to_sigma import read_Dmol3, write_sigma

def write_one(cosmo_path, *, num_profiles, averaging):
    # print(cosmo_path)
    try:
        dmol = read_Dmol3(inpath=cosmo_path, num_profiles=num_profiles, averaging=averaging)
        inchikey = os.path.split(cosmo_path)[1].split('.')[0]
        outpath = 'UD/sigma'+str(num_profiles)+'/'+inchikey+'.sigma'
        write_sigma(dmol, outpath, force = True)
    except BaseException as BE:
        print('ERROR:', BE, cosmo_path)
        # raise

if __name__ == '__main__':
    files = glob.glob('UD/cosmo/*.cosmo')

    # # SERIAL
    # for num_profiles in [3]:
    #     tic = timeit.default_timer()
    #     for ic, cosmo in enumerate(files):
    #         if ic % 50 == 0:
    #             print(ic,'out of', len(files),':', cosmo)
    #         write_one(cosmo, num_profiles = num_profiles, averaging='Hsieh')
    #     print(timeit.default_timer()-tic,'s')

    # PARALLEL
    f = functools.partial(write_one, num_profiles = 3, averaging='Hsieh')
    pool = multiprocessing.Pool(2)
    tic = timeit.default_timer()
    pool.map(f, files)
    print(timeit.default_timer()-tic,'s to write all profiles')