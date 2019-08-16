import cCOSMO
import os
import itertools
import glob
import pandas
import timeit

tic = timeit.default_timer()

# Construct the (empty) database
here = os.path.abspath(os.path.dirname(__file__))
db = cCOSMO.DelawareProfileDatabase(here+"/UD/complist.txt", 
                                    here+"/UD/sigma3/")

# Find all the sigma profiles, add them to database
inchis = []
for filepath in glob.glob(here+"/UD/sigma3/*.sigma"):
    inchi = os.path.split(filepath)[1].split('.')[0]
    inchis.append(inchi)
    # if len(inchis) > 20: break # You can uncomment this line to control the number of profiles that get loaded
    db.add_profile(inchi)

print(timeit.default_timer()-tic,'s to load all the profiles')

# Conditions at which the activity coefficients are to be evaluated
T = 298.15 # [K]
z = [0.3, 0.7]
# At temperature T and mole fractions of z, carry out the COSMO calculations
data = []
# Iterate over each pair of fluids, 
for ip, pair in enumerate(itertools.combinations(inchis, 2)):
    # Do only every 20th point to cut down on the number of calculations
    if ip % 20 != 0:
        continue
    # Print every 100th time an activity coefficient is calculated
    # to make sure the calculations are not stuck
    if ip % (100*20) == 0:
        print(ip)
    COSMO = cCOSMO.COSMO3(pair, db);
    COSMO.get_mutable_COSMO_constants().fast_Gamma = True # Turn on fast mode
    COSMO.get_mutable_COSMO_constants().Gamma_rel_tol = 1e-14
    lngamma = COSMO.get_lngamma(T, z)
    lngamma_comb = [COSMO.get_lngamma_comb(T, z, 0), COSMO.get_lngamma_comb(T, z, 1)]
    lngamma_disp = COSMO.get_lngamma_disp(z)
    data.append(dict(
        lngamma0 = lngamma[0], lngamma1 = lngamma[1], 
        lngamma0_comb = lngamma_comb[0], lngamma1_comb = lngamma_comb[1],
        lngamma0_disp = lngamma_disp[0], lngamma1_disp = lngamma_disp[1],
        T_K=T, InChIKey0 = pair[0], InChIKey1 = pair[1],
        z0_molar = z[0], z1_molar = z[1]
    ))
pandas.DataFrame(data).to_csv('validation_data.csv')