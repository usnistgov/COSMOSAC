import cCOSMO
import numpy as np
import os
import timeit

here = os.path.abspath(os.path.dirname(__file__))
dbVT = cCOSMO.VirginiaTechProfileDatabase(
    here+"/profiles/VT2005/Sigma_Profile_Database_Index_v2.txt", 
    here+"/profiles/VT2005/Sigma_Profiles_v2/")
db = cCOSMO.DelawareProfileDatabase(
    here+"/profiles/UD/complist.txt", 
    here+"/profiles/UD/sigma3/")

# Load the profiles we are going to consider
names = [ "METHANE", "ETHANE", "WATER" ]
for iden in names:
    db.add_profile(db.normalize_identifier(iden))
    dbVT.add_profile(dbVT.normalize_identifier(iden))

for pair in [['METHANE', 'ETHANE'], ['METHANE', 'WATER']]:
    o = []
    for COSMO in [cCOSMO.COSMO1(pair, dbVT), cCOSMO.COSMO3(pair, db)]:
        T = 400.15;
        z = np.array([0.5, 0.5])
        cConsts = COSMO.get_mutable_COSMO_constants()
        Nreps = 1000
        vals = []
        times = []
        for fast_Gamma in [False, True]:
            cConsts.fast_Gamma = fast_Gamma
            tic = timeit.default_timer()
            for i in range(100):
                lngamma = COSMO.get_lngamma(T, z)
            toc = timeit.default_timer()
            elap = (toc-tic)*1e6/Nreps
            vals.append(lngamma)
            times.append(elap)
        assert(np.sum(np.array(vals[0]) - np.array(vals[1])) < 1e-9)
        o.append(times)
    print(' + '.join(pair) + ' & ' + '{0:0.1f} & {1:0.1f}'.format(*o[0]) + ' & {0:0.1f} & {1:0.1f}'.format(*o[1]))

import matplotlib.pyplot as plt
for i in names:
    prof = db.get_profile(db.normalize_identifier(i))
    print(prof.name)
    nhb = prof.profiles.nhb
    plt.plot(nhb.sigma, nhb.psigmaA/np.sum(nhb.psigmaA))

plt.xlabel(r'$\sigma$ / e/$\AA^2$')
plt.ylabel(r'$p(\sigma)$ ')
plt.show()