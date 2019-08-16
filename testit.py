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

names = [ "METHANE", "ETHANE" ]
for iden in names:
    db.add_profile(db.normalize_identifier(iden))
    dbVT.add_profile(dbVT.normalize_identifier(iden))

for COSMO in [cCOSMO.COSMO1(names, dbVT), cCOSMO.COSMO3(names, db)]:
    T = 400.15;
    z = np.array([0, 1])
    # print(COSMO.get_lngamma(T, z))

    # Print the constants in use (debugging use, not generally necessary)
    # print(COSMO.get_mutable_combinatorial_constants())
    # print(COSMO.get_mutable_COSMO_constants())

    psigma_mix = 0
    xA = 0
    for _z, name in zip(z,names):
        prof = db.get_profile(db.normalize_identifier(name))
        psigma_mix += _z*prof.profiles.nhb.psigmaA
        xA += prof.A_COSMO_A2
    psigma_mix /= xA

    # Change a constant
    comb_consts = COSMO.get_mutable_combinatorial_constants()
    comb_consts.z_coordination = 10.0

    cConsts = COSMO.get_mutable_COSMO_constants()
    cConsts.fast_Gamma = False

    Nreps = 1000
    for fast_Gamma in [False,True]:
        cConsts.fast_Gamma = fast_Gamma
        print('fast_Gamma:', fast_Gamma)
        tic = timeit.default_timer()
        for i in range(100):
            lngamma = COSMO.get_lngamma(T, z)
        toc = timeit.default_timer()
        print(lngamma)
        print((toc-tic)*1e6/Nreps,'us/call')

    # tic = timeit.default_timer()
    # Gamma_mix = COSMO.get_Gamma(T, psigma_mix)
    # toc = timeit.default_timer()
    # print('t(get_Gamma):',toc-tic)
    # print('Gamma_mix:', Gamma_mix)
    print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')

import matplotlib.pyplot as plt
for i in names:
    prof = db.get_profile(db.normalize_identifier(i))
    print(prof.name)
    nhb = prof.profiles.nhb
    plt.plot(nhb.sigma, nhb.psigmaA/np.sum(nhb.psigmaA))

plt.xlabel(r'$\sigma$ / e/$\AA^2$')
plt.ylabel(r'$p(\sigma)$ ')
plt.show()