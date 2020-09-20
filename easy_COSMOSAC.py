import cCOSMO
import numpy as np
import os
import timeit

here = os.path.abspath(os.path.dirname(__file__))

def calc_LNAC(*, method, dbname, names, composition, T):

    if dbname == 'UD':
        db = cCOSMO.DelawareProfileDatabase(
            here+"/profiles/UD/complist.txt",
            here+"/profiles/UD/sigma3/")
    elif dbname == 'VT':
        if method != 'COSMOSAC-2002':
            raise ValueError("Not possible to use VT2005 database with this method; need 3 sigma profiles")
        db = cCOSMO.VirginiaTechProfileDatabase(
            here+"/profiles/VT2005/Sigma_Profile_Database_Index_v2.txt",
            here+"/profiles/VT2005/Sigma_Profiles_v2/")

    else:
        raise ValueError('Invalid database name: ' + dbname)

    # Add the fluids we want into the database
    for name in names:
        iden = db.normalize_identifier(name)
        db.add_profile(iden)

    assert(len(names) == len(composition))
    if method == 'COSMOSAC-2002':
        COSMO = cCOSMO.COSMO1(names, db)
    elif method in ['COSMOSAC-2010','COSMOSAC-dsp']:
        COSMO = cCOSMO.COSMO3(names, db)

        # Specialize for 2010, no dispersive contribution, but residual and combinatorial
        if method  == 'COSMOSAC-2010':
            return COSMO.get_lngamma_comb(T, composition) + COSMO.get_lngamma_resid(T, composition)
    else:
        raise ValueError('Invalid method: ' + method)

    # If we haven't already returned, then do the calculation of ln(gamma)
    return COSMO.get_lngamma(T, composition)

if __name__ == '__main__':
    for db in ['UD','VT']:
        for method in ['COSMOSAC-2002','COSMOSAC-2010','COSMOSAC-dsp']:
            try:
                names=['ETHANOL', 'METHANOL']
                composition=[0.3, 0.7]
                T=300
                print(db, method,names,T,'K',composition, calc_LNAC(T=T, dbname=db, names = names, composition = composition, method=method))
                #print(db, method, calc_LNAC(T=300, dbname=db, names=['ETHANOL','METHANOL'], composition = [0.3, 0.7], method=method))
            except BaseException as BE:
                print(BE)
