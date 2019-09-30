"""
A small script that shows how to extract and plot the 
NHB, OH, OT portions of the sigma profiles

A part of usnistgov/COSMOSAC
"""
import os
import timeit
import json

import numpy as np
import matplotlib.pyplot as plt
import cCOSMO

here = os.path.abspath(os.path.dirname(__file__))
dbVT = cCOSMO.VirginiaTechProfileDatabase(
    here+"/profiles/VT2005/Sigma_Profile_Database_Index_v2.txt", 
    here+"/profiles/VT2005/Sigma_Profiles_v2/")
dbUD = cCOSMO.DelawareProfileDatabase(
    here+"/profiles/UD/complist.txt", 
    here+"/profiles/UD/sigma3/")

# Load the fluids into both databases (they start empty)
names = [ "METHANOL", "ETHANOL" ]
for iden in names:
    dbUD.add_profile(dbUD.normalize_identifier(iden))
    dbVT.add_profile(dbVT.normalize_identifier(iden))

# Plot the profiles
import matplotlib.pyplot as plt

# Iterate over the databases
for db,dbname in zip([dbVT, dbUD],['VT','UD']):
    # Iterate over the names
    for name, dashes in zip(names,[[2,2], []]):
        # Get the sigma profiles
        prof = db.get_profile(db.normalize_identifier(name))
        # Iterate over the profiles to be plotted
        for key in ['nhb', 'oh', 'ot']:
            profile = getattr(prof.profiles, key)
            PA = np.sum(profile.psigmaA)
            if PA > 0:
               plt.plot(profile.sigma, profile.psigmaA/PA, dashes=dashes, label=dbname+':'+name+':'+key)
            print(dbname, name, key, PA)

# Labeling and saving
plt.legend(loc='best')
plt.xlabel(r'$\sigma$ / e/$\AA^2$')
plt.ylabel(r'$p(\sigma)$ ')
plt.savefig('methanol_ethanol_profiles.pdf')
plt.show()