"""
This module enhances the COSMO-SAC package with the `DirectImport` feature,
providing a more flexible and organized method for managing sigma profiles.
Unlike the standard import method requiring sigma profiles in a single folder,
`DirectImport` allows for storing profiles in separate directories and specifying
the path and name of the sigma file for each component. This approach facilitates
enhanced organization and flexible testing of different sigma profiles for the
same component.

Author: Ivan Antolovic
E-Mail: ivan.antolovic@tu-berlin.de

Note: The `DirectImport` method originated during the work on the paper:
https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.4c00342
Example profiles are taken from: https://github.com/ivanantolo/cosmopharm
"""

import cCOSMO
from pathlib import Path

# Get the directory where the script is located
script_dir = Path(__file__).resolve().parent

# Constants
PROFILES_API = script_dir / "profiles/pharmaceuticals"
PROFILES_POLY = script_dir / "profiles/polymers"
PATH_TO_SIGMAS = script_dir / "profiles/sigma3"
PATH_TO_COMPLIST = script_dir / "profiles/complist.txt"
PATH_TO_PROFILES = [PROFILES_API, PROFILES_POLY]

def import_delaware(names, path_to_sigmas, path_to_complist):
    """
    Imports sigma profiles using the DelawareProfileDatabase.

    Parameters:
    - names: A list of component names.
    - path_to_sigmas: The file path to the sigma profiles.
    - path_to_complist: The file path to the component list.

    Returns:
    A tuple containing the COSMO3 object and the database object.

    Considerations:
    - All sigma profiles must be located in the same directory.
    - All sigma profiles must be listed in complist.txt.
    - Using .xlsx files for adding new profiles and converting them to .txt is recommended.
    - Much of the information in complist.txt is not necessary for importing sigma profiles.
    - Names used in complist.txt do not correspond to unique identifiers.
    - One name corresponds to one identifier, which poses a challenge when dealing with modifications.
    - The same name may be associated with different .sigma files.
    """
    db = cCOSMO.DelawareProfileDatabase(str(path_to_complist), str(path_to_sigmas))
    for name in names:
        db.add_profile(name)
    return cCOSMO.COSMO3(names, db), db

def import_direct(names, paths):
    """
    Directly imports sigma profiles using DirectImport.

    Parameters:
    - names: A list of component names.
    - paths: A list of paths to the sigma profiles.

    Returns:
    A tuple containing the COSMO3 object and the database object.

    Advantages:
    - Profiles can be stored in separate directories.
    - No need for a .txt file listing all profiles.
    - Provides flexibility for testing different sigma profiles for the same component.
    - Enhances the organization of profiles.
    """
    db = cCOSMO.DirectImport()
    for name, path in zip(names, paths):
        db.add_profile(name, str(path))
    return cCOSMO.COSMO3(names, db), db

def display_profile(db, name):
    """
    Displays information about a single profile from the database object.

    Parameters:
    - db: The database object containing the profiles.
    - name: The name for the profile to display.
    """
    try:
        profile = db.get_profile(name)
        print(f"Name: {profile.name}")
        print(f"Surface Area (A^2): {profile.A_COSMO_A2}")
        print(f"Volume (A^3): {profile.V_COSMO_A3}")
        # print(f"Sigma Profile (non-hydrogen-bonding segments, first 5 elements): {profile.profiles.nhb.sigma[:5]}")
        # print(f"Probability (non-hydrogen-bonding segments, first 5 elements): {profile.profiles.nhb.psigmaA[:5]}")
    except Exception as e:
        print(f"Error retrieving profile for {name}: {e}")

if __name__ == "__main__":
    names = ['SIM', 'PLGA50']

    # Import profiles using both methods
    cosmo_delaware, db_delaware = import_delaware(names, PATH_TO_SIGMAS, PATH_TO_COMPLIST)
    cosmo_direct, db_direct = import_direct(names=names, paths=PATH_TO_PROFILES)

    # Display information about the imported profiles
    print("\nProfiles imported using DelawareProfileDatabase:")
    for name in names:
        display_profile(db_delaware, name)

    print("\nProfiles imported using DirectImport:")
    for name in names:
        display_profile(db_direct, name)
