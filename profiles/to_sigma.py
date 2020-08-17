# Standard library
from __future__ import division
import re
import os
from io import StringIO
from math import exp
from collections import namedtuple
import timeit
import json
import itertools

# Conda packages
import pandas
import scipy.spatial.distance
import numpy as np
import matplotlib.pyplot as plt

# From https://doi.org/10.1039/b801115j
# Covalent radii in angstrom, used to determine bonding
covalent_radius = {
    'H':0.31,
    'He':0.28,
    'Li':1.28,
    'Be':0.96,
    'B':0.84,
    'C':0.76, # sp3 hybridization, sp2: 0.73 sp: 0.69
    'N':0.71,
    'O':0.66,
    'F':0.57,
    'Ne':0.58,
    'Na':1.66,
    'Mg':1.41,
    'Al':1.21,
    'Si':1.11,
    'P':1.07,
    'S':1.05,
    'Cl':1.02,
    'Ar':1.06,
    'K':2.03,
    'Ca':1.76,
    'Sc':1.70,
    'Ti':1.60,
    'V':1.53,
    'Cr':1.39,
    'Mn':1.39, # l.s.; h.s.: 1.61
    'Fe':1.32, # l.s.; h.s.: 1.52
    'Co':1.26, # l.s.; h.s.: 1.50
    'Ni':1.24,
    'Cu':1.32,
    'Zn':1.22,
    'Ga':1.22,
    'Ge':1.20,
    'As':1.19,
    'Se':1.20,
    'Br':1.20,
    'Kr':1.16,
    'Rb':2.20,
    'Sr':1.95,
    'Y':1.90,
    'Zr':1.75,
    'Nb':1.64,
    'Mo':1.54,
    'Tc':1.47,
    'Ru':1.46,
    'Rh':1.42,
    'Pd':1.39,
    'Ag':1.45,
    'Cd':1.44,
    'In':1.42,
    'Sn':1.39,
    'Sb':1.39,
    'Te':1.38,
    'I':1.39,
    'Xe':1.40,
    'Cs':2.44,
    'Ba':2.15,
    'La':2.07,
    'Ce':2.04,
    'Pr':2.03,
    'Nd':2.01,
    'Pm':1.99,
    'Sm':1.98,
    'Eu':1.98,
    'Gd':1.96,
    'Tb':1.94,
    'Dy':1.92,
    'Ho':1.92,
    'Er':1.89,
    'Tm':1.90,
    'Yb':1.87,
    'Lu':1.87,
    'Hf':1.75,
    'Ta':1.70,
    'W':1.62,
    'Re':1.51,
    'Os':1.44,
    'Ir':1.41,
    'Pt':1.36,
    'Au':1.36,
    'Hg':1.32,
    'Tl':1.45,
    'Pb':1.46,
    'Bi':1.48,
    'Po':1.40,
    'At':1.50,
    'Rn':1.50,
    'Fr':2.60,
    'Ra':2.21,
    'Ac':2.15,
    'Th':2.06,
    'Pa':2.00,
    'U':1.96,
    'Np':1.90,
    'Pu':1.87,
    'Am':1.80,
    'Cm':1.69
}
# Construct the complete list of possible lengths of bonds. If a distance between atoms is less 
# than the sum of covalent radii of the two atoms forming the possible bond, it is considered
# to be covalently bonded.
bond_distances = {
    (k1, k2): covalent_radius[k1] + covalent_radius[k2] for k1, k2 in itertools.combinations(covalent_radius.keys(), 2)
}
# Also put in the backwards pair (O,H instead of H,O)
for keys, d in bond_distances.copy().items():
    bond_distances[tuple(reversed(keys))] = d
# Also put in the atoms with themselves (e.g., C,C)
for key in covalent_radius.keys():
    bond_distances[(key, key)] = 2*covalent_radius[key]

def get_seg_DataFrame(COSMO_contents):
    # Read in the information.  Look for (X, Y, Z), and search to the end of the line, then capture until 
    # you get to a pair of two end-of-line characters, or an eol character followed by the end of string
    if "DMol3/COSMO Results" in COSMO_contents:
        sdata = re.search(r"\(X, Y, Z\)[\sa-zA-Z0-9\[\]/]+\n(.+)(\n\n|$)", COSMO_contents, re.DOTALL).group(1).rstrip()
        table_assign = ['n','atom','x / a.u.','y / a.u.','z / a.u.','charge / e','area / A^2','charge/area / e/A^2','potential']
    elif "Gaussian COSMO output" in COSMO_contents:
        # Gaussian09: same unit (Bohr in G09) as Dmol3 but format is different
        sdata = re.search(r"\(X, Y, Z\)[\sa-zA-Z0-9\[\]/\#\n]+\#\n([\s\S]+)", COSMO_contents, re.DOTALL).group(1).rstrip()
        table_assign = ['n','atom','x / a.u.','y / a.u.','z / a.u.','charge / e','area / A^2','charge/area / e/A^2','potential']
    elif "GAMESS/COSab RESULTS" in COSMO_contents:
        # GAMESS: same unit (Bohr in GAMESS/COSab) as Dmol3 but format is different
        sdata = re.search(r"\(X, Y, Z\)[\sa-zA-Z0-9\(\)\./\*]+\n([\s0-9\-\n.]+)(=+)", COSMO_contents, re.DOTALL).group(1).rstrip()
        table_assign = ['n','atom','x / a.u.','y / a.u.','z / a.u.','charge / e','area / A^2','charge/area / e/A^2','potential']
    # Annotate the columns appropriately with units(!)
    return pandas.read_csv(StringIO(sdata), names=table_assign, sep=r'\s+',engine= 'python')

def get_atom_DataFrame(COSMO_contents):
    # Read in the information
    if "DMol3/COSMO Results" in COSMO_contents:
        sdata = re.search(r"!DATE[a-zA-Z0-9:\s]+\n(.+)end\s+\nend", COSMO_contents, re.DOTALL).group(1)
        table_assign = ['atomidentifier','x / A','y / A','z / A','?1','?2','?3','atom','?4']
    elif "Gaussian COSMO output" in COSMO_contents:
        # Gaussian09: same unit as Dmol3 (Angstrom in G09) but format is slightly different
        sdata = re.search(r"!DATE[a-zA-Z0-9:\s]*\n(.+)end\s*\nend", COSMO_contents, re.DOTALL).group(1)
        table_assign = ['atomidentifier','x / A','y / A','z / A','?1','?2','?3','atom','?4']
    elif "GAMESS/COSab RESULTS" in COSMO_contents:
        # GAMESS: same unit (Angstrom in GAMESS) as Dmol3 but format is different
        sdata = re.search(r"EQUILIBRIUM GEOMETRY[\sa-zA-Z0-9\(\)\./\*\n]+\-+\n(\s[\s\S]+)\n\n\n", COSMO_contents, re.DOTALL).group(1)
        table_assign = ['atom','charge','x / A','y / A','z / A']
    # Annotate the columns appropriately with units(!)
    return pandas.read_csv(StringIO(sdata), names=table_assign, sep=r'\s+',engine = 'python')

def get_area_volume(COSMO_contents):
    # Read in the information
    if "DMol3/COSMO Results" in COSMO_contents:
        area = float(re.search(r"Total surface area of cavity \(A\*\*2\)     =(.+)\n", COSMO_contents).group(1).strip())
        volume = float(re.search(r"Total volume of cavity \(A\*\*3\)           =(.+)\n", COSMO_contents).group(1).strip())
    elif "Gaussian COSMO output" in COSMO_contents:
        # Gaussian09: "area" = area of the solute cavity surface (in Bohr^2)
        # Gaussian09: "volume" = volume enclosed by the solute cavity surface (in Bohr^3)
        area = float(re.search(r"area  =(.+)\n", COSMO_contents).group(1).strip())*(0.52917721067)**2
        volume = float(re.search(r"volume=(.+)\n", COSMO_contents).group(1).strip())*(0.52917721067)**3
    elif "GAMESS/COSab RESULTS" in COSMO_contents:
        # GAMESS: units of area and volume are same as DMol3
        area = float(re.search(r"Total surface area of cavity \(A\*\*2\)\s+=(.+)\n", COSMO_contents).group(1).strip())
        volume = float(re.search(r"Total volume of cavity \(A\*\*3\)\s+=(.+)\n", COSMO_contents).group(1).strip())

    return area, volume

def weightbin_sigmas(sigmavals, sigmas_grid):
    """
    """
    # Regular grid, so every bin is of same width
    bin_width = sigmas_grid[1] - sigmas_grid[0]

    psigmaA = np.zeros_like(sigmas_grid)
    for sigma, area in sigmavals:
        # Check sigma
        if sigma < np.min(sigmas_grid):
            raise ValueError('Sigma [{0:g}] is less than minimum of grid [{1}]'.format(sigma, np.min(sigmas_grid)))
        if sigma > np.max(sigmas_grid):
            raise ValueError('Sigma [{0:g}] is greater than maximum of grid [{1}]'.format(sigma, np.max(sigmas_grid)))
        # The index to the left of the point in sigma
        left = int((sigma-sigmas_grid[0])/bin_width)
        # Weighted distance from the left edge of the cell to the right side
        w_left = (sigmas_grid[left+1]-sigma)/bin_width
        # Add the area into the left and right nodes of the bin, each part weighted by the value 
        # of sigma, if equal to the left edge of the bin, then the w_left=1, if the right side, 
        # w_left = 0
        psigmaA[left] += area*w_left
        psigmaA[left+1] += area*(1.0-w_left)
    return psigmaA

Dmol3COSMO = namedtuple('Dmol3COSMO',['sigmas','psigmaA_nhb','psigmaA_OH','psigmaA_OT','df_atom','meta'])
DispersiveValues = namedtuple('DispersiveValues',['dispersion_flag', 'dispersive_molecule', "Nbonds", 'has_COOH'])

class Dmol3COSMOParser(object):

    def __init__(self, inpath, num_profiles = 3, averaging = None):


        # Open the COSMO file and read in its contents
        COSMO_contents = open(inpath).read()

        # Parse the COSMO file and get metadata
        self.df = get_seg_DataFrame(COSMO_contents)
        self.df_atom = get_atom_DataFrame(COSMO_contents)
        self.area_A2, self.volume_A3 = get_area_volume(COSMO_contents)

        averaging_options = ['Hsieh','Mullins']
        if averaging not in averaging_options:
            raise ValueError('averaging[' + averaging + '] not in '+str(averaging_options))

        self.num_profiles = num_profiles
        self.averaging = averaging

        # Convert coordinates in a.u. (actually, bohr) to Angstroms
        for field in ['x','y','z']:
            self.df[field + ' / A'] = self.df[field + ' / a.u.']*0.52917721067 # https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
        # Calculate the effective circular radius for this segment patch from its area
        self.df['rn / A']  = (self.df['area / A^2']/np.pi)**0.5
        self.df['rn^2 / A^2'] = self.df['rn / A']**2

        # Recalculate the charge density here because the vaues in COSMO file are unneccssarily truncated
        self.sigma = np.array(self.df['charge / e']/self.df['area / A^2'])
        self.rn2 = np.array(self.df['rn^2 / A^2'])

        assert(int(self.df.iloc[-1].n) == len(self.rn2))

        # Calculate the distances between each pair of segments in a euclidean sense
        XA = np.c_[self.df['x / A'], self.df['y / A'], self.df['z / A']]
        self.dist_mat_squared = scipy.spatial.distance.cdist(XA, XA, 'euclidean')**2

        # Calculate the distances between each atom center in a euclidean sense
        XA = np.c_[self.df_atom['x / A'], self.df_atom['y / A'], self.df_atom['z / A']]
        self.dist_mat_atom = scipy.spatial.distance.cdist(XA, XA, 'euclidean')

        # Set a flag if the molecule is water; water is treated specially in some cases
        self.is_water = self.df_atom.atom.tolist().count('H') == 2 and self.df_atom.atom.tolist().count('O') == 1 and len(self.df_atom) == 3

        # Calculate the dispersive values
        self.disp = self.get_dispersive_values()

        # Tag each atom with its hydrogen bonding class
        self.df_atom['hb_class'] = self.get_HB_classes_per_atom()

        # Store the number of bonds in the DataFrame
        self.df_atom['Nbonds'] = self.disp.Nbonds

        # Average the charge densities on each segment
        self.sigma_averaged = self.average_sigmas(self.sigma)

        # Split up the profiles
        self.sigma_nhb, self.sigma_OH, self.sigma_OT = self.split_profiles(self.sigma_averaged, self.num_profiles)

    def get_bonds(self, i):
        """
        Get the list of tuples of (j, atom_name_j) of things that are bonded to this atom
        """
        bonds = []
        atom_name_i = self.df_atom.atom.iloc[i]

        for j, distance in enumerate(self.dist_mat_atom[i]):
            if i == j: continue # atom cannot bond to itself
            atom_name_j = self.df_atom.atom.iloc[j]
            if distance < 1.15*bond_distances[(atom_name_i, atom_name_j)]:
                bonds.append((j, atom_name_j))
        return bonds

    def get_HB_classes_per_atom(self):
        # Determine the hydrogen bonding class of each atom
        # Options are:
        #
        # 1) nhb: the atom is not a candidate to hydrogen bond
        # 2) OH: the atom is either the oxygen or the hydrogen in a OH hydrogen-bonding group
        # 3) OT: the atom is N,F, or an oxygen that is not part of an OH bonding group
        #
        # Loop over each atom-atom pair, and determine if it is a candidate for hydrogen bonding
        hydrogen_bonding_atoms = []
        for i in range(len(self.df_atom)):
            atom_name_i = self.df_atom.atom.iloc[i]
            if atom_name_i in ['N','F']:
                # Definite hydrogen bonding atom, cannot be in OH class,
                # so we are done, set it to OT
                hydrogen_bonding_atoms.append('OT')

            elif atom_name_i in ['O','H']:
                # If hydrogen is bonded to O, N, or F, then it is hydrogen-bonding,
                # otherwise not.  The determination of OT or OH hydrogen-bonding class
                # depends on the atom that the hydrogen is bonded to.
                #
                atom_type = 'NHB'
                for j, distance in enumerate(self.dist_mat_atom[i]):
                    if i == j:
                        continue

                    js, atom_name_js = zip(*self.get_bonds(i))

                    # It's OH because an OH bond is found
                    if (atom_name_i == 'H' and 'O' in atom_name_js) or (atom_name_i == 'O' and 'H' in atom_name_js):
                        atom_type = 'OH'
                        break

                    # It's OT because it's an oxygen bonded to something other than H
                    if atom_name_i == 'O':
                        atom_type = 'OT'
                        break

                    # It's OT because H is bonded to F or N
                    if atom_name_i == 'H' and ('F' in atom_name_js or 'N' in atom_name_js):
                        atom_type = 'OT'
                        break     

                hydrogen_bonding_atoms.append(atom_type)

            else:
                # Definite non-hydrogen-bonding
                hydrogen_bonding_atoms.append('NHB')
        return hydrogen_bonding_atoms

    def get_dispersive_values(self):
        """
        Calculate the dispersive parameters needed for the COSMO-SAC-dsp
        model.

        Returns:
            vals (DispersiveValues): Instance of DispersiveValues namedtuple
        """

        # Determine the dispersive energy parameter of each atom
        dispersive_parameter_lib = {
            'C(sp3)': 115.7023,
            'C(sp2)': 117.4650,
            'C(sp)': 66.0691,
            '-O-': 95.6184,
            '=O': -11.0549,
            'N(sp3)': 15.4901,
            'N(sp2)': 84.6268,
            'N(sp)': 109.6621,
            'F': 52.9318,
            'Cl': 104.2534,
            'H(OH)': 19.3477,
            'H(NH)': 141.1709,
            'H(water)': 58.3301,
            'H(other)': 0
        }
        
        dispersive_molecule = 0
        invalid_atom = False
        has_COOH = False
        Natom_nonzero = 0
        dispersion_flag = 'NHB'
        Nbonds = []
        for i in range(len(self.df_atom)):
            dispersive_parameter = None
            atom_name_i = self.df_atom.atom.iloc[i]
            bonds = self.get_bonds(i)
            if len(bonds) == 0 and len(self.df_atom):
                raise ValueError("Atom "+atom_name_i+" in a polyatomic molecule has no bonds, this is not good.")
            # print(i, atom_name_i, bonds)
            Nbonds.append(len(bonds))

            if atom_name_i == 'C':
                if len(bonds) == 4:
                    dispersive_parameter = dispersive_parameter_lib['C(sp3)']
                elif len(bonds) == 3:
                    dispersive_parameter = dispersive_parameter_lib['C(sp2)']
                elif len(bonds) == 2:
                    dispersive_parameter = dispersive_parameter_lib['C(sp)']
                # If C is double-bonded to an O and single-bonded to the other O, then a candidate for 
                # being a carboxyl group, let's dig deeper
                #
                # First, we see if the carbon is bonded to two oxygens
                js, atom_name_js = zip(*bonds)
                if len(bonds) == 3 and atom_name_js.count('O') == 2:
                    # If either O is bonded to C and an H, and nothing else, then we are in a carboxyl group
                    for j, atom_name_j in bonds:
                        if not atom_name_j == 'O': continue
                        ks, atom_name_ks = zip(*self.get_bonds(j))
                        if sorted(atom_name_ks) == ['C','H']:
                            has_COOH = True
                            break
            elif atom_name_i == 'N':
                if len(bonds) == 3:
                    dispersive_parameter = dispersive_parameter_lib['N(sp3)']
                elif len(bonds) == 2:
                    dispersive_parameter = dispersive_parameter_lib['N(sp2)']
                elif len(bonds) == 1:
                    dispersive_parameter = dispersive_parameter_lib['N(sp)']
                else:
                    raise ValueError('N with a funny number of bonds: '+len(bonds))
            elif atom_name_i == 'O':
                if len(bonds) == 2:
                    dispersive_parameter = dispersive_parameter_lib['-O-']
                elif len(bonds) == 1:
                    dispersive_parameter = dispersive_parameter_lib['=O']
                else:
                    raise ValueError('O with a funny number of bonds: '+len(bonds))
            elif atom_name_i == 'F':
                dispersive_parameter = dispersive_parameter_lib['F']
            elif atom_name_i == 'Cl':
                dispersive_parameter = dispersive_parameter_lib['Cl']
            elif atom_name_i == 'H':
                js, atom_name_js = zip(*bonds)
                if self.is_water:
                    dispersive_parameter = dispersive_parameter_lib['H(water)']
                elif 'O' in atom_name_js:
                    dispersive_parameter = dispersive_parameter_lib['H(OH)']
                elif 'N' in atom_name_js:
                    dispersive_parameter = dispersive_parameter_lib['H(NH)']
                else:
                    dispersive_parameter = None
            else:
                # The atom is unmatched, therefore we will return an invalid number
                # because we don't know what to do about this, but we still
                # need to calculate the number of atoms, so don't exit
                invalid_atom = True
                dispersive_parameter = np.nan
                
            if dispersive_parameter is not None:
                Natom_nonzero += 1
                dispersive_molecule += dispersive_parameter

        # If the number of atoms with non-zero dispersion parameters is greater than
        # zero, then we divide by the number of nonzero parameters, otherwise not
        # because that would be a division by zero
        if Natom_nonzero > 0:
            dispersive_molecule /= Natom_nonzero

        # If at least one atom is invalid, return an invalid number
        if invalid_atom:
            dispersive_molecule = np.nan

        # Determine the HB-DONOR-ACCEPTOR or HB-ACCEPTOR flags
        possble_Hbonders = ['O', 'N', 'F']
        if any([(atom in set(self.df_atom.atom)) for atom in possble_Hbonders]):
            Hbond = False
            for i in range(len(self.df_atom)):
                atom_name_i = self.df_atom.atom.iloc[i]
                if atom_name_i in possble_Hbonders:
                    bonds = self.get_bonds(i)
                    js, atom_name_js = zip(*bonds)
                    if 'H' in atom_name_js:
                        Hbond = True
            if not Hbond:
                dispersion_flag = 'HB-ACCEPTOR'
            else:
                dispersion_flag = 'HB-DONOR-ACCEPTOR'

        return DispersiveValues(dispersion_flag, dispersive_molecule, Nbonds, has_COOH)

    def average_sigmas(self, sigmavals):
        """
        Calculate the averaged charge densities on each segment, in e/\AA^2
        """

        # This code also works, kept for demonstration purposes, but the vectorized
        # numpy code is much faster, if a bit harder to understand
        # def calc_sigma_m(m):
        #     increments = rn2*r_av2/(rn2+r_av2)*np.exp(-f_decay*dist_mat_squared[m,:]/(rn2+r_av2))
        #     return np.sum(increments*sigma)/np.sum(increments)
        # def calc_sigma_m_mat(m):
        #     increments = THEMAT[m,:]
        #     return np.sum(increments*sigma)/np.sum(increments)

        if self.averaging == 'Mullins':
            self.r_av2 = 0.8176300195**2  # [A^2]  # Also equal to ((7.5/pi)**0.5*0.52917721092)**2
            self.f_decay = 1.0
        elif self.averaging == 'Hsieh':
            self.r_av2 = (7.25/np.pi) # [A^2]
            self.f_decay = 3.57
        else:
            raise ValueError("??")
            
        THEMAT = np.exp(-self.f_decay*self.dist_mat_squared/(self.rn2+self.r_av2))*self.rn2*self.r_av2/(self.rn2+self.r_av2)
        return np.sum(THEMAT*sigmavals, axis=1)/np.sum(THEMAT,axis=1)

    def split_profiles(self, sigmavals, num_profiles):
        """
        Split the samples into 
        """

        if num_profiles == 1:
            sigma_nhb = [(s,a) for s,a in zip(sigmavals,  df['area / A^2'])]
            sigma_OH = None
            sigma_OT = None
        elif num_profiles == 3:

            indices = np.array(self.df.atom.astype(int)-1)
            self.df['atom_name'] = self.df_atom.atom[indices].reset_index(drop=True)
            self.df['hb_class'] = self.df_atom.hb_class[indices].reset_index(drop=True)

            # N.B. : The charge density in sigmavals is the averaged charge density, 
            # not the charge density of the segment itself!!
            mask_OH = (
                ((self.df.atom_name == 'O') & (sigmavals > 0.0) & (self.df.hb_class == 'OH'))
                |
                ((self.df.atom_name == 'H') & (sigmavals < 0.0) & (self.df.hb_class == 'OH'))
            )
            mask_OT = (
                (self.df.atom_name.isin(['O', 'N', 'F']) & (sigmavals > 0.0) & (self.df.hb_class == 'OT'))
                |
                ((self.df.atom_name == 'H') & (sigmavals < 0.0) & (self.df.hb_class == 'OT'))
            )
            mask_nhb = ~(mask_OT | mask_OH)
            sigma_nhb = np.c_[sigmavals[mask_nhb], self.df[mask_nhb]['area / A^2']]
            sigma_OH = np.c_[sigmavals[mask_OH], self.df[mask_OH]['area / A^2']]
            sigma_OT = np.c_[sigmavals[mask_OT], self.df[mask_OT]['area / A^2']]
        else:
            raise ValueError('Number of profiles [{0}] is invalid'.format(num_profiles))
        return sigma_nhb, sigma_OH, sigma_OT

    def get_meta(self):
        """
        Return a dictionary with the metadata about this set of sigma profiles
        """
        return {
            'name': '?',
            'CAS': '?',
            'area [A^2]': self.area_A2,
            'volume [A^3]': self.volume_A3,
            'r_av [A]': self.r_av2**0.5,
            'f_decay': self.f_decay,
            'sigma_hb [e/A^2]': 0.0084,
            'averaging': self.averaging
        }

    def get_outputs(self):

        # The points where the surface-charge density 
        # will be evaluated, -0.025 to 0.025, in increments of 0.001
        bin_width = 0.001
        sigmas = np.arange(-0.025, 0.025+0.0001, bin_width) # [e/A^2]

        if self.num_profiles == 1:
            psigmaA = weightbin_sigmas(self.sigma_nhb, sigmas)
            # print('cumulative time after reweighting', timeit.default_timer()-tic,'s')
            assert(abs(sum(psigmaA)-meta['area [A^2]']) < 0.001)
            return Dmol3COSMO(sigmas, psigmaA, None, None, self.df_atom, self.get_meta())

        elif self.num_profiles == 3:
            psigmaA_nhb = weightbin_sigmas(self.sigma_nhb, sigmas)
            psigmaA_OH = weightbin_sigmas(self.sigma_OH, sigmas)
            psigmaA_OT = weightbin_sigmas(self.sigma_OT, sigmas)
            sigma_0 = 0.007 # [e/A^2]

            psigmaA_hb = psigmaA_OT + psigmaA_OH

            P_hb = 1 - np.exp(-sigmas**2/(2*sigma_0**2))
            psigmaA_OH *= P_hb
            psigmaA_OT *= P_hb

            psigmaA_nhb = psigmaA_nhb + psigmaA_hb*(1-P_hb)

            dispersion_flag = 'NHB'

            meta = self.get_meta()
            
            # Determine dispersion flag for the molecule
            if self.is_water:
                dispersion_flag = 'H2O'
            elif self.disp.has_COOH:
                dispersion_flag = 'COOH'
            else:
                dispersion_flag = self.disp.dispersion_flag

            meta['disp. flag'] = dispersion_flag
            meta['disp. e/kB [K]'] = self.disp.dispersive_molecule
            # print('cumulative time after reweighting', timeit.default_timer()-tic,'s')
            return Dmol3COSMO(sigmas, psigmaA_nhb, psigmaA_OH, psigmaA_OT, self.df_atom, meta)

def read_Dmol3(**kwargs):
    """
    A convenience function that will pass all arguments along to class and then return the outputs
    """
    return Dmol3COSMOParser(**kwargs).get_outputs()

def overlay_profile(sigmas, psigmaA, path):
    for profile in psigmaA:
        plt.plot(sigmas, profile, 'o-', mfc='none')
        print(np.c_[sigmas, profile])
    df = pandas.read_csv(path, names= ['charge/area / e/A^2', 'p(sigma)*A / A^2'], skiprows = 4, sep = ' ')
    plt.plot(df['charge/area / e/A^2'], df['p(sigma)*A / A^2'], '^-', mfc='none')
    plt.xlabel(r'$\sigma$ / e/A$^2$')
    plt.ylabel(r'$p(\sigma)A_i$ / A$^2$')
    plt.show()

Delaware_template = """# meta: {meta:s}
# Rows are given as: sigma [e/A^2] followed by a space, then psigmaA [A^2]
# In the case of three sigma profiles, the order is NHB, OH, then OT
"""

def write_sigma(dmol, ofpath, header = 'Delaware', force = True):

    if header == 'Delaware':
        out = Delaware_template.format(meta = json.dumps(dmol.meta))
    else:
        raise ValueError("Bad header option")
    summ = 0
    for profile in [dmol.psigmaA_nhb, dmol.psigmaA_OH, dmol.psigmaA_OT]:
        if profile is not None:
            for ir in range(profile.shape[0]):
                out += '{0:0.3f} {1:17.14e}\n'.format(dmol.sigmas[ir], profile[ir])
                summ += profile[ir]

    if os.path.exists(ofpath) and not force:
        raise ValueError("File [{0:s}] already exists and force has not been requested".format(ofpath))

    with open(ofpath,'w') as fp:
        fp.write(out)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate a sigma profile from a COSMO file')
    parser.add_argument('--inpath', type=str, required=True, nargs=1, help='The path to the cosmo file that is to be processed')
    parser.add_argument('--outpath', type=str, required=True, nargs=1, help='The path to the sigma profile file that is to be generated')
    parser.add_argument('--n', type=int, nargs=1, required=True, choices=[1,3], help='The number of profiles to generate, either 1 or 3')
    parser.add_argument('--averaging', type=str, nargs=1, required=True, choices=['Hsieh','Mullins'], help="The scheme used to do averaging of profiles, either 'Mullins' to use f_decay = 1 and r_av = 0.8176300195 A or 'Hsieh' to use f_decay = 3.57 and r_av = sqrt(7.25/pi)")

    args = parser.parse_args()
    # args = parser.parse_args('--n 3 --averaging Mullins --inpath UD/cosmo/ZPTVNYMJQHSSEA-UHFFFAOYSA-N.cosmo --outpath ZPTVNYMJQHSSEA-UHFFFAOYSA-N.sigma'.split(' '))  # For testing


    # For DMol3 .cosmo file test (ethanol)
    #args = parser.parse_args('--n 3 --averaging Hsieh --inpath DMol3_TEST/LFQSCWFLJHTTHZ-UHFFFAOYSA-N.cosmo --outpath DMol3_TEST/LFQSCWFLJHTTHZ-UHFFFAOYSA-N.sigma'.split(' '))  # For testing
    # For gaussian09 .cosmo file test
    #args = parser.parse_args('--n 3 --averaging Hsieh --inpath GAUSSIAN09_TEST/ethanol.cosmo --outpath GAUSSIAN09_TEST/ethanol.sigma'.split(' '))
    # For GAMESS .gout (including cosmo data) file test
    #args = parser.parse_args('--n 3 --averaging Hsieh --inpath GAMESS_TEST/ETHANOL.gout --outpath GAMESS_TEST/ETHANOL.sigma'.split(' '))

    try:
        dmol = read_Dmol3(inpath = args.inpath[0], num_profiles=args.n[0], averaging=args.averaging[0])
        write_sigma(dmol, args.outpath[0])
    except BaseException as BE:
        print(BE)
        print(args)
