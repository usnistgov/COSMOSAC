"""
Backwards compatibility tests for existing COSMO formats.

This test validates that modifications to to_sigma.py do not affect
the existing support for Gaussian09, DMol3, and GAMESS formats.
"""
import os
import sys
import numpy as np

# Script location is ORCA_TEST; parent is profiles
script_dir = os.path.dirname(os.path.abspath(__file__))
profiles_dir = os.path.dirname(script_dir)
sys.path.insert(0, profiles_dir)

from to_sigma import read_Dmol3, write_sigma

# For pytest compatibility - make it optional
try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False


def read_sigma_file(path):
    """
    Read a .sigma file and return sigma values and profile data as numpy arrays.
    """
    sigmas = []
    psigmaA_values = []
    
    with open(path, 'r') as fp:
        for line in fp:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    sigma = float(parts[0])
                    psigmaA = float(parts[1])
                    sigmas.append(sigma)
                    psigmaA_values.append(psigmaA)
                except ValueError:
                    pass
    
    return np.array(sigmas), np.array(psigmaA_values)


def test_gaussian09_test_suite():
    """Test that Gaussian09 COSMO files still parse and regenerate correctly."""
    test_cases = [
        ('GAUSSIAN09_TEST', 'ethanol.cosmo', 'ethanol.sigma'),
    ]
    
    for test_dir, cosmo_file, sigma_file in test_cases:
        test_path = os.path.join(profiles_dir, test_dir)
        
        if not os.path.exists(test_path):
            print(f"⊘ Skipping {test_dir}/{cosmo_file} (directory not found)")
            continue
        
        input_file = os.path.join(test_path, cosmo_file)
        reference_sigma = os.path.join(test_path, sigma_file)
        
        if not os.path.exists(input_file) or not os.path.exists(reference_sigma):
            print(f"⊘ Skipping {test_dir}/{cosmo_file} (files not found)")
            continue
        
        temp_output = os.path.join(test_path, f'{cosmo_file.split(".")[0]}_test.sigma')
        
        try:
            # Parse and regenerate
            dmol = read_Dmol3(inpath=input_file, num_profiles=3, averaging='Hsieh')
            write_sigma(dmol, temp_output, force=True)
            
            # Read both files
            gen_sigmas, gen_psigmaA = read_sigma_file(temp_output)
            ref_sigmas, ref_psigmaA = read_sigma_file(reference_sigma)
            
            # Check sigma grid
            np.testing.assert_array_almost_equal(
                gen_sigmas, ref_sigmas, decimal=3,
                err_msg=f"Sigma grid mismatch in {test_dir}/{cosmo_file}"
            )
            
            # Check profile values with reasonable tolerance
            for i, (gen_val, ref_val) in enumerate(zip(gen_psigmaA, ref_psigmaA)):
                if abs(ref_val) > 1e-14:
                    rel_error = abs(gen_val - ref_val) / abs(ref_val)
                    assert rel_error < 1e-6, (
                        f"Gaussian mismatch in {test_dir}/{cosmo_file} at sigma={gen_sigmas[i]:.3f}: "
                        f"rel_error={rel_error:.6e}"
                    )
                else:
                    abs_error = abs(gen_val - ref_val)
                    assert abs_error < 1e-14
            
            print(f"✓ Gaussian09: {test_dir}/{cosmo_file}")
            
        finally:
            if os.path.exists(temp_output):
                os.remove(temp_output)


def test_dmol3_test_suite():
    """Test that DMol3 COSMO files still parse and regenerate correctly."""
    test_path = os.path.join(profiles_dir, 'DMol3_TEST')
    
    # DMol3 test uses IUPAC code as filename
    cosmo_file = 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N.cosmo'
    sigma_file = 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N.sigma'
    
    input_file = os.path.join(test_path, cosmo_file)
    reference_sigma = os.path.join(test_path, sigma_file)
    
    if not os.path.exists(input_file) or not os.path.exists(reference_sigma):
        print("⊘ Skipping DMol3 test (files not found)")
        return
    
    temp_output = os.path.join(test_path, f'{cosmo_file.split(".")[0]}_test.sigma')
    
    try:
        # Parse and regenerate
        dmol = read_Dmol3(inpath=input_file, num_profiles=3, averaging='Hsieh')
        write_sigma(dmol, temp_output, force=True)
        
        # Read both files
        gen_sigmas, gen_psigmaA = read_sigma_file(temp_output)
        ref_sigmas, ref_psigmaA = read_sigma_file(reference_sigma)
        
        # Check sigma grid
        np.testing.assert_array_almost_equal(
            gen_sigmas, ref_sigmas, decimal=3,
            err_msg="Sigma grid mismatch in DMol3_TEST"
        )
        
        # Check profile values
        for i, (gen_val, ref_val) in enumerate(zip(gen_psigmaA, ref_psigmaA)):
            if abs(ref_val) > 1e-14:
                rel_error = abs(gen_val - ref_val) / abs(ref_val)
                assert rel_error < 1e-6
            else:
                abs_error = abs(gen_val - ref_val)
                assert abs_error < 1e-14
        
        print("✓ DMol3: DMol3_TEST/ethanol")
        
    finally:
        if os.path.exists(temp_output):
            os.remove(temp_output)


def test_gamess_test_suite():
    """Test that GAMESS .gout files still parse and regenerate correctly."""
    test_path = os.path.join(profiles_dir, 'GAMESS_TEST')
    
    cosmo_file = 'ETHANOL.gout'
    sigma_file = 'ETHANOL.sigma'
    
    input_file = os.path.join(test_path, cosmo_file)
    reference_sigma = os.path.join(test_path, sigma_file)
    
    if not os.path.exists(input_file) or not os.path.exists(reference_sigma):
        print("⊘ Skipping GAMESS test (files not found)")
        return
    
    temp_output = os.path.join(test_path, f'{cosmo_file.split(".")[0]}_test.sigma')
    
    try:
        # Parse and regenerate
        dmol = read_Dmol3(inpath=input_file, num_profiles=3, averaging='Hsieh')
        write_sigma(dmol, temp_output, force=True)
        
        # Read both files
        gen_sigmas, gen_psigmaA = read_sigma_file(temp_output)
        ref_sigmas, ref_psigmaA = read_sigma_file(reference_sigma)
        
        # Check sigma grid
        np.testing.assert_array_almost_equal(
            gen_sigmas, ref_sigmas, decimal=3,
            err_msg="Sigma grid mismatch in GAMESS_TEST"
        )
        
        # Check profile values
        for i, (gen_val, ref_val) in enumerate(zip(gen_psigmaA, ref_psigmaA)):
            if abs(ref_val) > 1e-14:
                rel_error = abs(gen_val - ref_val) / abs(ref_val)
                assert rel_error < 1e-6
            else:
                abs_error = abs(gen_val - ref_val)
                assert abs_error < 1e-14
        
        print("✓ GAMESS: GAMESS_TEST/ETHANOL")
        
    finally:
        if os.path.exists(temp_output):
            os.remove(temp_output)


if __name__ == '__main__':
    print("Running backwards compatibility tests...\n")
    test_gaussian09_test_suite()
    test_dmol3_test_suite()
    test_gamess_test_suite()
    print("\n✓ All backwards compatibility tests passed!")
