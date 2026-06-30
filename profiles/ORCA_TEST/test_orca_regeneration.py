"""
Test for ORCA CPCM sigma profile regeneration.

This test validates that ethanol.sigma can be regenerated from ORCA CPCM data.
It serves as minimal reproducible validation required by the PR review process.
"""
import os
import sys
import numpy as np

# Script location is ORCA_TEST; parent is profiles
orca_test_dir = os.path.dirname(os.path.abspath(__file__))
profiles_dir = os.path.dirname(orca_test_dir)
sys.path.insert(0, profiles_dir)

from to_sigma import read_Dmol3, write_sigma, Dmol3COSMOParser

# For pytest compatibility - make it optional
try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False


def read_sigma_file(path):
    """
    Read a .sigma file and return sigma values and profile data as numpy arrays.
    
    Format:
    # meta: {...}
    # Rows are given as: sigma [e/A^2] followed by a space, then psigmaA [A^2]
    # In the case of three sigma profiles, the order is NHB, OH, then OT
    sigma1 psigmaA1
    sigma2 psigmaA2
    ...
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
                    # Skip lines that can't be parsed
                    pass
    
    return np.array(sigmas), np.array(psigmaA_values)


def test_orca_ethanol_regeneration():
    """
    Test that ethanol.sigma regenerates correctly from ORCA CPCM data.
    
    This is the minimal reproducible validation showing that:
    1. ORCA ethanol.cpcm can be parsed successfully
    2. The companion ethanol.cpcm_corr file is found and loaded
    3. The regenerated .sigma file matches the reference closely
    """
    input_file = os.path.join(orca_test_dir, 'ethanol.cpcm')
    reference_sigma = os.path.join(orca_test_dir, 'ethanol.sigma')
    
    # Generate a temporary output file
    temp_output = os.path.join(orca_test_dir, 'ethanol_regenerated.sigma')
    
    try:
        # Parse the ORCA CPCM input
        # This will raise an error if the .cpcm_corr file is missing or invalid
        parser = Dmol3COSMOParser(
            inpath=input_file,
            num_profiles=3,
            averaging='Hsieh'
        )
        
        # Get the outputs
        outputs = parser.get_outputs()
        
        # Write the regenerated .sigma file
        write_sigma(outputs, temp_output, force=True)
        
        # Read both the generated and reference .sigma files
        gen_sigmas, gen_psigmaA = read_sigma_file(temp_output)
        ref_sigmas, ref_psigmaA = read_sigma_file(reference_sigma)
        
        # Check that the sigma grids match (should be identical)
        np.testing.assert_array_almost_equal(
            gen_sigmas, ref_sigmas,
            decimal=3,
            err_msg="Sigma grid values don't match between generated and reference"
        )
        
        # Check that the profile values are close
        # Use a relative tolerance that accounts for numerical precision
        # Allow some variation due to floating-point differences
        relative_tolerance = 1e-6
        
        # For very small values, use absolute tolerance
        absolute_tolerance = 1e-14
        
        for i, (gen_val, ref_val) in enumerate(zip(gen_psigmaA, ref_psigmaA)):
            if abs(ref_val) > absolute_tolerance:
                rel_error = abs(gen_val - ref_val) / abs(ref_val)
                assert rel_error < relative_tolerance, (
                    f"Profile value mismatch at index {i} (sigma={gen_sigmas[i]:.3f}): "
                    f"generated={gen_val:.6e}, reference={ref_val:.6e}, "
                    f"relative error={rel_error:.6e}"
                )
            else:
                # For values close to zero, use absolute tolerance
                abs_error = abs(gen_val - ref_val)
                assert abs_error < absolute_tolerance, (
                    f"Profile value mismatch at index {i} (sigma={gen_sigmas[i]:.3f}): "
                    f"generated={gen_val:.6e}, reference={ref_val:.6e}, "
                    f"absolute error={abs_error:.6e}"
                )
        
        print("✓ ORCA ethanol.sigma regenerated successfully and matches reference")
        
    finally:
        # Clean up temporary file
        if os.path.exists(temp_output):
            os.remove(temp_output)


def test_orca_cpcm_detection():
    """Test that ORCA CPCM format is correctly detected."""
    input_file = os.path.join(orca_test_dir, 'ethanol.cpcm')
    
    with open(input_file, 'r') as fp:
        contents = fp.read()
    
    # Check that ORCA markers are present
    assert '# CARTESIAN COORDINATES (A.U.) + RADII (A.U.) + ATOMIC NUMBER' in contents
    assert '# SURFACE POINTS (A.U.)' in contents
    print("✓ ORCA CPCM format markers detected correctly")


def test_orca_cpcm_corr_file_required():
    """Test that missing .cpcm_corr file raises an error."""
    input_file = os.path.join(orca_test_dir, 'ethanol.cpcm')
    
    # Temporarily rename the .cpcm_corr file
    corr_file = input_file + '_corr'
    corr_file_backup = corr_file + '.backup'
    
    try:
        if os.path.exists(corr_file):
            os.rename(corr_file, corr_file_backup)
        
        # Attempting to parse without .cpcm_corr should raise FileNotFoundError
        error_raised = False
        error_message = None
        try:
            parser = Dmol3COSMOParser(
                inpath=input_file,
                num_profiles=3,
                averaging='Hsieh'
            )
        except FileNotFoundError as e:
            error_raised = True
            error_message = str(e)
        
        assert error_raised, "FileNotFoundError was not raised for missing .cpcm_corr"
        assert "Missing ORCA corrected-charge file" in error_message, \
            f"Error message doesn't match. Got: {error_message}"
        
        print("✓ Missing .cpcm_corr file correctly raises FileNotFoundError")
        
    finally:
        # Restore the .cpcm_corr file
        if os.path.exists(corr_file_backup):
            os.rename(corr_file_backup, corr_file)


if __name__ == '__main__':
    print("Running ORCA regeneration tests...\n")
    test_orca_cpcm_detection()
    test_orca_cpcm_corr_file_required()
    test_orca_ethanol_regeneration()
    print("\n✓ All ORCA regeneration tests passed!")
