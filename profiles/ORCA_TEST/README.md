# Profile Validation (ORCA + Backward Compatibility)

This folder contains the ORCA fixture and scripts to validate `to_sigma.py`.

## Run validation

From this folder (`profiles/ORCA_TEST`):

```bash
python test_orca_regeneration.py
python test_backwards_compatibility.py
```

Both commands should finish with all checks passing.

## What is validated

1. ORCA fixture regeneration (`ORCA_TEST/ethanol.cpcm` + `ethanol.cpcm_corr` -> `ethanol.sigma`)
2. Strict `.cpcm_corr` requirement (missing file must fail)
3. No regressions for existing formats in the parent `profiles` folder:
   - DMol3: `DMol3_TEST`
   - GAMESS: `GAMESS_TEST`
   - Gaussian09: `GAUSSIAN09_TEST/ethanol.cosmo -> ethanol.sigma`

## Manual ORCA regeneration

```bash
python ..\to_sigma.py --inpath ORCA_TEST/ethanol.cpcm --outpath ORCA_TEST/ethanol_regenerated.sigma --n 3 --averaging Hsieh
```

Then compare `ORCA_TEST/ethanol_regenerated.sigma` with `ORCA_TEST/ethanol.sigma`.
