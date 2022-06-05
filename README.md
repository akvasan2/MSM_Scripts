# MSM_Scripts

## General overview

Scripts used to construct an MSM for the conformational transitions of an outer membrane beta barrel.

## Description of scripts

1. Featurization.tcl:  Determine distance features from multiple simulations. In our case, this is the minimum distance between common hydrogen bonding residue pairs.
2. TICA.py: Determine optimal number of TICA eigenvectors to use for dimensionality reduction. Uses a scheme from Razavi et. al. to use TICs with non-Gaussian distribution.
3. VAMP2.py: Determine the optimal number of microstates to use when building MSM
4. ITS.py: Determines MSM lag time
5. MSM.py: Builds MSM using num of TICS/Microstates and lag time. Can save MSM and clusters for easy access in later commands
6. Save_data.py: Saves TICA and data information as .pickle files for quick loading
