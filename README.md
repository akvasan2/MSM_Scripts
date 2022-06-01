# MSM_Scripts

## General overview

Scripts used to construct an MSM for the conformational transitions of an outer membrane beta barrel.

## Description of scripts

1. Featurization.tcl:  Determine distance features from multiple simulations. In our case, this is the minimum distance between common hydrogen bonding residue pairs.
2. TICA.py: Determine optimal number of TICA eigenvectors to use for dimensionality reduction. Uses a scheme from Razavi et. al. to use TICs with non-Gaussian distribution.
3. 
