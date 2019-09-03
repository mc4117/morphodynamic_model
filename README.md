Morphodynamic Thetis model
================

This repository contains the model described in the paper

*Mariana C. A. Clare, James Percival, Athanasios Angeloudis, Colin J. Cotter and Matthew D. Piggott*, **Hydro-morphodynamics 2D modelling using a discontinuous Galerkin discretisation**, Computers & Geosciences.


Software requirements
-------------------------

1. Firedrake (www.firedrakeproject.org)
    * The version of the Firedrake model we use has been stored and can be downloaded from the following site: https://doi.org/10.5281/zenodo.3385061.
1. Thetis (https://thetisproject.org/download.html)
    * We use the Thetis model at the following commit id: bab2234b4f84de6c3e5244ca0e466ec15a62d991, available at: https://github.com/thetisproject/thetis/commit/bab2234b4f84de6c3e5244ca0e466ec15a62d991
    
3. Python 3.5 or later



Simulation scripts
------------------

* Section 4: Migrating Trench Test Case
    
    Reproduce with:

```
#!bash
    $ python trench.py
```


* Section 5: Meander Test Case
    
    Reproduce with:
```
#!bash
    $ python meander.py
```
