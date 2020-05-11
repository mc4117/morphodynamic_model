Morphodynamic Thetis model
================

This repository contains the model described in the paper

*Mariana C. A. Clare, James Percival, Athanasios Angeloudis, Colin J. Cotter and Matthew D. Piggott*, **Hydro-morphodynamics 2D modelling using a discontinuous Galerkin discretisation**, Computers & Geosciences.


Software requirements
-------------------------

1. Firedrake (www.firedrakeproject.org)
    * The version of the Firedrake model we use has been stored and can be downloaded from the following site: https://doi.org/10.5281/zenodo.3385061.
1. Thetis (https://thetisproject.org/download.html)
    * The version of the Thetis model we use has been stored and can be downloaded from the following site: https://zenodo.org/record/3385804    
3. Python 3.5 or later



Simulation scripts
------------------

* Section 4: Migrating Trench Test Case
    
    Reproduce Thetis results with:

```
#!bash
    $ python trench.py
```

    The plot created by this script is Figure 5 in the paper.


* Section 5: Meander Test Case
    
    Reproduce Thetis results with:
```
#!bash
    $ python meander.py
```

    The two plots created by this script are Figures 14a and 14b in the paper.

Additional notes
-------------------------

To verify our hydro-morphodynamics model in Thetis, we used the Telemac-Mascaret model (http://opentelemac.com/index.php/download).

The steering files for the Migrating Trench Test Case were obtained via private communication from the authors of 
*Villaret, Catherine, et al. "First-order uncertainty analysis using Algorithmic Differentiation of morphodynamic models." Computers & Geosciences 90 (2016): 144-151.*

The steering files for the Meander Test case are those for the "Yen" validation example in Sisyphe. The steering files come as part of the Telemac-Mascaret download.
