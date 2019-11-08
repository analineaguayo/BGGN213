Class12\_InClass
================
Analine Aguayo
11/8/2019

Prepare protein structure for dosking
-------------------------------------

We want to download the 1HSG PDB structure and then produce a "protein-only" and "ligand-only" new separate PDB files.

``` r
library(bio3d)

get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

    ## [1] "./1hsg.pdb"
