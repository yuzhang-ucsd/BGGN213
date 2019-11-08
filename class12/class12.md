Class12: Structural Bioinformatics
================

## Prepare protein structure for docking

We want to download the 1HSG PDB structure and produce a “protein-only”
and “ligand-only” new separate PDB files.

``` r
library("bio3d")
get.pdb("1HSG")
```

    ## Warning in get.pdb("1HSG"): ./1HSG.pdb exists. Skipping download

    ## [1] "./1HSG.pdb"

Produce a “1hsg\_protein.pdb” and “1hsg\_ligand.pdb” file

``` r
pdb <- read.pdb("1HSG.pdb")
ligand <- atom.select(pdb, "ligand", value = T)
write.pdb(ligand, file = "1hsg_ligand.pdb")
protein <- atom.select(pdb, "protein", value = T)
write.pdb(protein, file = "1hsg_protein.pdb")
```
