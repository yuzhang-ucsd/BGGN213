Untitled
================

## What is in the PDB database?

Download a CSV file from the PDB site (PDB ID: 1IEB)

``` r
file <- read.csv("data/Data Export Summary.csv")
```

> Q1: Determine the percentage of structures solved by X-Ray and
> Electron Microscopy.

``` r
ans <- file$Total/sum(file$Total)*100
names(ans) <- file$Experimental.Method
round(ans,2)
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##               89.06                8.13                2.51 
    ##               Other        Multi Method 
    ##                0.19                0.10

> Q2: Determine what proportion of structures are protein.

``` r
round(sum(file$Proteins)/sum(file$Total)*100,2)
```

    ## [1] 92.71
