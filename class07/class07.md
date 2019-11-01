class07
================
Yu Zhang
2019/10/23

## Revisit

# how to make a function read-only? like oo rescale

``` r
source("http://tinyurl.com/rescale-R")
```

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale(c(3,10,NA,8,2,3,4))
```

    ## [1] 0.125 1.000    NA 0.750 0.000 0.125 0.250

\#\#IMPORTANT\!\!\!

``` r
aa <- c(T,T,F) #to simplify TRUE and False
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
which(is.na(x))
```

    ## [1] 3 5

``` r
which(is.na(x) & is.na(y))
```

    ## [1] 3

``` r
both_na_location <- function(x, y) {
 which( is.na(x) & is.na(y) )
}
```

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA, 2, 2)
both_na <- function(x, y) {
 sum( is.na(x) & is.na(y) )
}
both_na(x,y2)
```

    ## [1] 3

``` r
both_na2 <- function(x, y) {
 if(length(x) != length(y)) {
 stop("not the same length")
 }
 sum( is.na(x) & is.na(y) )
}
```

## Grading

``` r
# student 1
s1 <- c(100, 100, 100, 100, 100, 100, 100, 90)
# student 2
s2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
min(s1)
```

    ## [1] 90

``` r
min(s2)
```

    ## [1] NA

``` r
grade <- function(s){
  return(mean(s[-which.min(s)], na.rm = TRUE))  ## VERY IMPROTANT
}
grade(s2)
```

    ## [1] 92.83333

``` r
#install.packages("bio3d")
library("bio3d")
#install.packages("BiocManager")
#install.packages("ggplot2")
#paletter
```

``` r
#install.packages("rgl")
library(rgl)
plot3d(rnorm(100), rnorm(100), rnorm(100))
abclines3d(0, 0, 0, a = diag(3), col = "gray")

open3d()
```

    ## wgl 
    ##   2

``` r
cols <- rainbow(7)
layout3d(matrix(1:16, 4,4), heights=c(1,3,1,3))
text3d(0,0,0,"tetrahedron3d"); next3d()
shade3d(tetrahedron3d(col=cols[1])); next3d()
text3d(0,0,0,"cube3d"); next3d()
shade3d(cube3d(col=cols[2])); next3d()
text3d(0,0,0,"octahedron3d"); next3d()
shade3d(octahedron3d(col=cols[3])); next3d()
text3d(0,0,0,"dodecahedron3d"); next3d()
shade3d(dodecahedron3d(col=cols[4])); next3d()
text3d(0,0,0,"icosahedron3d"); next3d()
shade3d(icosahedron3d(col=cols[5])); next3d()
text3d(0,0,0,"cuboctahedron3d"); next3d()
shade3d(cuboctahedron3d(col=cols[6])); next3d()
text3d(0,0,0,"oh3d"); next3d()
shade3d(oh3d(col=cols[7]))
```
