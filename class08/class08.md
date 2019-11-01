class08
================
Yu Zhang
2019/10/25

## K-means

``` r
tmp <- c(rnorm(30,-3), rnorm(30,3))
x <- cbind(x=tmp, y=rev(tmp))
plot(x)
```

![](class08_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# k-means algorithm with 2 centers, run 20 times
k=kmeans(x, centers=2, nstart=20)
# sum of squares: between_SS / total_SS =  91.8 %
k$size
```

    ## [1] 30 30

``` r
k$cluster
```

    ##  [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1
    ## [36] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

``` r
k$centers
```

    ##           x         y
    ## 1  3.482720 -2.975203
    ## 2 -2.975203  3.482720

``` r
table(k$cluster)
```

    ## 
    ##  1  2 
    ## 30 30

``` r
#plot(x, col = "red") 
#or col = 3
plot(x, col = k$cluster)
points(k$centers, col = "blue", pch  = 15)
```

![](class08_files/figure-gfm/unnamed-chunk-1-2.png)<!-- --> \#\# HC

``` r
# First we need to calculate point (dis)similarity
# as the Euclidean distance between observations
dist_matrix <- dist(x)
# The hclust() function returns a hierarchical
# clustering model
hc <- hclust(d = dist_matrix)
# the print method is not so useful here
# or
# hc <- hclust(dist(x))
hc 
```

    ## 
    ## Call:
    ## hclust(d = dist_matrix)
    ## 
    ## Cluster method   : complete 
    ## Distance         : euclidean 
    ## Number of objects: 60

``` r
plot(hc)
abline(h = 6, col = "red")
```

![](class08_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
table(cutree(hc, h=4))
```

    ## 
    ##  1  2  3  4  5  6 
    ##  9 20  1 20  9  1

``` r
grps <- cutree(hc, k =2)
plot(x,col=grps)
```

![](class08_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

## 

``` r
# Step 1. Generate some example data for clustering
x <- rbind(
 matrix(rnorm(100, mean=0, sd = 0.3), ncol = 2), # c1
 matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2), # c2
 matrix(c(rnorm(50, mean = 1, sd = 0.3), # c3
 rnorm(50, mean = 0, sd = 0.3)), ncol = 2))
colnames(x) <- c("x", "y")
# Step 2. Plot the data without clustering
plot(x)
```

![](class08_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# Step 3. Generate colors for known clusters
# (just so we can compare to hclust results)
col <- as.factor( rep(c("c1","c2","c3"), each=50) )
plot(x, col=col)
```

![](class08_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
k2 <- kmeans(x,centers = 3, nstart = 20)
table(col,k2$cluster)
```

    ##     
    ## col   1  2  3
    ##   c1  2 45  3
    ##   c2 48  0  2
    ##   c3  3  0 47

``` r
hc2 <- hclust(dist(x))
plot(hc2)
```

![](class08_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
grps2 <- cutree(hc, k =7)
table(grps2)
```

    ## grps2
    ##  1  2  3  4  5  6  7 
    ##  9 20  1 17  9  3  1

## PCA Examples

``` r
mydata <- read.csv("https://tinyurl.com/expression-CSV",
 row.names=1)
head(mydata) 
```

    ##        wt1 wt2  wt3  wt4 wt5 ko1 ko2 ko3 ko4 ko5
    ## gene1  439 458  408  429 420  90  88  86  90  93
    ## gene2  219 200  204  210 187 427 423 434 433 426
    ## gene3 1006 989 1030 1017 973 252 237 238 226 210
    ## gene4  783 792  829  856 760 849 856 835 885 894
    ## gene5  181 249  204  244 225 277 305 272 270 279
    ## gene6  460 502  491  491 493 612 594 577 618 638

``` r
## lets do PCA
pca <- prcomp(t(mydata), scale=TRUE)
## See what is returned by the prcomp() function
attributes(pca) 
```

    ## $names
    ## [1] "sdev"     "rotation" "center"   "scale"    "x"       
    ## 
    ## $class
    ## [1] "prcomp"

``` r
## A basic PC1 vs PC2 2-D plot
plot(pca$x[,1], pca$x[,2]) 
```

![](class08_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
## Variance captured per PC
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
pca.var.per
```

    ##  [1] 92.6  2.3  1.1  1.1  0.8  0.7  0.6  0.4  0.4  0.0

``` r
summary(pca)
```

    ## Importance of components:
    ##                           PC1    PC2     PC3     PC4     PC5     PC6
    ## Standard deviation     9.6237 1.5198 1.05787 1.05203 0.88062 0.82545
    ## Proportion of Variance 0.9262 0.0231 0.01119 0.01107 0.00775 0.00681
    ## Cumulative Proportion  0.9262 0.9493 0.96045 0.97152 0.97928 0.98609
    ##                            PC7     PC8     PC9      PC10
    ## Standard deviation     0.80111 0.62065 0.60342 3.348e-15
    ## Proportion of Variance 0.00642 0.00385 0.00364 0.000e+00
    ## Cumulative Proportion  0.99251 0.99636 1.00000 1.000e+00

``` r
barplot(pca.var.per, main="Scree Plot",
 xlab="Principal Component", ylab="Percent Variation")
```

![](class08_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
colvec <- as.factor( substr( colnames(mydata), 1, 2) )
plot(pca$x[,1], pca$x[,2], col=colvec, pch=16,
 xlab=paste0("PC1 (", pca.var.per[1], "%)"),
 ylab=paste0("PC2 (", pca.var.per[2], "%)"))
identify(pca$x[,1], pca$x[,2], labels=colnames(mydata)) 
```

![](class08_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

    ## integer(0)

``` r
plot(pca$x[,1:2], col = c(rep("red",5),rep("blue",5)))
```

![](class08_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->
