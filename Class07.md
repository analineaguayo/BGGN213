Class 7: R Functions and Packages
================
Analine Aguayo
10/23/2019

\#\#Revisit our functions from last day \# control alt I puts in an R
code

``` r
source("http://tinyurl.com/rescale-R")
```

Let’s try our rescale() function from last
    day

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale2(c(3, 10, NA, 7, 8))
```

    ## [1] 0.0000000 1.0000000        NA 0.5714286 0.7142857

\#\#Write a function both\_()na We want to write a function, called
both\_na(), that counts how many positions in two input vectors, x and
y, both have a missing value

``` r
#from stackOverFlow
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
which( is.na(x))
```

    ## [1] 3 5

Working snippet of code

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
both.na <-function (x,y) {sum(is.na(x) & is.na(y))}
```

``` r
both.na(x,y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)

both.na(x,y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
y3 <-c(1, NA, NA, NA, NA, NA)

both.na(x,y3)
```

    ## [1] 5

``` r
both.na2<-function (x,y) {
  if(length(x) != length(y)){
    stop("Input vectors should be the same length")
  }
  sum(is.na(x) & is.na(y))
  }
```

``` r
# student 1
c(100, 100, 100, 100, 100, 100, 100, 90)
```

    ## [1] 100 100 100 100 100 100 100  90

``` r
# student 2
c(100, NA, 90, 90, 90, 90, 97, 80)
```

    ## [1] 100  NA  90  90  90  90  97  80

``` r
x <- c(100, 100, 100, 100, 100, 100, 100, 90)
y <- c(100, NA, 90, 90, 90, 90, 97, 80)
xgrade <- sum(x[-which.min(x)], na.rm=TRUE) / (length(x)-1)
ygrade <- sum(y[-which.min(y)], na.rm=TRUE) / (length(y)-1)
xgrade
```

    ## [1] 100

``` r
ygrade
```

    ## [1] 79.57143

``` r
grades <- c(xgrade,ygrade)
grades
```

    ## [1] 100.00000  79.57143
