# `LogConcDEAD`

**Computes a log-concave (maximum likelihood) estimator for i.i.d. data in any number of dimensions**

<hr>

Original Package Link: https://github.com/cran/LogConcDEAD

Deployed Pkgdown Website:
https://github.com/jhu-statprogramming-fall-2022/biostat840-project3-pkgdown-JungiinChoi

Original Package Author: Yining Chen

This is a list of customizations in my pkgdown website.

1. The bootswatch theme (by changing it to "cyborg")
2. The theme for syntax highlighting (by changing it to "arrow-dark")
3. Background, foreground and navbar&sidebar color (""#202123", "#B8BCC2" and "#306cc9" for each)
4. Reorganized the navbar without changing the default contents
5. Reorganized the footer (the authors’ information on the right, and puts advertising pkgdown on the left.)

<hr>

## Description

This package contains a function to compute the maximum likelihood estimator of a log-concave density in any number of dimensions using Shor's r-algorithm.

Functions to plot (for 1- and 2-d data), evaluate and draw samples from the maximum likelihood estimator are provided.

## Details

This package contains a selection of functions for maximum likelihood estimation under the constraint of log-concavity.

`mlelcd` computes the maximum likelihood estimator (specified via its value at data points). Output is a list of class "LogConcDEAD" which is used as input to various auxiliary functions.

`hatA` calculates the difference between the sample covariance and the fitted covariance.

`dlcd` evaluates the estimated density at a particular point.

`dslcd` evaluates the smoothed version of estimated density at a particular point.

`rlcd` draws samples from the estimated density.

`rslcd` draws samples from the smoothed version of estimated density.

`interplcd` interpolates the estimated density on a grid for plotting purposes.

`dmarglcd` evaluates the estimated marginal density by integrating the estimated density over an appropriate subspace.

`interpmarglcd` evaluates a marginal density estimate at equally spaced points along the axis for plotting purposes. This is done by integrating the estimated density over an appropriate subspace.

`plot.LogConcDEAD` produces plots of the maximum likelihood estimator, optionally using the rgl package.

`print` and `summary` methods are also available.

## Note

The authors gratefully acknowledge the assistance of Lutz Duembgen at the University of Bern for his insight into the objective function in `mlelcd`.

For one dimensional data, the active set algorithm in `logcondens` is much faster.


## Author(s)

Yining Chen (maintainer) y.chen101@lse.ac.uk

Madeleine Cule

Robert Gramacy

Richard Samworth

## Examples

```

## Some simple normal data, and a few plots

x <- matrix(rnorm(200),ncol=2)
lcd <- mlelcd(x)
g <- interplcd(lcd)
par(mfrow=c(2,2), ask=TRUE)
plot(lcd, g=g, type="c")
plot(lcd, g=g, type="c", uselog=TRUE)
plot(lcd, g=g, type="i")
plot(lcd, g=g, type="i", uselog=TRUE)

## Some plots of marginal estimates
par(mfrow=c(1,1))
g.marg1 <- interpmarglcd(lcd, marg=1)
g.marg2 <- interpmarglcd(lcd, marg=2)
plot(lcd, marg=1, g.marg=g.marg1)
plot(lcd, marg=2, g.marg=g.marg2) 

## generate some points from the fitted density
generated <- rlcd(100, lcd)
genmean <- colMeans(generated)

## evaluate the fitted density
mypoint <- c(0, 0)
dlcd(mypoint, lcd, uselog=FALSE)
mypoint <- c(10, 0)
dlcd(mypoint, lcd, uselog=FALSE)

## evaluate the marginal density
dmarglcd(0, lcd, marg=1)
dmarglcd(1, lcd, marg=2)
```