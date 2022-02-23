# Replication Code for Instrumental Variable Treatment of Nonclassical Measurement Error Models [Hu and Schennach (2008)](https://doi.org/10.1111/j.0012-9682.2008.00823.x)
## Summary
This repository is my attempt at reproducing the estimation procedure from Yingyao Hu's and Susanne Schennach's paper titled Instrumental Variable Treatment of Nonclassical Measurement Error Models published in Econometrica (2008). The original code is written in GAUSS and translated to the R programming language. The estimator is a sieve maximum likelihood estimator, whose densities are conditioned on the unobserved covariates. I make a slight departure from their original implementation by using a Monte Carlo EM algorithm to integrate out the unobserved variable from the objective function. The maximization step of the algorithm estimates the densities corresponding to the response function, measurement error, and instrumental variable separately as opposed to jointly.

There are a number of computational improvements that need to be made before the estimator can be generalized to other models. Currently, the estimator performs well when each density is specified parametrically (e.g. normally distributed). Further updates to this repository will include the case for a nonparametric model, in which the densities can be approximated using the Method of Sieves (e.g. tensor product of Hermite polynomial basis functions). Doing so requires several computational constraints on model parameters to guarantee valid densities and a location normalization as required in their paper.




