# d2wlasso: data-driven weighted lasso
The R package `d2wlasso` implements structured variable selection with q-values.

The reference for the original pliable lasso can be found at:
* [Structured variable selection with q-values](https://doi.org/10.1093/biostatistics/kxt012) by Tanya P. Garcia et al (2013).
* [Influence of Measures of Significance based Weights in the Weighted Lasso](https://www.statindex.org/articles/285259) by Tanya P. Garcia and Samuel Müller (2014).

## Installation

```
devtools::install_github("rakheon/d2wlasso", force = TRUE)
```

## Example

```
# data generation
x = matrix(rnorm(100*5, 0, 1),100,5)
z = matrix(rbinom(100, 1, 0.5),100,1)
y = matrix(z[,1] + 2*x[,1] - 2*x[,2] + rnorm(100, 0, 1), 100)

# variable selection with d2wlasso
dwl0 <- d2wlasso(x,z,y)
dwl1 <- d2wlasso(x,z,y,delta=2)
dwl2 <- d2wlasso(x,z,y,include.z=FALSE,delta=2)
dwl3 <- d2wlasso(x,z,y,weight_fn = "sqrt")
dwl4 <- d2wlasso(x,z,y,wt="adapt")
dwl5 <- d2wlasso(x,z,y,wt="t_val")
dwl6 <- d2wlasso(x,z,y,wt="q_parcor")

# select delta with cross-validation
dwlcv0 <- d2wlasso(x,z,y,lasso.delta.cv.mult = TRUE, ncv = 3)
dwlcv1 <- d2wlasso(x,z,y,lasso.delta.cv.mult = TRUE, ncv = 3, delta.cv.seed = 1)
dwlcv2 <- d2wlasso(x,z,y,weight_fn = "square",lasso.delta.cv.mult = TRUE, ncv = 3, delta.cv.seed = 1)


```
