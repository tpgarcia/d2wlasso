# d2wlasso: data-driven weighted lasso
The R package `d2wlasso` implements structured variable selection with q-values.

The reference for the original pliable lasso can be found at:
* [Structured variable selection with q-values](https://doi.org/10.1093/biostatistics/kxt012) by Tanya P. Garcia et al (2013).

## Installation

```
devtools::install_github("rakheon/d2wlasso", force = TRUE)
```

## Example

```
# data generation
x=matrix(rnorm(100*5, 0, 1),100,5)

# variable selection with d2wlasso


```
