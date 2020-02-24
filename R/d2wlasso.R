#' d2wlasso package
#'
#' This package provides functions to perform variable selection with weighted lasso for both linear regression
#' and the Cox proportional hazards regression. The weights are chosen to direct the variable selection procedure
#' so that covariates that are highly
#' associated with the response are likely to be selected and covariates that are weakly associated with the response
#' are less likely to be selected.  Association between the response and the covariates is based on
#' results from simpler linear/Cox regressions between the response and each covariate, and include, for example,
#' q-values, partial correlation coefficients, t-statistics of regression coefficients, and exclusion frequency statistics.
#'
#'
#'@references
#' Garcia, T.P. and M¨uller, S. (2016). Cox regression with exclusion frequency-based weights to
#' identify neuroimaging markers relevant to Huntington’s disease onset. Annals of Applied Statistics, 10, 2130-2156.
#'
#' Garcia, T.P. and M¨uller, S. (2014). Influence of measures of significance-based weights in the weighted Lasso.
#' Journal of the Indian Society of Agricultural Statistics (Invited paper), 68, 131-144.
#'
#' Garcia, T.P., Mueller, S., Carroll, R.J., Dunn, T.N., Thomas, A.P., Adams, S.H., Pillai, S.D., and Walzem, R.L.
#' (2013). Structured variable selection with q-values. Biostatistics, DOI:10.1093/biostatistics/kxt012.
#'
#'
#' @docType package
#' @name d2wlasso
NULL
