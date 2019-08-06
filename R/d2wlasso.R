#' Implement structured variables selection with q-values
#'
#' @param x (n by m1) matrix of main covariates where m1 is the number of covariates and n is the sample size
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable
#' @param y (n by 1) matrix of response variable
#'
#' @return cor.out
#' @return parcor.out
#' @export
#'
#' @examples
#' x=matrix(rnorm(100*5, 0, 1),100,5)
#' z <- matrix(rbinom(100, 1, 0.5),100,1)
#' y=matrix(z[1,] + 2*x[1,] - 2*x[2,] + rnorm(100, 0, 1), 100)
#' d2wlasso(x,z,y)
d2wlasso <- function(x,z,y,ttest=FALSE){

    # dimension setting
    n <- nrow(x)
    m0 <- ncol(z)
    m1 <- ncol(x)

    # arranging input data
    X <- cbind(z, x)
    colnames(X) <- c("Diet",paste("X_",seq(1,m1),sep=""))
    microbes <- t(X)
    phenotypes <- t(y)

    # compute correlation and partial correlation (for taking into account z) between x and y
    cor.out <- correlations(microbes,phenotypes,partial=FALSE,ttest=ttest)
    parcor.out <- correlations(microbes,phenotypes,partial=TRUE,ttest=ttest)

    return(list("cor.out"=cor.out, "parcor.out"=parcor.out))
}
