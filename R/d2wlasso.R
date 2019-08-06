#' Implement structured variables selection with q-values
#'
#' @param x (n by m1) matrix of main covariates
#' @param z (n by m0) matrix of additional fixed covariates affecting response variable
#' @param y vector of response variable
#'
#' @return sumx: sum of main covariates
#' @return sumz: sum of additional fixed covariates
#' @export
#'
#' @examples
#' x=matrix(rnorm(100*5, 0, 1),100,5)
#' z=matrix(rnorm(100*1, 0, 1),100,1)
#' y=2*x[,1] - (2+2*z[,1])*x[,2] + rnorm(100, 0, 1)
#' d2wlasso(x,z,y)
d2wlasso <- function(x,z,y){
    print(x)
    return(list("sumx"=sum(x), "sumz"=sum(z)))
}
