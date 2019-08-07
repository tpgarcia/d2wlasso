#' Implement structured variables selection with q-values
#'
#' @param x (n by m) matrix of main covariates where m is the number of covariates and n is the sample size
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
d2wlasso <- function(x,z,y,ttest=FALSE,method=c("bootstrap","smoother")[1],plots=FALSE,pi0.true=FALSE,pi0.val=0.9,
                     weight_fn=c("identity","sqrt","inverse_abs","square")[1],
                     include.diet=TRUE,diet.wt=1000,thresh.q=TRUE,delta.range=2){

    # dimension setting
    n <- nrow(x)
    m0 <- ncol(z)
    m <- ncol(x)

    # arranging input data
    X <- cbind(z, x)
    colnames(X) <- c("Diet",paste("X_",seq(1,m),sep=""))
    microbes <- t(X)
    phenotypes <- t(y)
    out.nrow <- nrow(microbes)
    out.rownames <- rownames(microbes)

    # compute correlation and partial correlation (for taking into account z) between x and y
    cor.out <- correlations(microbes,phenotypes,partial=FALSE,ttest=ttest)
    parcor.out <- correlations(microbes,phenotypes,partial=TRUE,ttest=ttest)

    # compute q-value as used by JD Storey with some adjustments made
    microbe.parcor.out.qvalues <- q.computations(parcor.out,method=method,
                                                 plots=FALSE,file="parcor",
                                                 pi0.true=pi0.true,pi0.val=pi0.val)
    out.qvalue <- c(0,t(microbe.parcor.out.qvalues$qval.mat))
    out.pvalue <- c(0,t(parcor.out$pvalues))

    ## Weights set to q-values after taking into account diet
    weights <- microbe.parcor.out.qvalues$qval.mat

    ## Weight functions
    g1 <- function(x){
        return(x)
    }
    g2 <- function(x){
        return(sqrt(x))
    }
    g3 <- function(x){
        return(1/abs(x))
    }
    g4 <- function(x){
        return(x^2)
    }
    if (weight_fn=="sqrt"){
        g <- g2
    } else if (weight_fn=="inverse_abs"){
        g <- g3
    } else if (weight_fn=="square"){
        g <- g4
    } else {
        g <- g1
    }

    ## Lasso calculations: With Diet forced in model ##
    out.w <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(delta.range),
                                   dimnames = list(out.rownames,paste("w.delta.",delta.range,sep=""))))
    print(out.w)
    for(d in 1:length(delta.range)){
        print(d)
        lasso.w <- lasso.computations(weights,microbes,phenotypes,g,plots=FALSE,file="weight_",
                                       include.diet=include.diet,diet.wt=diet.wt,thresh.q=thresh.q,
                                       delta=delta.range[d])
        out.w[,d] <- out.w[,d] + as.matrix(lasso.w$interest)
    }

    return(list("cor.out"=cor.out, "parcor.out"=parcor.out, "qval"=out.qvalue, "pval"=out.pvalue, "out.w"=out.w))
}
