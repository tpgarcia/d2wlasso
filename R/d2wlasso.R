#' Implement structured variables selection with q-values
#'
#' @param x (n by m) matrix of main covariates where m is the number of covariates and n is the sample size
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable
#' @param y (n by 1) matrix of response variable
#' @param ttest logical. If TRUE, p-value for each covariate is computed from the linear regression and this does not require normality of the covariates. If FALSE, p-value is computed as the p-value of the correlation coefficient. Default is FALSE.
#' @param q_method indicates the method for choosing optimal tuning parameter in the q-value computation as proposed in Storey and Tibshirani (2003). One of "bootstrap" or "smoother". Default is "smoother" (smoothing spline).
#' @param plots logical. If TRUE, figures are plotted. Default is FALSE.
#' @param pi0.true logical. If TRUE, the estimate of the true proportion of the null hypothesis is set to the value of pi0.val which is given by the user. If FALSE, the estimate of the true proportion of the null hypothesis is computed by bootstrap or smoothing spline. Default is FALSE.
#' @param pi0.val A user supplied estimate of the true proportion of the null hypothesis. Used only when pi0.true is TRUE. Default is 0.9.
#' @param wt The weights to be used for the weighted lasso. One of "one","adapt","q_cor" or "q_parcor". "one" gives no weight. "adapt" gives adaptive lasso weights, that is, the inverse of the absolute value of regression coefficients. "q_cor" gives weights set to q-values BEFORE taking into account diet. "q_parcor" gives weights set to q-values AFTER taking into account diet.
#' @param weight_fn The function applied to the weights for the weighted lasso. One of "identity","sqrt","inverse_abs","square". "identity" is the identity function, "sqrt" is the square root function, "inverse_abs" is the inverse of the absolute value and "square" is the square function. Not used if wt is set to "adapt". Default is "identity".
#' @param include.z logical. If TRUE, the additional covariate z is forced to be included in the model. Default is TRUE.
#' @param z.wt constant for forcing z in the model. If z is not included in the model even if include.z is TRUE, try different value. Default is 1000.
#' @param thresh.q logical. If TRUE, remove excessively small weights by using threshold to maintain stability. The threshold is set to 0.0001.
#' @param alpha indicates cut-off for q-values (thresholding). That is, the covariate with q-value less than this cut-off is included in the model.
#' @param alpha.bh indicates cut-off for Benjamini-Hochberg adjusted p-value (thresholding). That is, the covariate with BH-adjusted p-value less than this cut-off is included in the model.
#' @param delta Among the lasso solution path, the best descriptive model is the one which minimizes the loss function: (residual sum of squares)/(estimator of the model error variance) - (sample size) + delta*(number of predictors in the selected model). If delta = 2, this loss function is Mallows' Cp.
#' @param lasso.delta.cv.mult logical. If TRUE, we run vfold cross-validation to select optimal delta multiple times (ncv times). Default is FALSE.
#' @param vfold indicates the number of folds of the cross-validation for selecting delta.
#' @param ncv indicates the number of cross-validation runs for selecting delta.
#'
#' @return
#' \itemize{
#'    \item {qval:} {q-value as proposed in Storey and Tibshirani (2003)}
#'    \item {BH-pval:} {Benjamini-Hochberg adjusted p-value as proposed in Benjamini and Hochberg (1995)}
#'    \item {pval:} {p-value for each covariate}
#'    \item {out.cor:} {p-value for each covariate}
#'    \item {out.benhoch.cor:} {p-value for each covariate}
#'    \item {out.parcor:} {p-value for each covariate}
#'    \item {out.benhoch:} {p-value for each covariate}
#'    \item {out.w:} {p-value for each covariate}
#'    \item {alpha:} {p-value for each covariate}
#'    \item {alpha.bh:} {p-value for each covariate}
#'    \item {delta:} {p-value for each covariate}
#' }
#' @return out.cor
#' @return out.benhoch.cor
#' @return out.parcor
#' @return out.benhoch
#' @return qval
#' @export
#'
#' @examples
#' x=matrix(rnorm(100*5, 0, 1),100,5)
#' z <- matrix(rbinom(100, 1, 0.5),100,1)
#' y=matrix(z[1,] + 2*x[1,] - 2*x[2,] + rnorm(100, 0, 1), 100)
#' d2wlasso(x,z,y)
d2wlasso <- function(x,z,y,ttest=FALSE,q_method=c("bootstrap","smoother")[2],plots=FALSE,pi0.true=FALSE,pi0.val=0.9,
                     wt=c("one","adapt","q_cor","q_parcor")[4],weight_fn=c("identity","sqrt","inverse_abs","square")[1],
                     include.z=TRUE,z.wt=1000,thresh.q=TRUE,alpha=0.15,alpha.bh=0.05,delta=2,
                     lasso.delta.cv.mult=FALSE,vfold=10,ncv=100,percents.range=c(50,60,70,80,90,100),delta.cv.seed=NULL){

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

    ## Results for testing if a microbe has an effect on phenotype, but NOT
    ##            accounting for diet
    ## That is, we test : H_0 : \beta_{x_j}=0
    microbe.cor.out.qvalues <- q.computations(cor.out, method=q_method,
                                              plots=plots,file="cor",
                                              pi0.true=pi0.true,pi0.val=pi0.val)
    microbe.cor.out <- q.interest(microbe.cor.out.qvalues$qval.mat,alpha=alpha,criteria="less")
    out.cor   <- c(0,t(microbe.cor.out$interest))

    ## Results from Benjamini-Hochberg adjusted p-values when p-values do not account for diet
    benhoch.cor.results <- ben.hoch.interest(cor.out$pvalues,alpha=alpha.bh)
    out.benhoch.cor <- c(0,t(benhoch.cor.results$interest))

    ## Results for testing if a microbe has an effect on phenotype, but AFTER
    ##            accounting for diet
    ## That is, we test : H_0 : \beta_{x_j|z}=0
    # compute q-value as used by JD Storey with some adjustments made
    microbe.parcor.out.qvalues <- q.computations(parcor.out,method=q_method,
                                                 plots=plots,file="parcor",
                                                 pi0.true=pi0.true,pi0.val=pi0.val)
    out.qvalue <- c(0,t(microbe.parcor.out.qvalues$qval.mat))
    out.pvalue <- c(0,t(parcor.out$pvalues))
    q.out <- q.interest(microbe.parcor.out.qvalues$qval.mat,alpha=alpha,criteria="less")
    out.parcor <- c(1,t(q.out$interest))

    ## Results for Benjamini-Hochberg Method ##
    benhoch.results <- ben.hoch.interest(parcor.out$pvalues,alpha=alpha.bh)
    out.benhoch <- c(1,t(benhoch.results$interest))
    out.benhoch.pval.adjust <- c(0,t(benhoch.results$pval.adjust))

    if (wt == "one"){
        ## No weights
        weights <- matrix(1,nrow=nrow(microbes)-1,ncol=nrow(phenotypes))
    } else if (wt == "adapt"){
        ## Weights set to absolute value of partial correlations
        weights <- parcor.out$estimate
    } else if (wt == "q_cor"){
        ## Weights set to q-values BEFORE taking into account diet
        weights <- microbe.cor.out.qvalues$qval.mat
    } else {
        ## Weights set to q-values after taking into account diet
        weights <- microbe.parcor.out.qvalues$qval.mat
    }

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

    out.w <- as.data.frame(matrix(0,nrow=out.nrow,ncol=1,
                                  dimnames = list(out.rownames,paste("w.delta.",delta,sep=""))))

    if (wt == "adapt"){
        lasso.w <- lasso.computations(weights,microbes,phenotypes,g3,plots=plots,file="weight_",
                                      include.diet=include.z,diet.wt=z.wt,corr.g=TRUE,
                                      delta=delta)
    } else {
        lasso.w <- lasso.computations(weights,microbes,phenotypes,g,plots=plots,file="weight_",
                                      include.diet=include.z,diet.wt=z.wt,thresh.q=thresh.q,
                                      delta=delta)
    }
    out.w <- as.matrix(lasso.w$interest)

    ## mult.cv.delta.out.w5 : stores results from weighted lasso when weights are set to q-values AFTER taking into account diet,
    ##           and weight function g1

    nsimu = 1; j = 1

    mult.cv.delta.out.w5 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=nsimu,
                                                 dimnames = list(out.rownames,paste("w5.mult.nsimu.",seq(1,nsimu),sep=""))))

    mult.delta.w5 <- as.data.frame(matrix(0, nrow = 1, ncol = ncv,
                                          dimnames = list("delta", seq(1,ncv))))

    mult.cv.delta.out.w5.summary <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(percents.range),
                                                         dimnames = list(out.rownames,paste("w5.mult.cv.",percents.range,sep=""))))

    ## mult.cv.delta.out.w6 : stores results from weighted lasso when weights absolute value of partial correlations,
    ##           and weight function g3

    mult.cv.delta.out.w6 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=nsimu,
                                                 dimnames = list(out.rownames,paste("w6.mult.nsimu.",seq(1,nsimu),sep=""))))

    mult.delta.w6 <- as.data.frame(matrix(0, nrow = 1, ncol = ncv,
                                          dimnames = list("delta", seq(1,ncv))))


    mult.cv.delta.out.w6.summary <- as.data.frame(matrix(0,nrow=out.nrow,ncol=length(percents.range),
                                                         dimnames = list(out.rownames,paste("w6.mult.cv.",percents.range,sep=""))))

    if(lasso.delta.cv.mult==TRUE){
        if (!is.null(delta.cv.seed)){
            set.seed(delta.cv.seed)
        }
        include.diet <- TRUE

        ## Weights set to q-values after taking into account diet
        weights <- microbe.parcor.out.qvalues$qval.mat

        for(v in 1:ncv){
            mult.cv.delta.lasso.w5 <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,file="weight5_",
                                                         include.diet=include.diet,diet.wt=z.wt,thresh.q=thresh.q,delta=delta,
                                                         cv.criterion="delta_cv",vfold=vfold)
            mult.cv.delta.out.w5[,j] <- mult.cv.delta.out.w5[,j] + as.matrix(mult.cv.delta.lasso.w5$interest)
            mult.delta.w5[,v] <- mult.delta.w5[,v] + mult.cv.delta.lasso.w5$delta.out
        }

        ## Weights set to absolute value of partial correlations
        weights <- parcor.out$estimate

        for(v in 1:ncv){
            mult.cv.delta.lasso.w6 <- lasso.computations(weights,microbes,phenotypes,g3,plots=FALSE,file="weight6_",
                                                         include.diet=include.diet,diet.wt=z.wt,corr.g=TRUE,delta=delta,
                                                         cv.criterion="delta_cv",vfold=vfold)
            mult.cv.delta.out.w6[,j] <- mult.cv.delta.out.w6[,j] + as.matrix(mult.cv.delta.lasso.w6$interest)
            mult.delta.w6[,v] <- mult.delta.w6[,v] + mult.cv.delta.lasso.w6$delta.out
        }

        for(v in 1:length(percents.range)){

            ## Weights set to q-values after taking into account diet
            mult.cv.delta.out.w5.summary[,v] <- org.mult.cv.delta(mult.cv.delta.out.w5,percents.range[v],ncv)

            ## Weights set to absolute value of partial correlations
            mult.cv.delta.out.w6.summary[,v] <- org.mult.cv.delta(mult.cv.delta.out.w6,percents.range[v],ncv)
        }

    }

    return(list("qval"=out.qvalue,"BH-pval"=out.benhoch.pval.adjust, "pval"=out.pvalue, "out.cor"=out.cor, "out.benhoch.cor"=out.benhoch.cor, "out.parcor"=out.parcor, "out.benhoch"=out.benhoch, "out.w"=out.w, "alpha"=alpha, "alpha.bh"=alpha.bh, "delta"=delta, "mult.delta.w5"=mult.delta.w5, "mult.delta.w6"=mult.delta.w6, "mult.cv.delta.out.w5.summary"=mult.cv.delta.out.w5.summary, "mult.cv.delta.out.w6.summary"=mult.cv.delta.out.w6.summary, "mult.cv.delta.out.w5"=mult.cv.delta.out.w5, "mult.cv.delta.out.w6"=mult.cv.delta.out.w6))
}
