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
d2wlasso <- function(x,z,y,ttest=TRUE,method=c("bootstrap","smoother")[2],plots=FALSE,pi0.true=FALSE,pi0.val=0.9,
                     wt=c("one","adapt","q_cor","q_parcor")[4],weight_fn=c("identity","sqrt","inverse_abs","square")[1],
                     include.z=TRUE,z.wt=1000,thresh.q=TRUE,alpha=0.15,delta=2,
                     vfold=10,lasso.delta.cv.mult=FALSE,delta.cv.seed=NULL,ncv=100,percents.range=c(50,60,70,80,90,100)){

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
    microbe.cor.out.qvalues <- q.computations(cor.out, method=method,
                                              plots=FALSE,file="cor",
                                              pi0.true=pi0.true,pi0.val=pi0.val)
    microbe.cor.out <- q.interest(microbe.cor.out.qvalues$qval.mat,alpha=alpha,criteria="less")
    out.cor   <- c(0,t(microbe.cor.out$interest))

    ## Results from Benjamini-Hochberg adjusted p-values when p-values do not account for diet
    benhoch.cor.results <- ben.hoch.interest(cor.out$pvalues,alpha=alpha)
    out.benhoch.cor <- c(0,t(benhoch.cor.results$interest))

    ## Results for testing if a microbe has an effect on phenotype, but AFTER
    ##            accounting for diet
    ## That is, we test : H_0 : \beta_{x_j|z}=0
    # compute q-value as used by JD Storey with some adjustments made
    microbe.parcor.out.qvalues <- q.computations(parcor.out,method=method,
                                                 plots=FALSE,file="parcor",
                                                 pi0.true=pi0.true,pi0.val=pi0.val)
    out.qvalue <- c(0,t(microbe.parcor.out.qvalues$qval.mat))
    out.pvalue <- c(0,t(parcor.out$pvalues))
    q.out <- q.interest(microbe.parcor.out.qvalues$qval.mat,alpha=alpha,criteria="less")
    out.parcor <- c(1,t(q.out$interest))

    ## Results for Benjamini-Hochberg Method ##
    benhoch.results <- ben.hoch.interest(parcor.out$pvalues,alpha=alpha)
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
        lasso.w <- lasso.computations(weights,microbes,phenotypes,g3,plots=FALSE,file="weight_",
                                      include.diet=include.z,diet.wt=z.wt,corr.g=TRUE,
                                      delta=delta)
    } else {
        lasso.w <- lasso.computations(weights,microbes,phenotypes,g,plots=FALSE,file="weight_",
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

    return(list("cor.out"=cor.out, "parcor.out"=parcor.out, "qval"=out.qvalue,"out.benhoch.pval.adjust"=out.benhoch.pval.adjust, "pval"=out.pvalue, "out.cor"=out.cor, "out.benhoch.cor"=out.benhoch.cor, "out.parcor"=out.parcor, "out.benhoch"=out.benhoch, "out.w"=out.w, "alpha"=alpha, "delta"=delta, "mult.delta.w5"=mult.delta.w5, "mult.delta.w6"=mult.delta.w6, "mult.cv.delta.out.w5.summary"=mult.cv.delta.out.w5.summary, "mult.cv.delta.out.w6.summary"=mult.cv.delta.out.w6.summary, "mult.cv.delta.out.w5"=mult.cv.delta.out.w5, "mult.cv.delta.out.w6"=mult.cv.delta.out.w6))
}
