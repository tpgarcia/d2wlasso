#' Implement structured variables selection with q-values
#'
#' @param x (n by m) matrix of main covariates where m is the number of covariates and n is the sample size
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable
#' @param y (n by 1) matrix of response variable
#' @param cox.delta (n by 1) matrix of status for survival analysis
#' @param factor.z logical. If TRUE, the additional fixed variable z is used as factor
#' @param reg.type indicates the model for fitting. Either "linear" or "cox". Default is "linear".
#' @param ttest logical. If TRUE, p-value for each covariate is computed from the linear regression and this does not require normality of the covariates. If FALSE, p-value is computed as the p-value of the correlation coefficient. Default is FALSE.
#' @param q_method indicates the method for choosing optimal tuning parameter in the q-value computation as proposed in Storey and Tibshirani (2003). One of "bootstrap" or "smoother". Default is "smoother" (smoothing spline).
#' @param plots logical. If TRUE, figures are plotted. Default is FALSE.
#' @param pi0.true logical. If TRUE, the estimate of the true proportion of the null hypothesis is set to the value of pi0.val which is given by the user. If FALSE, the estimate of the true proportion of the null hypothesis is computed by bootstrap or smoothing spline. Default is FALSE.
#' @param pi0.val A user supplied estimate of the true proportion of the null hypothesis. Used only when pi0.true is TRUE. Default is 0.9.
#' @param wt The weights to be used for the weighted lasso. One of "one","t_val","parcor","p_val","bhp_val","adapt","q_cor" or "q_parcor". "one" gives no weight. "t_val" gives weight of the inverse absolute t-statistics of the regression coefficients. "parcor" gives weight of the inverse absolute partial correlation between the main covariate and the response after accounting for z. "p_val" gives p-value of each predictor's coefficient as weights. "bhp_val" gives Benjamini-Hochberg adjusted p-value of each predictor's coefficient as weights. "adapt" gives adaptive lasso weights, that is, the inverse of the absolute value of regression coefficients. "q_cor" gives weights set to q-values BEFORE taking into account diet. "q_parcor" gives weights set to q-values AFTER taking into account diet.
#' @param weight_fn The function applied to the weights for the weighted lasso. One of "identity","sqrt","inverse_abs","square". "identity" is the identity function, "sqrt" is the square root function, "inverse_abs" is the inverse of the absolute value and "square" is the square function. Not used if wt is set to "adapt". Default is "identity".
#' @param include.z logical. If TRUE, the additional covariate z is forced to be included in the model. Default is TRUE.
#' @param z.wt constant for forcing z in the model. If z is not included in the model even if include.z is TRUE, try different value. Default is 1000.
#' @param thresh.q logical. If TRUE, remove excessively small weights by using threshold to maintain stability. The threshold is set to 0.0001.
#' @param alpha indicates cut-off for q-values (thresholding). That is, the covariate with q-value less than this cut-off is included in the model.
#' @param alpha.bh indicates cut-off for Benjamini-Hochberg adjusted p-value (thresholding). That is, the covariate with BH-adjusted p-value less than this cut-off is included in the model.
#' @param delta Among the lasso solution path, the best descriptive model is the one which minimizes the loss function: (residual sum of squares)/(estimator of the model error variance) - (sample size) + delta*(number of predictors in the selected model). If delta = 2, this loss function is Mallows' Cp.
#' @param robust indicates whether it is desired to make the estimate more robust for small p-values.
#' @param lasso.delta.cv.mult logical. If TRUE, we run vfold cross-validation to select optimal delta multiple times (ncv times). Default is FALSE.
#' @param vfold indicates the number of folds of the cross-validation for selecting delta.
#' @param ncv indicates the number of cross-validation runs for selecting delta.
#' @param delta.cv.seed For reproducible cross-validation result, the seed can be fixed.
#' @param run.aic.bic If TRUE, the Cox regression with exclusion frequency-based weights is performed with randomly partitioning the index
#' @param run.kmeans.aic.bic If TRUE, the Cox regression with exclusion frequency-based weights is performed with partitioning the index using k-means
#' @param run.kquart.aic.bic If TRUE, the Cox regression with exclusion frequency-based weights is performed with partitioning the index using k-quartile
#' @param run.sort.aic.bic If TRUE, the Cox regression with exclusion frequency-based weights is performed with partitioning the index using sorted partition
#' @param nboot indicates the number of bootstrap samples for the Cox regression with exclusion frequency-based weights
#' @param k indicates the number of partitions for the Cox regression with exclusion frequency-based weights. Default is 4.
#' @param direction indicates the direction of stepwise regression for the Cox regression with exclusion frequency-based weights. One of "both", "forward" or "backward". Default is "backward".
#'
#' @return
#' \itemize{
#'    \item \strong{qval:} {q-value as proposed in Storey and Tibshirani (2003)}
#'    \item \strong{bh.pval:} {Benjamini-Hochberg adjusted p-value as proposed in Benjamini and Hochberg (1995)}
#'    \item \strong{pval:} {p-value for each covariate}
#'    \item \strong{out.cor:} {variable selection results for testing if a main covariate has an effect on the response variable, but NOT accounting for the additional fixed covariate z}
#'    \item \strong{out.parcor:} {variable selection results for testing if a main covariate has an effect on the response variable, but AFTER accounting for the additional fixed covariate z}
#'    \item \strong{out.benhoch.cor:} {variable selection results from Benjamini-Hochberg adjusted p-values when p-values do not account for the additional fixed covariate z}
#'    \item \strong{out.benhoch.parcor:} {variable selection results from Benjamini-Hochberg adjusted p-values when p-values account for the additional fixed covariate z}
#'    \item \strong{out.w:} {variable selection results from weighted lasso}
#'    \item \strong{alpha:} {level of significance to compare with the q-values}
#'    \item \strong{alpha.bh:} {level of significance to compare with the Benjamini-Hochberg adjusted p-values}
#'    \item \strong{delta:} {the multiplier to the number of predictors in the penalized loss function for variable selection. The loss function is defined as: \code{ SSE/(sigma^2) - n + delta*p } where SSE denotes the residual sum of squares, sigma denotes the estimator of the model error variance, n is the sample size and p is the number of predictors in the selected model. }
#'    \item \strong{cv.delta.w:} {the selected delta from the cross-validation for the weighted lasso with q-values}
#'    \item \strong{cv.delta.adapt:} {the selected delta from the cross-validation for the adaptive lasso}
#'    \item \strong{cv.out.w:} {the aggregated result of the cross-validation for the weighted lasso with q-values}
#'    \item \strong{cv.out.adapt:} {the aggregated result of the cross-validation for the adaptive lasso}
#'    \item \strong{w.aic.boot:} {variable selection results from AIC and the Cox regression with exclusion frequency-based weights from random partitioning}
#'    \item \strong{w.bic.boot:} {variable selection results from BIC and the Cox regression with exclusion frequency-based weights from random partitioning}
#'    \item \strong{w.kmeans.aic.boot:} {variable selection results from AIC and the Cox regression with exclusion frequency-based weights from k-means partitioning}
#'    \item \strong{w.kmeans.bic.boot:} {variable selection results from BIC and the Cox regression with exclusion frequency-based weights from k-means partitioning}
#'    \item \strong{w.kquart.aic.boot:} {variable selection results from AIC and the Cox regression with exclusion frequency-based weights from k-quartile partitioning}
#'    \item \strong{w.kquart.bic.boot:} {variable selection results from BIC and the Cox regression with exclusion frequency-based weights from k-quartile partitioning}
#'    \item \strong{w.sort.aic.boot:} {variable selection results from AIC and the Cox regression with exclusion frequency-based weights from sorted partitioning}
#'    \item \strong{w.sort.bic.boot:} {variable selection results from BIC and the Cox regression with exclusion frequency-based weights from sorted partitioning}
#' }
#' @export
#'
#' @examples
#' x = matrix(rnorm(100*5, 0, 1),100,5)
#' z = matrix(rbinom(100, 1, 0.5),100,1)
#' y = matrix(z[,1] + 2*x[,1] - 2*x[,2] + rnorm(100, 0, 1), 100)
#' dwl0 <- d2wlasso(x,z,y)
#' dwl1 <- d2wlasso(x,z,y,delta=2)
#' dwl2 <- d2wlasso(x,z,y,include.z=FALSE)
#' dwl3 <- d2wlasso(x,z,y,weight_fn = "sqrt")
#' dwl4 <- d2wlasso(x,z,y,wt="adapt")
#' dwl5 <- d2wlasso(x,z,y,wt="t_val")
#' dwl6 <- d2wlasso(x,z,y,wt="q_parcor")
#' dwlcv0 <- d2wlasso(x,z,y,lasso.delta.cv.mult = TRUE, ncv = 3)
#' dwlcv1 <- d2wlasso(x,z,y,lasso.delta.cv.mult = TRUE, ncv = 3, delta.cv.seed = 1)
#' dwlcv2 <- d2wlasso(x,z,y,weight_fn = "square",lasso.delta.cv.mult = TRUE, ncv = 3, delta.cv.seed = 1)
#'
#' x <- matrix(rnorm(100*5, 0, 1),100,5)
#' z <- matrix(rbinom(100, 1, 0.5),100,1)
#' y <- matrix(exp(z[,1] + 2*x[,1] - 2*x[,2] + rnorm(100, 0, 2)), 100)
#' cox.delta <- matrix(1,nrow=length(y),ncol=1)
#' dwlcox1 <- d2wlasso(x,z,y,cox.delta = cox.delta, reg.type = "cox")
#' dwlcox2 <- d2wlasso(x,z,y,cox.delta = cox.delta, reg.type = "cox", nboot = 50)
#' dwlcox3 <- d2wlasso(x,z,y,cox.delta = cox.delta, reg.type = "cox", wt="t_val")
#' dwlcoxcv1 <- d2wlasso(x,z,y,cox.delta = cox.delta,reg.type = "cox",lasso.delta.cv.mult = TRUE, ncv = 3, nboot = 50)
d2wlasso <- function(x,z,y,cox.delta=NULL,factor.z=TRUE,reg.type=c("linear","cox")[1],ttest=TRUE,q_method=c("bootstrap","smoother")[2],plots=FALSE,pi0.true=FALSE,pi0.val=0.9,
                     wt=c("one","t_val","parcor","p_val","bhp_val","adapt","q_cor","q_parcor")[7],weight_fn=c("identity","sqrt","inverse_abs","square")[1],
                     include.z=TRUE,z.wt=1000,thresh.q=TRUE,alpha=0.15,alpha.bh=0.05,delta=2,robust=TRUE,q.old=FALSE,
                     lasso.delta.cv.mult=FALSE,vfold=10,ncv=100,delta.cv.seed=NULL,#cv.criterion=FALSE,
                     run.aic.bic=TRUE,
                     #run.fixed.aic.bic=TRUE,
                     run.kmeans.aic.bic=TRUE,
                     run.kquart.aic.bic=TRUE,
                     run.sort.aic.bic=TRUE,
                     #run.aic=TRUE,
                     #run.bic=TRUE,
                     nboot=100,k=4,#k.split=k,
                     direction="backward"){
    real_data = TRUE
    k.split=k
    run.aic = TRUE
    run.bic = TRUE
    cv.criterion=FALSE

    # dimension setting
    n <- nrow(x)
    m0 <- ncol(z)
    m <- ncol(x)

    # arranging input data
    X <- cbind(z, x)
    colnames(X) <- c("Fixed",paste("X_",seq(1,m),sep=""))
    microbes <- t(X)
    if (reg.type=="linear"){
        phenotypes <- list(yy=t(y),delta=NULL)
    } else if (length(cox.delta)!=length(y)){
        stop("length of y should be equal to length of cox.delta")
    } else {
        phenotypes <- list(yy=t(y),delta=t(cox.delta))
    }
    out.nrow <- nrow(microbes)
    out.rownames <- rownames(microbes)
    #print("here")

    # compute correlation and partial correlation (for taking into account z) between x and y
    cor.out <- correlations(factor.z,microbes,phenotypes,partial=FALSE,ttest=ttest,format.data=FALSE,reg.type=reg.type)
    parcor.out <- correlations(factor.z,microbes,phenotypes,partial=TRUE,ttest=ttest,format.data=FALSE,reg.type=reg.type)
    fstat.out <- ftests(factor.z,microbes)
    #print(cor.out)

    #####################
    ## Get Weights: t-statistic, p-values, Benjamin-Hochberg p-values, q-values, partial correlations ##
    #####################

    ## t-values with \beta_k in y=diet + \beta_k * microbe_k
    weight.tvalue <- parcor.out$tvalues
    #out.tvalue[,j] <- parcor.out$tvalues

    ## p-values with \beta_k in y=diet + \beta_k * microbe_k
    weight.pvalue.noadj <- parcor.out$pvalues
    #out.pvalue.noadj[,j] <- parcor.out$pvalues

    ## Benjamini-Hochberg adjusted p-values from y=diet + \beta_k * microbe_k
    benhoch.results <- ben.hoch.interest(parcor.out$pvalues,alpha=alpha)
    weight.pvalue.benhoch <- benhoch.results$pval.adjust
    #out.pvalue.benhoch[,j] <- benhoch.results$pval.adjust

    ## partial correlations
    weight.parcor <- parcor.out$estimate
    #out.parcor.value[,j] <- parcor.out$estimate

    ## Results for testing if a microbe has an effect on phenotype, but NOT
    ##            accounting for diet
    ## That is, we test : H_0 : \beta_{x_j}=0
    microbe.cor.out.qvalues <- q.computations(cor.out, method=q_method,
                                              plots=plots,file="cor",robust=robust,q.old=q.old,
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
                                                 plots=plots,file="parcor",robust=robust,q.old=q.old,
                                                 pi0.true=pi0.true,pi0.val=pi0.val)
    out.qvalue <- c(0,t(microbe.parcor.out.qvalues$qval.mat))
    out.pvalue <- c(0,t(parcor.out$pvalues))
    q.out <- q.interest(microbe.parcor.out.qvalues$qval.mat,alpha=alpha,criteria="less")
    out.parcor <- c(1,t(q.out$interest))

    ## Results for Benjamini-Hochberg Method ##
    benhoch.results <- ben.hoch.interest(parcor.out$pvalues,alpha=alpha.bh)
    out.benhoch <- c(1,t(benhoch.results$interest))
    out.benhoch.pval.adjust <- c(0,t(benhoch.results$pval.adjust))

    ## Results for Cox with Exclusion-frequency weights ##
    out.aic.boot <- NULL; out.bic.boot <- NULL#; out.fixed.aic.boot <- NULL; out.fixed.bic.boot <- NULL
    out.kmeans.aic.boot <- NULL; out.kmeans.bic.boot <- NULL; out.kquart.aic.boot <- NULL; out.kquart.bic.boot <- NULL
    out.sort.aic.boot <- NULL; out.sort.bic.boot <- NULL
    w.aic.boot <- NULL; w.bic.boot <- NULL#; w.fixed.aic.boot <- NULL; w.fixed.bic.boot <- NULL
    w.kmeans.aic.boot <- NULL; w.kmeans.bic.boot <- NULL; w.kquart.aic.boot <- NULL; w.kquart.bic.boot <- NULL
    w.sort.aic.boot <- NULL; w.sort.bic.boot <- NULL

    if (wt == "one"){
        ## No weights
        weights <- matrix(1,nrow=nrow(microbes)-1,ncol=nrow(phenotypes))
    } else if (wt == "t_val"){
        weights <- weight.tvalue
    } else if (wt == "parcor"){
        weights <- weight.parcor
    } else if (wt == "p_val"){
        weights <- weight.pvalue.noadj
    } else if (wt == "bhp_val"){
        weights <- weight.pvalue.benhoch
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

    if ((wt == "t_val")||(wt == "parcor")||(wt == "adapt")){
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

    ## mult.cv.delta.out.w6 : stores results from weighted lasso when weights absolute value of partial correlations,
    ##           and weight function g3

    mult.cv.delta.out.w6 <- as.data.frame(matrix(0,nrow=out.nrow,ncol=nsimu,
                                                 dimnames = list(out.rownames,paste("w6.mult.nsimu.",seq(1,nsimu),sep=""))))

    mult.delta.w6 <- as.data.frame(matrix(0, nrow = 1, ncol = ncv,
                                          dimnames = list("delta", seq(1,ncv))))

    if(lasso.delta.cv.mult==TRUE){

        include.diet <- TRUE

        ## Weights set to q-values after taking into account diet
        weights <- microbe.parcor.out.qvalues$qval.mat
        if (!is.null(delta.cv.seed)){
            set.seed(delta.cv.seed)
        }

        for(v in 1:ncv){
            mult.cv.delta.lasso.w5 <- lasso.computations(weights,microbes,phenotypes,g,plots=FALSE,file="weight5_",
                                                         include.diet=include.diet,diet.wt=z.wt,thresh.q=thresh.q,delta=delta,
                                                         cv.criterion=FALSE,vfold=vfold)
            mult.cv.delta.out.w5[,j] <- mult.cv.delta.out.w5[,j] + as.matrix(mult.cv.delta.lasso.w5$interest)
            mult.delta.w5[,v] <- mult.delta.w5[,v] + mult.cv.delta.lasso.w5$delta.out
        }

        ## Weights set to absolute value of partial correlations
        weights <- parcor.out$estimate
        if (!is.null(delta.cv.seed)){
            set.seed(delta.cv.seed)
        }

        for(v in 1:ncv){
            mult.cv.delta.lasso.w6 <- lasso.computations(weights,microbes,phenotypes,g3,plots=FALSE,file="weight6_",
                                                         include.diet=include.diet,diet.wt=z.wt,corr.g=TRUE,delta=delta,
                                                         cv.criterion=FALSE,vfold=vfold)
            mult.cv.delta.out.w6[,j] <- mult.cv.delta.out.w6[,j] + as.matrix(mult.cv.delta.lasso.w6$interest)
            mult.delta.w6[,v] <- mult.delta.w6[,v] + mult.cv.delta.lasso.w6$delta.out
        }

    }

    ## exclusion frequency weights

    if (reg.type=="cox"){

    ## Store weights for each bootstrap
    tmp2.store <- as.data.frame(matrix(0, nrow = (out.nrow-1), ncol = nboot,
                                       dimnames = list(out.rownames[-1],
                                                       seq(1,nboot))))

    weight.aic.boot <- tmp2.store
    weight.bic.boot <- tmp2.store

    weight.fixed.aic.boot <- tmp2.store
    weight.fixed.bic.boot <- tmp2.store

    weight.kmeans.aic.boot <- tmp2.store
    weight.kmeans.bic.boot <- tmp2.store

    weight.kquart.aic.boot <- tmp2.store
    weight.kquart.bic.boot <- tmp2.store

    weight.sort.aic.boot <- tmp2.store
    weight.sort.bic.boot <- tmp2.store

    if(run.kmeans.aic.bic==TRUE | run.kquart.aic.bic==TRUE |  run.sort.aic.bic==TRUE){

        ## Get parameter estimates from ridge regression
        beta.values <- ridge.regression(microbes,phenotypes)

        if(run.kmeans.aic.bic==TRUE){
            ## K-means clustering
            kmeans.out <- kmeans(beta.values,centers=k.split,iter.max=100)
            index.group.kmeans <- kmeans.out$cluster
        }

        if(run.kquart.aic.bic==TRUE){

            ## K-quart clustering
            index.group.kquart <- cut(beta.values, breaks=quantile(beta.values,
                                                                   probs=seq(0,1, by=1/k.split)),include.lowest=TRUE)
            index.group.kquart <- as.numeric(factor(index.group.kquart, labels=1:k.split))
        }

        if(run.sort.aic.bic==TRUE){

            ## K-sort clustering
            beta.sort <- sort(abs(beta.values),decreasing=TRUE,index.return=TRUE)
            sort.beta.index <- beta.sort$ix

            ## index of ordering
            index.group.sort <- index.sort.partition(n=ncol(microbes),k=k,sort.beta.index)
        }
    }

    for(b in 1:nboot){
        ##print(b)
        if(run.aic.bic==TRUE){
            ## Randomly partition the index
            rand.index <- random.partition(n=ncol(microbes),p=nrow(microbes)-1,k=k)
        }

        #if(run.fixed.aic.bic==TRUE){
            ## Ensure fixed covariates are in the partition + randomly partition the rest
        #    rand.fixed.index <- fixed.plus.random.partition(fixed.covariates,n=ncol(microbes),p=nrow(microbes)-1,k=k)
        #}

        if(run.kmeans.aic.bic==TRUE){
            ## Partition the index using k-means
            kmeans.rand.index <- designed.partition(index.group.kmeans,k=k)
        }

        if(run.kquart.aic.bic==TRUE){
            ## Partition the index using k-quartile
            kquart.rand.index <- designed.partition(index.group.kquart,k=k)
        }

        if(run.sort.aic.bic==TRUE){
            ## Partition the index using k-quartile
            sort.rand.index <- designed.partition(index.group.sort,k=k)

            ## Measures how often the largest beta from ridge regression is in the true,
            ##  non-zero beta coefficients
            #beta.index.sort[sort.beta.index[1]+1,j] <- as.numeric(sort.beta.index[1]%in%fixed.covariates)
        }

        ## Apply stepwise AIC to each group
        for(l in 1:k){
            ##print(l)
            if(run.aic.bic==TRUE){
                ## Random partitioning
                index <- as.numeric(unlist(rand.index[l]))
                if(length(index)!=0){

                    if(run.aic==TRUE){
                        weight.aic.boot[,b] <- weight.aic.boot[,b] +
                            step.selection(factor.z,index,microbes,phenotypes,type="AIC",
                                           direction=direction,
                                           real_data=real_data)
                    }

                    if(run.bic==TRUE){
                        weight.bic.boot[,b] <- weight.bic.boot[,b] +
                            step.selection(factor.z,index,microbes,phenotypes,type="BIC",
                                           direction=direction,
                                           real_data=real_data)
                    }
                }
            }
            if (FALSE){
            if(run.fixed.aic.bic==TRUE){
                ## Fixed + Random partitioning
                index <- as.numeric(unlist(rand.fixed.index[l]))
                if(length(index)!=0){
                    if(run.aic==TRUE){
                        weight.fixed.aic.boot[,b] <- weight.fixed.aic.boot[,b] +
                            step.selection(factor.z,index,
                                           microbes,phenotypes,type="AIC",
                                           direction=direction,
                                           real_data=real_data)
                    }

                    if(run.bic==TRUE){
                        weight.fixed.bic.boot[,b] <- weight.fixed.bic.boot[,b] + step.selection(factor.z,index,
                                                                                                microbes,phenotypes,type="BIC",
                                                                                                direction=direction,
                                                                                                real_data=real_data)
                    }
                }
            }
            }
            if(run.kmeans.aic.bic==TRUE){
                ## k-means partitioning
                index <- as.numeric(unlist(kmeans.rand.index[l]))
                if(length(index)!=0){
                    if(run.aic==TRUE){
                        weight.kmeans.aic.boot[,b] <- weight.kmeans.aic.boot[,b] +
                            step.selection(factor.z,index,microbes,phenotypes,
                                           type="AIC",
                                           direction=direction,
                                           real_data=real_data)
                    }

                    if(run.bic==TRUE){
                        weight.kmeans.bic.boot[,b] <- weight.kmeans.bic.boot[,b] +
                            step.selection(factor.z,index,microbes,phenotypes,
                                           type="BIC",
                                           direction=direction,
                                           real_data=real_data)
                    }
                }
            }

            if(run.kquart.aic.bic==TRUE){
                ## k-quartile partitioning
                index <- as.numeric(unlist(kquart.rand.index[l]))

                if(length(index)!=0){
                    if(run.aic==TRUE){
                        weight.kquart.aic.boot[,b] <- weight.kquart.aic.boot[,b] +
                            step.selection(factor.z,index,microbes,phenotypes,
                                           type="AIC",
                                           direction=direction,
                                           real_data=real_data)
                    }

                    if(run.bic==TRUE){
                        weight.kquart.bic.boot[,b] <- weight.kquart.bic.boot[,b] +
                            step.selection(factor.z,index,microbes,phenotypes,
                                           type="BIC",
                                           direction=direction,
                                           real_data=real_data)
                    }
                }
            }

            if(run.sort.aic.bic==TRUE){
                ## sorted partitioning
                index <- as.numeric(unlist(sort.rand.index[l]))

                if(length(index)!=0){
                    if(run.aic==TRUE){
                        weight.sort.aic.boot[,b] <- weight.sort.aic.boot[,b] +
                            step.selection(factor.z,index,microbes,phenotypes,
                                           type="AIC",
                                           direction=direction,
                                           real_data=real_data)
                    }
                    if(run.bic==TRUE){
                        weight.sort.bic.boot[,b] <- weight.sort.bic.boot[,b] +
                            step.selection(factor.z,index,microbes,phenotypes,
                                           type="BIC",
                                           direction=direction,
                                           real_data=real_data)
                    }
                }
            }
        }
    }

    if(run.aic.bic==TRUE){
        out.aic.boot <- apply(weight.aic.boot,1,sum)/nboot
        out.bic.boot <- apply(weight.bic.boot,1,sum)/nboot
    }

    #if(run.fixed.aic.bic==TRUE){
    #    out.fixed.aic.boot <- apply(weight.fixed.aic.boot,1,sum)/nboot
    #    out.fixed.bic.boot <- apply(weight.fixed.bic.boot,1,sum)/nboot
    #}

    if(run.kmeans.aic.bic==TRUE){
        out.kmeans.aic.boot <- apply(weight.kmeans.aic.boot,1,sum)/nboot
        out.kmeans.bic.boot <- apply(weight.kmeans.bic.boot,1,sum)/nboot
    }

    if(run.kquart.aic.bic==TRUE){
        out.kquart.aic.boot <- apply(weight.kquart.aic.boot,1,sum)/nboot
        out.kquart.bic.boot <- apply(weight.kquart.bic.boot,1,sum)/nboot
    }

    if(run.sort.aic.bic==TRUE){
        out.sort.aic.boot <- apply(weight.sort.aic.boot,1,sum)/nboot
        out.sort.bic.boot <- apply(weight.sort.bic.boot,1,sum)/nboot
    }

    # Lasso fitting with exclusion frequency weights

    if(run.aic.bic==TRUE){
        ## weights are exclusion frequency (random partitioning)

        #if(run.aic==TRUE){
            weights <- data.frame(out.aic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.aic.bvalue <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                                   file="weight_pval_aic_boot_",
                                                   include.diet=include.z,format.data=format.data,
                                                   diet.wt=diet.wt,
                                                   thresh.q=thresh.q,delta=delta,
                                                   std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.aic.boot <- lasso.aic.bvalue$interest
        #}

        #if(run.bic==TRUE){
            weights <- data.frame(out.bic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.bic.bvalue <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                                   file="weight_pval_bic_boot_",
                                                   include.diet=include.z,format.data=format.data,
                                                   diet.wt=diet.wt,
                                                   thresh.q=thresh.q,delta=delta,
                                                   std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.bic.boot <- lasso.bic.bvalue$interest
        #}
    }

    if(run.kmeans.aic.bic==TRUE){
        ## weights are exclusion frequency (designed partitioning-kmeans)
        #if(run.aic==TRUE){
            weights <- data.frame(out.kmeans.aic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.kmeans.aic.bvalue <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                                          file="weight_pval_kmeans_aic_boot_",
                                                          include.diet=include.z,format.data=format.data,
                                                          diet.wt=diet.wt,
                                                          thresh.q=thresh.q,delta=delta,
                                                          std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.kmeans.aic.boot <- lasso.kmeans.aic.bvalue$interest
        #}

        #if(run.bic==TRUE){
            weights <- data.frame(out.kmeans.bic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.kmeans.bic.bvalue <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                                          file="weight_pval_kmeans_bic_boot_",
                                                          include.diet=include.z,format.data=format.data,
                                                          diet.wt=diet.wt,
                                                          thresh.q=thresh.q,delta=delta,
                                                          std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.kmeans.bic.boot <- lasso.kmeans.bic.bvalue$interest
        #}
    }

    if(run.kquart.aic.bic==TRUE){
        ## weights are exclusion frequency (designed partitioning-k quartiles)
        #if(run.aic==TRUE){
            weights <- data.frame(out.kquart.aic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.kquart.aic.bvalue <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                                          file="weight_pval_kquart_aic_boot_",
                                                          include.diet=include.z,format.data=format.data,
                                                          diet.wt=diet.wt,
                                                          thresh.q=thresh.q,delta=delta,
                                                          std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.kquart.aic.boot <- lasso.kquart.aic.bvalue$interest
        #}

        #if(run.bic==TRUE){
            weights <- data.frame(out.kquart.bic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.kquart.bic.bvalue <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                                          file="weight_pval_kquart_bic_boot_",
                                                          include.diet=include.z,format.data=format.data,
                                                          diet.wt=diet.wt,
                                                          thresh.q=thresh.q,delta=delta,
                                                          std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.kquart.bic.boot <- lasso.kquart.bic.bvalue$interest
        #}
    }

    if(run.sort.aic.bic==TRUE){
        ## weights are exclusion frequency (designed partitioning-k quartiles)
        #if(run.aic==TRUE){
            weights <- data.frame(out.sort.aic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.sort.aic.bvalue <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                                        file="weight_pval_sort_aic_boot_",
                                                        include.diet=include.z,format.data=format.data,
                                                        diet.wt=diet.wt,
                                                        thresh.q=thresh.q,delta=delta,
                                                        std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.sort.aic.boot <- lasso.sort.aic.bvalue$interest
        #}

        #if(run.bic==TRUE){
            weights <- data.frame(out.sort.bic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.sort.bic.bvalue <- lasso.computations(weights,microbes,phenotypes,g1,plots=FALSE,
                                                        file="weight_pval_sort_bic_boot_",
                                                        include.diet=include.z,format.data=format.data,
                                                        diet.wt=diet.wt,
                                                        thresh.q=thresh.q,delta=delta,
                                                        std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.sort.bic.boot <- lasso.sort.bic.bvalue$interest
        #}
    }

    }

    return(list("qval"=out.qvalue,"bh.pval"=out.benhoch.pval.adjust, "pval"=out.pvalue, "out.cor"=out.cor, "out.parcor"=out.parcor, "out.benhoch.cor"=out.benhoch.cor, "out.benhoch.parcor"=out.benhoch, "out.w"=out.w, "alpha"=alpha, "alpha.bh"=alpha.bh, "delta"=delta, "cv.delta.w"=mult.delta.w5, "cv.delta.adapt"=mult.delta.w6, "cv.out.w"=mult.cv.delta.out.w5, "cv.out.adapt"=mult.cv.delta.out.w6,
                #"out.aic.boot"=out.aic.boot, "out.bic.boot"=out.bic.boot,# "out.fixed.aic.boot"=out.fixed.aic.boot, "out.fixed.bic.boot"=out.fixed.bic.boot,
                #"out.kmeans.aic.boot"=out.kmeans.aic.boot, "out.kmeans.bic.boot"=out.kmeans.bic.boot, "out.kquart.aic.boot"=out.kquart.aic.boot, "out.kquart.bic.boot"=out.kquart.bic.boot,
                #"out.sort.aic.boot"=out.sort.aic.boot, "out.sort.bic.boot"=out.sort.bic.boot,
                "w.aic.boot"=w.aic.boot, "w.bic.boot"=w.bic.boot, "w.kmeans.aic.boot"=w.kmeans.aic.boot, "w.kmeans.bic.boot"=w.kmeans.bic.boot,
                "w.kquart.aic.boot"=w.kquart.aic.boot, "w.kquart.bic.boot"=w.kquart.bic.boot, "w.sort.aic.boot"=w.sort.aic.boot, "w.sort.bic.boot"=w.sort.bic.boot))
}
