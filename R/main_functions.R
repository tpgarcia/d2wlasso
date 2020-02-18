####################################
####################################
##
##
## Functions for the main method
##
####################################
####################################

#' Implement structured variables selection with q-values
#'
#' @param x (n by m) matrix of main covariates where m is the number of covariates and n is the sample size.
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable. This covariate should
#' always be selected. Can be NULL.
#' @param y (n by K) a matrix corresponding to K response variables. If \code{regression.type} is "cox",
#' \code{y} contains the observed event times for each of the k response types, k=1,...,K.
#' @param cox.delta (n by K) a matrix that denotes censoring when \code{regression.type} is "cox" (1 denotes
#' survival event is observed, 0 denotes the survival event is censored). Can be NULL.
#' @param factor.z logical. If TRUE, the fixed variable z is a factor variable.
#' @param regression.type a character indicator that is either "linear" for linear regression
#' or "cox" for Cox proportional hazards regression. Default is "linear".
#' @param ttest logical. If TRUE, p-value for each covariate is computed from univariate
#' linear/cox regression of the response on each covariate. If FALSE, the
#' p-value is computed from correlation coefficients. Default is FALSE.
#' @param q_method indicates the method for choosing optimal tuning parameter
#' in the q-value computation as proposed in Storey and Tibshirani (2003).
#' One of "bootstrap" or "smoother". Default is "smoother" (smoothing spline).
#' @param plots logical. If TRUE, figures are plotted. Default is FALSE.
#' @param pi0.true logical. If TRUE, the estimate of the true proportion of the null hypothesis is set to the value of pi0.val which is given by the user. If FALSE, the estimate of the true proportion of the null hypothesis is computed by bootstrap or smoothing spline. Default is FALSE.
#' @param pi0.val A user supplied estimate of the true proportion of the null hypothesis. Used only when pi0.true is TRUE. Default is 0.9.
#' @param wt The weights to be used for the weighted lasso. One of "one","t_val","parcor","p_val","bhp_val","adapt","q_cor" or "q_parcor".
#' "one" gives no weight. "t_val" gives weight of the inverse absolute t-statistics of the regression coefficients. "parcor" gives weight of the inverse absolute partial correlation between the main covariate and the response after accounting for z. "p_val" gives p-value of each predictor's coefficient as weights. "bhp_val" gives Benjamini-Hochberg adjusted p-value of each predictor's coefficient as weights. "adapt" gives adaptive lasso weights, that is, the inverse of the absolute value of regression coefficients. "q_cor" gives weights set to q-values BEFORE taking into account diet. "q_parcor" gives weights set to q-values AFTER taking into account diet.
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
#' dwlcox1 <- d2wlasso(x,z,y,cox.delta = cox.delta, regression.type = "cox")
#' dwlcox2 <- d2wlasso(x,z,y,cox.delta = cox.delta, regression.type = "cox", nboot = 50)
#' dwlcox3 <- d2wlasso(x,z,y,cox.delta = cox.delta, regression.type = "cox", wt="t_val")
#' dwlcoxcv1 <- d2wlasso(x,z,y,cox.delta = cox.delta,regression.type = "cox",lasso.delta.cv.mult = TRUE, ncv = 3, nboot = 50)
d2wlasso <- function(x,z,y,
                     cox.delta=NULL,
                     factor.z=TRUE,
                     regression.type=c("linear","cox")[1],

                     ttest=TRUE,
                     q_method=c("bootstrap","smoother")[2],
                     plots=FALSE,
                     pi0.true=FALSE,
                     pi0.val=0.9,
                     wt=c("one","t_val","parcor","p_val","bhp_val","adapt","q_cor","q_parcor")[7],
                     weight_fn=c("identity","sqrt","inverse_abs","square")[1],
                     include.z=TRUE,
                     z.wt=1000,
                     thresh.q=TRUE,
                     alpha=0.15,
                     alpha.bh=0.05,
                     delta=2,
                     robust=TRUE,
                     q.old=FALSE,
                     lasso.delta.cv.mult=FALSE,
                     vfold=10,
                     ncv=100,
                     delta.cv.seed=NULL,
                     run.aic.bic=TRUE,
                     run.kmeans.aic.bic=TRUE,
                     run.kquart.aic.bic=TRUE,
                     run.sort.aic.bic=TRUE,
                     nboot=100,
                     k=4,
                     direction="backward"){

    cv.criterion=FALSE

    #################################
    # Store the dimensions of x, z ##
    #################################
    n <- nrow(x)
    m <- ncol(x)

    if(!is.null(z)){
        m0 <- ncol(z)
    } else {
        m0 <- 0
    }

    #########################
    # arranging input data ##
    #########################
    if(!is.null(z)){
        XX <- cbind(z, x)
    } else {
        XX <- x
    }

    ## allocate names to X
    colnames_use <- NULL
    if(!is.null(z)){
        if(is.null(colnames(z))){
            colnames_use <- c(colnames_use,"Fixed")
            colnames(z) <- "Fixed"
        } else {
            colnames_use <- c(colnames_use,colnames(z))
        }
    }

    if(is.null(colnames(x))){
        colnames_use <- c(colnames_use,paste0("X_",seq(1,m)))
        colnames(x) <- paste0("X_",seq(1,m))
    } else {
        colnames_use <- c(colnames_use,colnames(x))
    }

    colnames(XX) <- colnames_use

    if (regression.type=="linear"){
        response <- list(yy=y,delta=NULL)
    } else if (length(cox.delta)!=length(y)){
        stop("length of y should be equal to length of cox.delta")
    } else {
        response <- list(yy=y,delta=cox.delta)
    }

    ########################################################
    ## Determine number of covariates and covariate names ##
    ########################################################
    number.of.covariates <- m0+m
    covariate.names <- colnames(XX)


    ############################################################
    ##                                                        ##
    ##                                                        ##
    ##          COMPUTE THE WEIGHTS                           ##
    ##                                                        ##
    ##                                                        ##
    ############################################################
    if(weight.type=="one"){
        ## Unit weights
        weights <- matrix(1,nrow=m,ncol=ncol(response$yy))
    } else if(weight.type=="corr.estimate" |
       weight.type=="corr.pvalue" |
       weight.type=="corr.bh.pvalue" |
       weight.type=="corr.tstat"|
       weight.type=="corr.qvalue"){

        ####################################
        ## compute correlation of x and y ##
        ####################################
        cor.out <- correlations(factor.z,x,z,response,partial=FALSE,ttest=ttest,
                                regression.type=regression.type)

        #########################
        ## Extract the weights ##
        #########################
        if(weight.type=="corr.tstat"){
            ## t-values with \beta_k in y= \beta_k * x_k
            weights <- cor.out$tvalues
            threshold.selection <- NULL
        } else if(weight.type=="corr.pvalue"){
            ## p-values with \beta_k in y= \beta_k * x_k
            weights <- cor.out$pvalues

        } else if(weight.type=="corr.estimate"){
            ## partial correlations
            weights <- cor.out$estimate
            threshold.selection <- NULL
        } else if(weight.type=="corr.bh.pvalue"){
            ## Benjamini-Hochberg adjusted p-values from y=z + \beta_k * x_k
            benhoch.results <- ben.hoch.interest(cor.out$pvalues,alpha=alpha.bh)
            weights <- benhoch.results$pval.adjust
        } else if(weight.type=="corr.qvalue"){

            ## Results for testing if a x_k has an effect on y, but NOT
            ##            accounting for z
            ## That is, we test : H_0 : \beta_{x_k}=0
            qvalues.results <- q.computations(cor.out, method=q_method,
                                              plots=plots,file="cor",robust=robust,q.old=q.old,
                                              pi0.true=pi0.true,pi0.val=pi0.val)

            weights <- qvalues.results$qval.mat
            threshold.selection <- q.interest(weights,alpha=alpha,criteria="less")
            if(!is.null(z)){
                ## we ignore z
                threshold.selection <- c(0,t(threshold.selection$interest))
            }

        }


    } else if(weight.type=="parcor.estimate" |
              weight.type=="parcor.pvalue" |
              weight.type=="parcor.bh.pvalue" |
              weight.type=="parcor.tstat"|
              weight.type=="parcor.qvalue"){
        ##################################################################
        ## compute partial correlation of x and y after adjusting for z ##
        ##################################################################

        if(!is.null(z)){
            parcor.out <- correlations(factor.z,x,z,response,partial=TRUE,ttest=ttest,
                                       regression.type=regression.type)

            #########################
            ## Extract the weights ##
            #########################
            if(weight.type=="parcor.tstat"){
                ## t-values with \beta_k in y=z + \beta_k * x_k
                weights <- parcor.out$tvalues
            } else if(weight.type=="parcor.pvalue"){
                ## p-values with \beta_k in y=z + \beta_k * x_k
                weights <- parcor.out$pvalues
            } else if(weight.type=="parcor.estimate"){
                ## partial correlations
                weights <- parcor.out$estimate
            } else if(weight.type=="parcor.bh.pvalue"){
                ## Benjamini-Hochberg adjusted p-values from y=z + \beta_k * x_k
                benhoch.results <- ben.hoch.interest(parcor.out$pvalues,alpha=alpha.bh)
                weights <- benhoch.results$pval.adjust
                threshold.selection <- c(1,t(benhoch.results$interest))

            } else if(weight.type=="parcor.qvalue"){
                ## Results for testing if a microbe has an effect on phenotype, but AFTER
                ##            accounting for diet
                ## That is, we test : H_0 : \beta_{x_j|z}=0
                # compute q-value as used by JD Storey with some adjustments made
                qvalues.results <- q.computations(parcor.out,method=q_method,
                                                             plots=plots,file="parcor",robust=robust,q.old=q.old,
                                                             pi0.true=pi0.true,pi0.val=pi0.val)
                weights <- qvalues.results$qval.mat
                threshold.selection <- q.interest(microbe.parcor.out.qvalues$qval.mat,alpha=alpha,criteria="less")
                threshold.selection <- c(1,t(threshold.selection$interest))

            }
        } else {
            warning("Your choice of weight.type requires z to be non-NULL.")
        }
    } else if(weight.type=="exfrequency.kmeans" |
              weight.type=="exfrequency.kquartiles"|
              weight.type=="exfrequency.ksorted"){

        k.split=k
        run.aic = TRUE
        run.bic = TRUE


        ## Store weights for each bootstrap
        if(!is.null(z)){
            tmp2.store <- as.data.frame(matrix(0, nrow = (number.of.covariates-1), ncol = nboot,
                                           dimnames = list(covariate.names[-1],
                                                           seq(1,nboot))))
        } else {
            tmp2.store <- as.data.frame(matrix(0, nrow = number.of.covariates, ncol = nboot,
                                               dimnames = list(covariate.names,
                                                               seq(1,nboot))))
        }
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

        run.kmeans.aic.bic <- FALSE
        run.kquart.aic.bic <- FALSE
        run.sort.aic.bic <- FALSE

        if(weight.type=="exfrequency.kmeans"){
            run.kmeans.aic.bic <- TRUE
        }

        if(weight.type=="exfrequency.kquartiles"){
            run.kquart.aic.bic <- TRUE
        }

        if(weight.type=="exfrequency.ksorted"){
            run.sort.aic.bic <- TRUE
        }

        #####################################
        ## Compute intermediary quantities ##
        #####################################
        if(run.kmeans.aic.bic==TRUE | run.kquart.aic.bic==TRUE |  run.sort.aic.bic==TRUE){

            ## Get parameter estimates from ridge regression
            beta.values <- ridge.regression(XX,response)

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
                index.group.sort <- index.sort.partition(n=nrow(XX),k=k,sort.beta.index)
            }
        }

        for(b in 1:nboot){
            ##print(b)
            if(run.aic.bic==TRUE){
                ## Randomly partition the index
                rand.index <- random.partition(n=nrow(XX),p=ncol(XX)-1,k=k)
            }

            #if(run.fixed.aic.bic==TRUE){
            ## Ensure fixed covariates are in the partition + randomly partition the rest
            #    rand.fixed.index <- fixed.plus.random.partition(fixed.covariates,n=nrow(XX),p=ncol(XX)-1,k=k)
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
                                step.selection(factor.z,index,XX,response,type="AIC",
                                               direction=direction,
                                               real_data=real_data)
                        }

                        if(run.bic==TRUE){
                            weight.bic.boot[,b] <- weight.bic.boot[,b] +
                                step.selection(factor.z,index,XX,response,type="BIC",
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
                                                   XX,response,type="AIC",
                                                   direction=direction,
                                                   real_data=real_data)
                            }

                            if(run.bic==TRUE){
                                weight.fixed.bic.boot[,b] <- weight.fixed.bic.boot[,b] + step.selection(factor.z,index,
                                                                                                        XX,response,type="BIC",
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
                                step.selection(factor.z,index,XX,response,
                                               type="AIC",
                                               direction=direction,
                                               real_data=real_data)
                        }

                        if(run.bic==TRUE){
                            weight.kmeans.bic.boot[,b] <- weight.kmeans.bic.boot[,b] +
                                step.selection(factor.z,index,XX,response,
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
                                step.selection(factor.z,index,XX,response,
                                               type="AIC",
                                               direction=direction,
                                               real_data=real_data)
                        }

                        if(run.bic==TRUE){
                            weight.kquart.bic.boot[,b] <- weight.kquart.bic.boot[,b] +
                                step.selection(factor.z,index,XX,response,
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
                                step.selection(factor.z,index,XX,response,
                                               type="AIC",
                                               direction=direction,
                                               real_data=real_data)
                        }
                        if(run.bic==TRUE){
                            weight.sort.bic.boot[,b] <- weight.sort.bic.boot[,b] +
                                step.selection(factor.z,index,XX,response,
                                               type="BIC",
                                               direction=direction,
                                               real_data=real_data)
                        }
                    }
                }
            }
        }

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

    ############################################################
    ##                                                        ##
    ##                                                        ##
    ##          COMPUTE THE WEIGHTED LASSO RESULTS            ##
    ##                                                        ##
    ##                                                        ##
    ############################################################
    if(lasso.delta.cv.mult==FALSE){
        lasso.w <- lasso.computations(weights,XX,response,g,plots=plots,file="weight_",
                                      include.diet=include.z,diet.wt=z.wt,thresh.q=thresh.q,
                                      delta=delta)
        weighted.lasso <- as.matrix(lasso.w$interest)

    } else if(lasso.delta.cv.mult==TRUE){

        lasso.w <-0

        for(v in 1:ncv){
            lasso.w.tmp <- lasso.computations(weights,XX,response,g3,plots=FALSE,file="weight6_",
                                                         include.diet=include.diet,diet.wt=z.wt,corr.g=TRUE,delta=delta,
                                                         cv.criterion=FALSE,vfold=vfold)
            lasso.w <- lasso.w + as.matrix(lasso.w.tmp$interest)
            ##mult.delta.w6[,v] <- mult.delta.w6[,v] + mult.cv.delta.lasso.w6$delta.out
        }

        weighted.lasso <- lasso.w
    }


    ## exclusion frequency weights

    if (regression.type=="cox"){





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
            lasso.aic.bvalue <- lasso.computations(weights,XX,response,g1,plots=FALSE,
                                                   file="weight_pval_aic_boot_",
                                                   include.diet=include.z,
                                                   diet.wt=diet.wt,
                                                   thresh.q=thresh.q,delta=delta,
                                                   std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.aic.boot <- lasso.aic.bvalue$interest
            #}

            #if(run.bic==TRUE){
            weights <- data.frame(out.bic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.bic.bvalue <- lasso.computations(weights,XX,response,g1,plots=FALSE,
                                                   file="weight_pval_bic_boot_",
                                                   include.diet=include.z,
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
            lasso.kmeans.aic.bvalue <- lasso.computations(weights,XX,response,g1,plots=FALSE,
                                                          file="weight_pval_kmeans_aic_boot_",
                                                          include.diet=include.z,
                                                          diet.wt=diet.wt,
                                                          thresh.q=thresh.q,delta=delta,
                                                          std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.kmeans.aic.boot <- lasso.kmeans.aic.bvalue$interest
            #}

            #if(run.bic==TRUE){
            weights <- data.frame(out.kmeans.bic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.kmeans.bic.bvalue <- lasso.computations(weights,XX,response,g1,plots=FALSE,
                                                          file="weight_pval_kmeans_bic_boot_",
                                                          include.diet=include.z,
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
            lasso.kquart.aic.bvalue <- lasso.computations(weights,XX,response,g1,plots=FALSE,
                                                          file="weight_pval_kquart_aic_boot_",
                                                          include.diet=include.z,
                                                          diet.wt=diet.wt,
                                                          thresh.q=thresh.q,delta=delta,
                                                          std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.kquart.aic.boot <- lasso.kquart.aic.bvalue$interest
            #}

            #if(run.bic==TRUE){
            weights <- data.frame(out.kquart.bic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.kquart.bic.bvalue <- lasso.computations(weights,XX,response,g1,plots=FALSE,
                                                          file="weight_pval_kquart_bic_boot_",
                                                          include.diet=include.z,
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
            lasso.sort.aic.bvalue <- lasso.computations(weights,XX,response,g1,plots=FALSE,
                                                        file="weight_pval_sort_aic_boot_",
                                                        include.diet=include.z,
                                                        diet.wt=diet.wt,
                                                        thresh.q=thresh.q,delta=delta,
                                                        std.y=std.y,est.MSE=est.MSE,cv.criterion=cv.criterion)
            w.sort.aic.boot <- lasso.sort.aic.bvalue$interest
            #}

            #if(run.bic==TRUE){
            weights <- data.frame(out.sort.bic.boot)
            colnames(weights) <- "response"
            rownames(weights) <- out.rownames[-1]
            lasso.sort.bic.bvalue <- lasso.computations(weights,XX,response,g1,plots=FALSE,
                                                        file="weight_pval_sort_bic_boot_",
                                                        include.diet=include.z,
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



###############
## Libraries ##
###############

library(xtable)		# to create LaTeX tables
library(lars)		# for LASSO approach
library(plotrix)		# for computing standard errors of mean in simulation study
library(MASS) # for ridge regression
library(glmnet)     # for ridge regression


#########################################################
## Function for making data frames for storing results ##
#########################################################

## EXPORT
store.xy <- function(XX,response){
	out <- as.data.frame(matrix(0,nrow=ncol(XX),ncol=ncol(response$yy)))
	rownames(out) <- colnames(XX)
	colnames(out) <- colnames(response$yy)
	return(out)
}

## EXPORT
store.x <- function(XX){
	out <- as.data.frame(rep(0,ncol(XX)))
	rownames(out) <- colnames(XX)
	return(out)
}


##########################################################
## Functions to get partial correlation and its p-value #
##########################################################

# When pearson correlation is 0 (because std. deviation is 0), I am setting
# p-value  to 1.

## EXPORT
#' @importFrom stats lm
#' @importFrom survival coxph Surv
corr.pvalue <- function(x,y,delta,method="pearson",alternative="two.sided",ttest=FALSE,regression.type){
    x <- as.numeric(x)
    y <- as.numeric(y)
    delta.use <- as.numeric(delta)

    if(regression.type=="linear"){
        out <- stats::cor.test(x,y,alternative=alternative,method=method,na.action=na.omit)
        estimate <- out$estimate
    } else {
        estimate <- 0
    }

    if(ttest==FALSE & regression.type=="linear"){
        p.value <- out$p.value
        t.stat=NULL
    } else {
        y1 <- y
        x1 <- x
        delta1 <- delta.use

        if(regression.type=="linear"){
            #######################
            ## linear regression ##
            #######################
            summary.out <- summary(stats::lm(y1~  x1))
            p.value <- summary.out$coefficients["x1","Pr(>|t|)"]
            t.stat <- summary.out$coefficients["x1","t value"]
            ##summary.out <- summary(lm(y~x))
            ##p.value <- summary.out$coefficients["x","Pr(>|t|)"]
        } else {
            #########################
            ## survival regression ##
            #########################
            survival.data <- list(time=y1,status=delta1,x1=x1)
            summary.out <- summary(survival::coxph(survival::Surv(time,status)~ x1,data=survival.data))
            p.value <- summary.out$coefficients["x1","Pr(>|z|)"]
            t.stat  <- summary.out$coefficients["x1","z"]
        }

    }
    list(p.value=p.value,estimate=estimate,t.stat=t.stat)
}


## EXPORT
#' @importFrom stats lm
#' @importFrom survival coxph Surv
parcorr.pvalue <- function(factor.z,x,y,z,delta=NULL,method="pearson",alternative="two.sided",ttest=FALSE,regression.type){
    x <- as.numeric(x)
    y <- as.numeric(y)
    z <- as.numeric(z)
    delta.use <- as.numeric(delta)

    if(regression.type=="linear"){
        if(factor.z==TRUE){
            xres <- residuals(stats::lm(x~factor(z)))
            yres <- residuals(stats::lm(y~factor(z)))
        } else {
            xres <- residuals(stats::lm(x~z))
            yres <- residuals(stats::lm(y~z))
        }
        out <- corr.pvalue(xres,yres,delta.use,method,alternative,ttest=FALSE,regression.type=regression.type)
        estimate <- out$estimate
    } else {
        estimate <- 0
    }

    if(ttest==FALSE & regression.type=="linear"){
        p.value <- out$p.value
        t.stat <- NULL
    } else {
        y1 <- y
        x1 <- x
        delta1 <- delta.use

        if(regression.type=="linear"){
            #######################
            ## linear regression ##
            #######################

            if(factor.z==TRUE){
                summary.out <- summary(stats::lm(y1 ~ factor(z) +  x1))  ## may have to change due to nature of z!!
            } else {
                summary.out <- summary(stats::lm(y1 ~ z +  x1))  ## may have to change due to nature of z!!
            }

            p.value <- summary.out$coefficients["x1","Pr(>|t|)"]
            t.stat <- summary.out$coefficients["x1","t value"]

            ## Or using reduced sum of squares:
            ##model.full <- lm(y1~z + x1)
            ##model.red  <- lm(y1~z)
            ##anova.out <- anova(model.full,model.red)
            ##p.value <- anova.out[-1,"Pr(>F)"]

            ##summary.out <- summary(lm(y ~ z +  x))
            ##p.value <- summary.out$coefficients["x","Pr(>|t|)"]
        } else {
            #########################
            ## survival regression ##
            #########################
            survival.data <- list(time=y1,status=delta1,z=z,x1=x1)
            if(factor.z==TRUE){
                summary.out <- summary(survival::coxph(survival::Surv(time,status)~ factor(z) + x1,
                                                    data=survival.data))
            } else{
                summary.out <- summary(survival::coxph(survival::Surv(time,status)~ z + x1,
                                                       data=survival.data))
            }
            p.value <- summary.out$coefficients["x1","Pr(>|z|)"]
            t.stat  <- summary.out$coefficients["x1","z"]
        }
    }


    list(p.value=p.value,estimate=estimate,t.stat=t.stat)
}


## EXPORT
correlations <- function(factor.z,x,z,response,partial=FALSE,ttest=FALSE,regression.type){

    ## Formatting data
    data.response <- response$yy
    data.delta <- response$delta

	data.XX <- x
	data.z <- z

	# Setting up matrices to store Pearson correlations and p-values
	correlation <- store.xy(data.XX,data.response)
	pvalues <- store.xy(data.XX,data.response)
	tvalues <- store.xy(data.XX,data.response)

	# Computing pearson correlations and p-values
	for (i in 1: ncol(data.XX)) {
		for(j in 1:ncol(data.response)) {
			if(partial==TRUE){
				tmp <- parcorr.pvalue(factor.z=factor.z,x=data.XX[,i],y=data.response[,j],
				                      z=data.z,delta=data.delta[j,],ttest=ttest,regression.type=regression.type)
			} else {
				tmp <- corr.pvalue(x=data.XX[,i],y=data.response[,j],delta=data.delta[,j],
				                   ttest=ttest,regression.type=regression.type)
			}
			correlation[i,j] <- tmp$estimate
			pvalues[i,j] <- tmp$p.value
			tvalues[i,j] <- tmp$t.stat
		}
	}
	list(estimate = correlation, pvalues = pvalues,tvalues=tvalues)
}



# function to compute partial correlations using pcor.R
# NOT USED.
correlations.pcor <- function(XX,response,partial="FALSE",ttest="FALSE"){
    ## Formatting data
    data.response <- response$yy

    data.XX <- XX[,-which(rownames(XX)=="Fixed")]
    diet <- XX[,"Fixed"]

    ## Setting up matrices to store Pearson correlations and p-values
    correlation <- store.xy(data.XX,data.response)
    pvalues <- store.xy(data.XX,data.response)

    ## Computing pearson correlations and p-values
    for (i in 1: ncol(data.XX)) {
        for(j in 1:ncol(data.response)) {
            tmp <- pcor.test(as.numeric(data.XX[i,]),
                             as.numeric(data.response[,j]),as.numeric(diet))

            correlation[i,j] <- tmp$estimate
            pvalues[i,j] <- tmp$p.value
        }
    }
    list(estimate = correlation, pvalues = pvalues)
}



##################################################
# Function to get F-test statistics and p-values #
##################################################
# EXPORT
fstat.pvalue <- function(factor.z,x,z){
    if(factor.z==TRUE){
        z <- factor(as.numeric(z))
    } else {
        z <- as.numeric(z)
    }
    x <- as.numeric(x)

    # Assuming variances in both diet groups are different
    result <- t.test(x ~ z)
    stat <- result$statistic
    p.value <- result$p.value

    #if(is.na(Fstat)){
    #	Fstat <- NA
    #	p.value <- 1
    #}
    return(list(Fstat = stat, p.value = p.value))
}

# EXPORT
ftests <- function(factor.z,x,z){
    # Formatting data
    data.z <- z
    data.XX <- x

    # Setting up matrices to store F-statistics and p-values
    fstat <- store.x(data.XX)
    pvalues <- store.x(data.XX)

    # Computing pearson correlations and p-values
    for (i in 1: ncol(data.XX)) {
        tmp <- fstat.pvalue(factor.z,data.XX[,i],z)
        fstat[i,] <- tmp$Fstat
        pvalues[i,] <- tmp$p.value
    }

    return(list(estimate = fstat, pvalues = pvalues))
}


####################
# Qvalues function #
####################

# modification of qplot so as to get individual plots

qplot2 <- function (qobj, rng = c(0, 0.1), smooth.df = 3, smooth.log.pi0 = FALSE, whichplot=1)
{
    q2 <- qobj$qval[order(qobj$pval)]
    if (min(q2) > rng[2]) {
        rng <- c(min(q2), quantile(q2, 0.1))
    }
    p2 <- qobj$pval[order(qobj$pval)]
    #par(mfrow = c(2, 2))
    lambda <- qobj$lambda
    if (length(lambda) == 1) {
        lambda <- seq(0, max(0.9, lambda), 0.05)
    }
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
        pi0[i] <- mean(p2 > lambda[i])/(1 - lambda[i])
    }
    if (smooth.log.pi0)
        pi0 <- log(pi0)
    spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
    if (smooth.log.pi0) {
        pi0 <- exp(pi0)
        spi0$y <- exp(spi0$y)
    }
    pi00 <- round(qobj$pi0, 3)
	if(whichplot==1){
	#TG change, 02/07/2012 ("gamma" is lambda in manuscript)
	par(mar=c(5, 4, 4, 2)+1)
    plot(lambda, pi0, xlab = expression(gamma), ylab = expression(hat(pi)(gamma)),pch=19,cex.lab=1.5,cex.axis=1.5)
    mtext(substitute(hat(pi) == that, list(that = pi00)),cex=2)
    lines(spi0,lwd=2)
	} else if (whichplot==2){
    plot(p2[q2 >= rng[1] & q2 <= rng[2]], q2[q2 >= rng[1] & q2 <=
        rng[2]], type = "l", xlab = "p-value", ylab = "q-value")
	} else if (whichplot==3){
    plot(q2[q2 >= rng[1] & q2 <= rng[2]], (1 + sum(q2 < rng[1])):sum(q2 <=
        rng[2]), type = "l", xlab = "q-value cut-off", ylab = "Number of significant XX",cex.lab=1.5,cex.axis=1.5)
	} else if (whichplot==4){
    plot((1 + sum(q2 < rng[1])):sum(q2 <= rng[2]), q2[q2 >= rng[1] &
        q2 <= rng[2]] * (1 + sum(q2 < rng[1])):sum(q2 <= rng[2]),
        type = "l", xlab = "significant tests", ylab = "expected false positives")
	}
    #par(mfrow = c(1, 1))
}

####################################################################
## qvalue.adj. is as used by JD Storey with some adjustments made ##
####################################################################

qvalue.adj<-function (p = NULL, lambda = seq(0, 0.9, 0.05), pi0.method = "smoother",
    fdr.level = NULL, robust = FALSE, gui = FALSE, smooth.df = 3,
    smooth.log.pi0 = FALSE,pi0.true=FALSE,pi0.val=0.9)
{
    if (is.null(p)) {
        qvalue.gui()
        return("Launching point-and-click...")
    }
    if (gui & !interactive())
        gui = FALSE
    if (min(p) < 0 || max(p) > 1) {
        if (gui)
            eval(expression(postMsg(paste("ERROR: p-values not in valid range.",
                "\n"))), parent.frame())
        else print("ERROR: p-values not in valid range.")
        return(0)
    }
    if (length(lambda) > 1 && length(lambda) < 4) {
        if (gui)
            eval(expression(postMsg(paste("ERROR: If length of lambda greater than 1, you need at least 4 values.",
                "\n"))), parent.frame())
        else print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
        return(0)
    }
    if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >=
        1)) {
        if (gui)
            eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).",
                "\n"))), parent.frame())
        else print("ERROR: Lambda must be within [0, 1).")
        return(0)
    }
    m <- length(p)
    if (length(lambda) == 1) {
        if (lambda < 0 || lambda >= 1) {
            if (gui)
                eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).",
                  "\n"))), parent.frame())
            else print("ERROR: Lambda must be within [0, 1).")
            return(0)
        }
        pi0 <- mean(p >= lambda)/(1 - lambda)
        pi0 <- min(pi0, 1)
    }
    else {
      # TG added pi0.true option
      if(pi0.true==TRUE){
        pi0 <- pi0.val
      } else {
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
          pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }

     # TG change: Remove any pi0 which have entry 0, and adjust lambda values

        if(sum(pi0==0)>0){
          ind.zero <- which(pi0==0)
          pi0 <- pi0[-ind.zero]
          lambda <- lambda[-ind.zero]
        }

        if (pi0.method == "smoother") {
          if (smooth.log.pi0)
            pi0 <- log(pi0)
          spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
          pi0 <- predict(spi0, x = max(lambda))$y
          if (smooth.log.pi0)
            pi0 <- exp(pi0)
          pi0 <- min(pi0, 1)
        }
        else if (pi0.method == "bootstrap") {
          minpi0 <- min(pi0)
          mse <- rep(0, length(lambda))
          pi0.boot <- rep(0, length(lambda))
          for (i in 1:1000) {
            p.boot <- sample(p, size = m, replace = TRUE)
            for (j in 1:length(lambda)) {
              pi0.boot[j] <- mean(p.boot > lambda[j])/(1 -
                                                       lambda[j])
            }
            mse <- mse + (pi0.boot - minpi0)^2
          }
          pi0 <- min(pi0[mse == min(mse)])
          pi0 <- min(pi0, 1)
        }
        else {
          print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
          return(0)
        }
      }
    }

    if (pi0 <= 0) {
        if (gui)
            eval(expression(postMsg(paste("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.",
                "\n"))), parent.frame())
        else print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
        return(0)
    }
    if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level >
        1)) {
        if (gui)
            eval(expression(postMsg(paste("ERROR: 'fdr.level' must be within (0, 1].",
                "\n"))), parent.frame())
        else print("ERROR: 'fdr.level' must be within (0, 1].")
        return(0)
    }
    u <- order(p)
    qvalue.rank <- function(x) {
        idx <- sort.list(x)
        fc <- factor(x)
        nl <- length(levels(fc))
        bin <- as.integer(fc)
        tbl <- tabulate(bin)
        cs <- cumsum(tbl)
        tbl <- rep(cs, tbl)
        tbl[idx] <- tbl
        return(tbl)
    }
    v <- qvalue.rank(p)
    qvalue <- pi0 * m * p/v
    if (robust) {
        qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for (i in (m - 1):1) {
        qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
    }
    if (!is.null(fdr.level)) {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue,
            pvalues = p, fdr.level = fdr.level, significant = (qvalue <=
                fdr.level), lambda = lambda)
    }
    else {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue,
            pvalues = p, lambda = lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
}


qvalue.old <- function(p, alpha=NULL, lam=NULL, robust=F,pi0.true=FALSE,pi0.val=0.9)
{
    #This is a function for estimating the q-values for a given set of p-values. The
    #methodology comes from a series of recent papers on false discovery rates by John
    #D. Storey et al. See http://www.stat.berkeley.edu/~storey/ for references to these
    #papers. This function was written by John D. Storey. Copyright 2002 by John D. Storey.
    #All rights are reserved and no responsibility is assumed for mistakes in or caused by
    #the program.
    #
    #Input
    #=============================================================================
    #p: a vector of p-values (only necessary input)
    #alpha: a level at which to control the FDR (optional)
    #lam: the value of the tuning parameter to estimate pi0 (optional)
    #robust: an indicator of whether it is desired to make the estimate more robust
    #        for small p-values (optional)
    #
    #Output
    #=============================================================================
    #remarks: tells the user what options were used, and gives any relevant warnings
    #pi0: an estimate of the proportion of null p-values
    #qvalues: a vector of the estimated q-values (the main quantity of interest)
    #pvalues: a vector of the original p-values
    #significant: if alpha is specified, and indicator of whether the q-value fell below alpha
    #    (taking all such q-values to be significant controls FDR at level alpha)

    #This is just some pre-processing
    m <- length(p)
    #These next few functions are the various ways to automatically choose lam
    #and estimate pi0
    if(!is.null(lam)) {
        pi0 <- mean(p>lam)/(1-lam)
        remark <- "The user prespecified lam in the calculation of pi0."
    }
    else{
        remark <- "A smoothing method was used in the calculation of pi0."
        #library(modreg)
        lam <- seq(0,0.95,0.01)
        pi0 <- rep(0,length(lam))
        for(i in 1:length(lam)) {
            pi0[i] <- mean(p>lam[i])/(1-lam[i])
        }
        spi0 <- smooth.spline(lam,pi0,df=3,w=(1-lam))
        pi0 <- predict(spi0,x=0.95)$y
    }

    # TG added, true pi0
    if(pi0.true==TRUE){
        pi0 <- pi0.val
    }

    #print(pi0)

    #The q-values are actually calculated here
    u <- order(p)
    v <- rank(p)
    qvalue <- pi0*m*p/v
    if(robust) {
        qvalue <- pi0*m*p/(v*(1-(1-p)^m))
        remark <- c(remark, "The robust version of the q-value was calculated. See Storey JD (2002) JRSS-B 64: 479-498.")
    }
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {
        qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
    }
    #Here the results are returned
    if(!is.null(alpha)) {
        list(remarks=remark, pi0=pi0, qvalue=qvalue, significant=(qvalue <= alpha), pvalue=p)
    }
    else {
        list(remarks=remark, pi0=pi0, qvalue=qvalue, pvalue=p)
    }
}


q.computations <- function(out, method=c("smoother","bootstrap")[2],
				plots=TRUE,file="name",robust=TRUE,q.old=FALSE,
				pi0.true=FALSE,pi0.val=0.9){

	qval.mat <- matrix(0,nrow=nrow(out$pvalues),ncol=ncol(out$pvalues))
	qval.mat <- as.data.frame(qval.mat)
	rownames(qval.mat) <- rownames(out$pvalues)
	colnames(qval.mat) <- colnames(out$pvalues)


	for(i in 1:ncol(qval.mat)){
		pvalues <- out$pvalues[,i]
		estimate <- out$estimate[,i]

		#qobj <- qvalue.adj(pvalues,pi0.method=method,lambda=seq(0,0.95,by=0.01),robust=FALSE,pi0.true=pi0.true,pi0.val=pi0.val)
		#qval <- qobj$qvalues
		if(q.old=="FALSE"){
		    qobj <- qvalue.adj(pvalues,pi0.method=method,lambda=seq(0,0.95,by=0.01),
		                       robust=robust,pi0.true=pi0.true,pi0.val=pi0.val)
		    qval <- qobj$qvalues
		} else {
		    qobj <- qvalue.old(pvalues,robust=robust,pi0.true=pi0.true,pi0.val=pi0.val)
		    qval <- qobj$qvalue
		}
		pi0 <- qobj$pi0
		qval.mat[,i] <- qval

		cnames <- colnames(out$pvalues)

		# Plots
		if(plots==TRUE){

			# Density histogram of p-values
			postscript(paste(file,"_histpval_",cnames[i],".eps",sep=""))
			hist(pvalues,main="",xlab="Density of observed p-values",ylab="",freq=FALSE,yaxt="n",ylim=c(0,5),cex.lab=2,cex.axis=2)
			abline(h=1,lty=1)
			abline(h=qobj$pi0,lty=2)
			axis(2,at = c(0,qobj$pi0,1,2,3,4),labels=c(0,round(qobj$pi0,2),1,2,3,4),las=2,cex.lab=2,cex.axis=2)
			box(lty=1,col="black") 	# for a box around plot
			dev.off()

			# Density histogram of estimate
			postscript(paste(file,"_histest_",cnames[i],".eps",sep=""))
			hist(estimate,main="",xlab="Density of observed statistic",ylab="",yaxt="n",freq=FALSE,ylim=c(0,5),cex.lab=2,cex.axis=2)
			axis(2,at = c(0,1,2,3,4),labels=c(0,1,2,3,4),las=2,cex.lab=2,cex.axis=2)
			dev.off()

			# Plot of \hat \pi vs. \lambda
			postscript(paste(file,"_pihat_vs_lambda_",cnames[i],".eps",sep=""))
			qplot2(qobj,rng=c(0,1),whichplot=1)
			dev.off()

			# Plot of significant tests vs. q cut-off
			postscript(paste(file,"_sigtest_vs_qval_",cnames[i],".eps",sep=""))
			qplot2(qobj,rng=c(0,1),whichplot=3)
			dev.off()
		}
	}
	list(qval.mat=qval.mat,pi0=pi0)
}





q.interest <- function(qval.mat,alpha=0.15, criteria = c("more","less")[2]){

	interest <- matrix(0,nrow=nrow(qval.mat),ncol=ncol(qval.mat))
	interest <- as.data.frame(interest)
	rownames(interest) <- rownames(qval.mat)
	colnames(interest) <- colnames(qval.mat)

	for(i in 1:ncol(interest)){
		qval <- qval.mat[,i]

		if(criteria == "less"){
			ind <- which(qval <= alpha)
		} else {
			ind <- which(qval > alpha)
		}
		if(length(ind)>0){
			interest[ind,i] <- 1
		}
	}
	list(interest=interest)
}

# Function to find minimum qvalues

minqval <- function(out,method=c("smoother","bootstrap")[2]){

    interest <- matrix(0,nrow=ncol(out$pvalues),ncol=1)
    interest <- as.data.frame(interest)
    rownames(interest) <- colnames(out$pvalues)

    for(i in 1:nrow(interest)){
        #print(i)
        pvalues <- out$pvalues[,i]
        qobj <- qvalue.adj(pvalues,pi0.method=method,lambda=seq(0,0.90,by=0.01))
        #qobj <- qvalue.adj(pvalues,pi0.method=method)
        qval <- qobj$qvalues
        interest[i,] <- min(qval)
    }
    return(interest)
}


##########################
# Bonferonni-Holm Method #
##########################


bon.holm.interest <- function(pvalues,alpha=0.05){

    interest <- matrix(0,nrow=nrow(pvalues),ncol=ncol(pvalues))
    interest <- as.data.frame(interest)
    rownames(interest) <- rownames(pvalues)
    colnames(interest) <- colnames(pvalues)

    pval.adjust <- matrix(0,nrow=nrow(pvalues),ncol=ncol(pvalues))
    pval.adjust <- as.data.frame(pval.adjust)
    rownames(pval.adjust) <- rownames(pval.adjust)
    colnames(pval.adjust) <- colnames(pval.adjust)

    for(i in 1:ncol(interest)){
        pval <- pvalues[,i]

        # adjust p-values using bonferroni-holm method
        pval.adjust[,i] <- p.adjust(pval,method="holm")

        ind <- which(pval.adjust[,i] <= alpha)

        if(length(ind)>0){
            interest[ind,i] <- 1
        }
    }
    list(interest=interest,pval.adjust=pval.adjust)
}

#############################
# Benjamini-Hochberg Method #
#############################
# EXPORT
#' @importFrom stats p.adjust
ben.hoch.interest <- function(pvalues,alpha=0.05){

  interest <- matrix(0,nrow=nrow(pvalues),ncol=ncol(pvalues))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(pvalues)
  colnames(interest) <- colnames(pvalues)

  pval.adjust <- matrix(0,nrow=nrow(pvalues),ncol=ncol(pvalues))
  pval.adjust <- as.data.frame(pval.adjust)
  rownames(pval.adjust) <- rownames(pval.adjust)
  colnames(pval.adjust) <- colnames(pval.adjust)

  for(i in 1:ncol(interest)){
    pval <- pvalues[,i]

    ## adjust p-values using Benjamini-Hochberg method
    pval.adjust[,i] <- stats::p.adjust(pval,method="BH")

    ind <- which(pval.adjust[,i] <= alpha)

    if(length(ind)>0){
      interest[ind,i] <- 1
    }
  }
  return(list(interest=interest,pval.adjust=pval.adjust))
}



#############################
## weighted LASSO Approach ##
#############################


# function to standardize a vector x
make.std <- function(x){
	N <- length(x)
	( x-mean(x) ) / ( sd(as.vector(x)) * sqrt( N / (N-1) ) )
}

make.center <- function(x){
	return(x-mean(x))
}

lasso <- function(weights,yy,XX,data.delta,g,file="file",plots=FALSE,include.diet=TRUE,
                  diet.wt=1000,thresh.q=FALSE,corr.g=FALSE,delta=2,std.y=TRUE,
                  est.MSE=c("TRUE","est.var","step")[2],
                  cv.criterion=c(FALSE,"delta_cv")[1],vfold=10){
    #print("lasso")
    #print(yy)
    y <- t(yy)
    #print(y)
    X <- t(XX)
    if(!is.null(data.delta)){
        data.delta <- t(data.delta)
    }
    N <- length(y)
    y1 <- y
    #  if(std.y==TRUE){
    #    ## Standardize y
    #    y1 <- make.std(y)
    #  } else {
    #    ## Centered response
    #    y1 <-  make.center(y)
    #  }
    # Standardized design matrix X
    X1 <- apply(X,2,make.std)
    weights.use <- g(weights)
    # Get rid of small weights
    if(thresh.q==TRUE){
        for(i in 1:length(weights)){
            weights.use[i] <- max(0.0001,abs(weights.use)[i])
        }
    }
    # Be sure to include diet?
    if(include.diet==TRUE){
        wts <- c(1e-5,weights.use)
    } else {
        wts <- c(1,weights.use)
    }
    X1 <- X1 %*% diag(1/wts)
    if(is.null(data.delta)){
        ## Fix use.Gram
        if(ncol(X1)>500){
            use.Gram <- FALSE
        } else {
            use.Gram <- TRUE
        }
        ## Run Lasso
        wLasso.out <-  lars(X1, y1, type = c("lasso"),
                            trace = FALSE, normalize = FALSE, intercept = FALSE, use.Gram = use.Gram)
        ## entry in Lasso
        entry.variables <- as.numeric(wLasso.out$entry)
        ##order.variables <- as.numeric(wLasso.out$actions)
        order.variables <- 0
        if(cv.criterion==TRUE){
            ## good implementation, mode="step"
            wLasso.cv <- cv.lars(X1, y1, type = c("lasso"), K = vfold,
                                 trace = FALSE, normalize = FALSE, intercept= FALSE,
                                 plot.it=FALSE,mode="step",use.Gram=use.Gram)
            bestindex <- wLasso.cv$index[which.min(wLasso.cv$cv)]
            ## final best descriptive model
            predict.out <- predict(wLasso.out, X1,s=bestindex, type = "coefficients", mode="step")
            delta.out <- delta
        } else if(cv.criterion=="delta_cv"){
            cv.out <- cv.delta(y1,X1,K=vfold,est.MSE=est.MSE)
            predict.out <- cv.out$predict.out
            delta.out <- cv.out$delta
        } else {
            ## Samuel's corrections 11/9/2011
            ## use Cp-like criterion to find best descriptive model
            p = dim(X1)[2]
            s = length(wLasso.out$df)
            p.pos = NULL
            RSS = NULL
            for (i in 1:s){
                RSS[i] = sum((y1-predict(wLasso.out, X1, s=i, type = c("fit"))$fit)**2)
                p.pre = predict(wLasso.out, X1, s=i, type = c("coefficients"))$coefficients
                p.pos = c(p.pos,length(p.pre[abs(p.pre)>0]))
            }
            ## Get estimated MSE
            if(est.MSE=="TRUE"){
                MSE <- 0.5
                ##print(MSE)
            } else if(est.MSE=="est.var") {
                MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
                MSE <- MSE^2
                ##print(MSE)
            } else {
                ## Get estimated MSE from best fit of forward stepwise regression with AIC as selection criterion
                X2 <- data.frame(X1)
                colnames(X2)[1] <- "Fixed"
                full.lm <- lm(y1~.,data=X2)
                start.lm <- lm(y1~-1 + Diet,data=X2)
                lowest.step.forward <- step(lm(y1 ~ -1+Diet, data=X2),
                                            list(lower=start.lm,upper=full.lm), direction='forward',trace=FALSE)
                MSE <- summary(lowest.step.forward)$sigma
                ##print(MSE)
                MSE <- summary(lowest.step.forward)$sigma^2
                ##	print(MSE)
                if(MSE < 1e-5){
                    MSE <-  sd(as.vector(y1)) * sqrt( N / (N-1) )
                    MSE <- MSE^2
                }
            }
            p.min = which.min(RSS/MSE+delta*p.pos)
            ## final best descriptive model
            predict.out <- predict(wLasso.out, X1, s=p.min, type = c("coefficients"))
            delta.out <- delta
        }
    } else {
        if(include.diet==TRUE){
            penalty <- c(0,rep(1,ncol(X1)-1))
        } else {
            penalty <- rep(1,ncol(X1))
        }
        ## Run Lasso
        #ytmp <- cbind(time=y1,status=data.delta)
        ytmp <- cbind(time=t(y1),status=t(data.delta))
        #print(ytmp)
        colnames(ytmp) <- c("time","status")
        ## entry in Lasso
        entry.variables <- 0 ## not available
        order.variables <- 0 ## not available

        if(cv.criterion==TRUE){
            wLasso.cv <- cv.glmnet(X1, ytmp, standardize=FALSE,family="cox",alpha=1,penalty.factor=penalty)
            lambda.opt <- wLasso.cv$lambda.min
        } else if(cv.criterion=="delta_cv"){
            ## not used
            # cv.out <- cv.delta(y1,X1,K=vfold,est.MSE=est.MSE)
            # predict.out <- cv.out$predict.out
            # delta.out <- cv.out$delta
        } else {
            wLasso.out <-  glmnet(X1, ytmp, standardize=FALSE,family="cox",alpha=1,penalty.factor=penalty)
            ## BIC/Deviance criterion: deviance + k*log(n)
            deviance <- deviance(wLasso.out)
            p.min <- which.min(deviance + wLasso.out$df*log(N))
            lambda.opt <- wLasso.out$lambda[p.min]
        }
        wLasso.fit <- glmnet(X1, ytmp, standardize=FALSE,family="cox",alpha=1,lambda=lambda.opt)
        ## final best descriptive model
        predict.out <- list(coefficients=wLasso.fit$beta)
        delta.out <- delta
    }
    ind <- which(as.logical(abs(predict.out$coefficients)>1e-10))
    ##ind <- which(predict.out$coefficients!=0)
    sig.variables <- rep(0,nrow(XX))
    sig.variables[ind] <- 1
    sign.of.variables <- rep(0,nrow(XX))
    ind.pos <- which(as.logical(predict.out$coefficients >0))
    sign.of.variables[ind.pos] <- 1
    ind.neg <- which(as.logical(predict.out$coefficients <0))
    sign.of.variables[ind.neg] <- -1
    if(plots==TRUE){
        postscript(paste(file,"_lasso1.eps",sep=""))
        plot(wLasso.out,cex.axis=1.5,cex.lab=1.5)
        dev.off()
        postscript(paste(file,"_lasso2.eps",sep=""))
        ##x11()
        ##plot(s,RSS+2*(s),type="l",cex.axis=1.5,cex.lab=1.5)
        ##abline(v=s.min,lty=2)
        ## Samuel's correction 11/9/2011
        par(mar=c(5, 4, 4, 2)+1)
        plot(1:s,RSS+2*(p.pos),type="l",cex.axis=1.5,cex.lab=1.5,ylab=substitute(M[n](that,p),
                                                                                 list(that=delta)), xlab="Steps")
        abline(v=p.min,lty=2)
        dev.off()
    }
    list(order.variables=order.variables,sig.variables=sig.variables,
         sign.of.variables=sign.of.variables,entry.variables=entry.variables,delta.out=delta.out)
}

lasso.computations <- function(weights,XX,response,g,plots=TRUE,file="name",include.diet=TRUE,
                               diet.wt=100,thresh.q=FALSE,corr.g=FALSE,delta=2,std.y="TRUE",
                               est.MSE=c("TRUE","est.var","step")[2],
                               cv.criterion=FALSE,vfold=10){
    #print(response)
    data.response <- response$yy
    data.delta <- response$delta

    interest <- matrix(0,nrow=nrow(XX),ncol=nrow(data.response))
    interest <- as.data.frame(interest)
    rownames(interest) <- rownames(XX)
    colnames(interest) <- rownames(data.response)

    interest.sign <- interest		# matrix to store sign of lasso coefficients

    R2.val <- matrix(0,nrow=nrow(data.response),ncol=1)
    R2.val <- as.data.frame(R2.val)
    colnames(R2.val) <- "R2"
    rownames(R2.val) <- rownames(data.response)

    order.var <- array(0,dim=c(nrow(data.response),300))

    entry.var <- array(0,dim=c(nrow(data.response),nrow(XX)))

    delta.out <- matrix(0,nrow=1,ncol=nrow(data.response))

    for(i in 1:ncol(interest)){
        #print(data.response)
        lasso.out <- lasso(weights[,i],data.response[i,],XX,data.delta[i,],g,
                           file=paste(file,rownames(data.response)[i],sep=""),
                           plots=plots,include.diet=include.diet,
                           diet.wt=diet.wt,thresh.q=thresh.q,corr.g=corr.g,delta=delta,std.y=std.y,
                           est.MSE=est.MSE,
                           cv.criterion=cv.criterion,vfold=vfold)
        interest[,i] <- lasso.out$sig.variables
        interest.sign[,i] <- lasso.out$sign.of.variables
        R2.val[i,]   <- get.R2(interest[,i],XX,data.response[i,])

        order.var[i,1:length(lasso.out$order.variables)] <- lasso.out$order.variables
        entry.var[i,] <- lasso.out$entry.variables

        delta.out[,i] <- lasso.out$delta.out

    }
    list(interest=interest,order.var=order.var,R2.val=R2.val,interest.sign=interest.sign,
         entry.var=entry.var,delta.out=delta.out)
}



###############################################
## Cross-Validation function to select delta ##
###############################################

cv.delta <- function(y1,X1,K=10,est.MSE=c("TRUE","est.var","step")[2]){
    ## sequence of delta values
    delta.cv <- seq(0.75,2,by=0.1)
    ##delta.cv <- seq(0.1,2,by=0.1)

    ## Randomly partition the data
    all.folds <- cv.folds(length(y1), K)

    ## Matrix to store residuals
    residmat <- matrix(0, length(delta.cv), K)

    for(j in seq(K)){
        ## data set to omit
        omit <- all.folds[[j]]

        ## Run Lasso with after removing omit data
        wLasso.out <- lasso.procedure(y1[-omit],X1[-omit,,drop=FALSE])$wLasso.out

        for(d in 1:length(delta.cv)){

            ## Find best-fitting model for specified delta
            beta.omit <- lasso.delta.choice(wLasso.out,y1[-omit],X1[-omit,,drop=FALSE],delta=delta.cv[d],est.MSE=est.MSE)

            ## Find final fit with data omitted
            fit <- X1[omit,,drop=FALSE] %*% beta.omit$predict.out$coefficients

            ## Store residual
            residmat[d,j] <-  apply((y1[omit] - fit)^2, 2, sum)

        }

    }

    cv <- apply(residmat,1,mean)

    ## Check which delta's lead to min(cv)
    delta.ind <- which(cv==min(cv))
    delta.opt <- delta.cv[delta.ind]

    delta <- mean(delta.opt)		## takes average of delta values

    wLasso.out <- lasso.procedure(y1,X1)$wLasso.out
    predict.out <- lasso.delta.choice(wLasso.out,y1,X1,delta=delta,est.MSE=est.MSE)$predict.out
    list(predict.out=predict.out,delta=delta)
}

lasso.procedure <- function(y1,X1){
	# Setup
	N <- length(y1)

	## adjust use.Gram
	if(ncol(X1)>500){
		use.Gram <- FALSE
	} else {
		use.Gram <- TRUE
	}


	# Run Lasso
	wLasso.out <-  lars(X1, y1, type = c("lasso"),
                trace = FALSE, normalize = FALSE, intercept = FALSE,use.Gram=use.Gram)

	list(wLasso.out=wLasso.out)
}

lasso.delta.choice <- function(wLasso.out,y1,X1,delta,est.MSE=c(TRUE,"est.var","step")[2]){
    # Setup
    N <- length(y1)

    p = dim(X1)[2]

    s = length(wLasso.out$df)
    p.pos = NULL

    RSS = NULL
    for (i in 1:s){
        RSS[i] = sum((y1-predict(wLasso.out, X1, s=i, type = c("fit"))$fit)**2)
        p.pre = predict(wLasso.out, X1, s=i, type = c("coefficients"))$coefficients
        p.pos = c(p.pos,length(p.pre[abs(p.pre)>0]))
    }

    ## Get estimated MSE
    if(est.MSE=="TRUE"){
        MSE <- 0.5
        ##print(MSE)
    } else if(est.MSE=="est.var") {
        MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
        MSE <- MSE^2
        ##print(MSE)
    } else {
        ## Get estimated MSE from best fit of forward stepwise regression with AIC as selection criterion
        X2 <- data.frame(X1)
        colnames(X2)[1] <- "Fixed"
        full.lm <- lm(y1~.,data=X2)
        start.lm <- lm(y1~-1 + Diet,data=X2)
        lowest.step.forward <- step(lm(y1 ~ -1+Diet, data=X2),
                                    list(lower=start.lm,upper=full.lm), direction='forward',trace=FALSE)
        MSE <- summary(lowest.step.forward)$sigma
        ##print(MSE)
        MSE <- summary(lowest.step.forward)$sigma^2
        ##	print(MSE)
        if(MSE < 1e-5){
            MSE <-  sd(as.vector(y1)) * sqrt( N / (N-1) )
            MSE <- MSE^2
        }
    }

    p.min = which.min(RSS/MSE+delta*p.pos)


    ## final best descriptive model
    predict.out <- predict(wLasso.out, X1, s=p.min, type = c("coefficients"))

    list(wLasso.out=wLasso.out,predict.out=predict.out)
}

# Function to get R^2 value
get.R2 <- function(sig,XX,response){
    if(sum(sig)==0){
        return("NA")
    } else {
        Xtmp <- t(XX[which(sig==1),])
        if(length(Xtmp)==dim(Xtmp)[2]){
            X <- make.std(as.numeric(Xtmp))
        } else {
            X <- apply(Xtmp,2,make.std)
        }

        y <- make.std(as.numeric(response))
        out <- summary(lm(y~X))
        R2 <- out$r.squared
        return(R2)
    }
}



#################################
## Exclusion-frequency weights ##
#################################

##############################################
## Function to determine size of partitions ##
##############################################

partition.size <- function(n,p,k){
    lo <- floor(p/k)
    hi <- ceiling(p/k)

    constraints <- c(lo,hi)

    ## We require that q_l < n for each l
    index <- which(constraints<n)

    if(length(index)< 2){
        return("k too low")
    }

    constraints <- constraints[index]


    ## Each q_l is such that : lo <= q_l <= hi
    #tmp <- t(matrix(rep(constraints,k),nrow=k,ncol=2,byrow=TRUE))

    ## Determining different q_l values
    #combinations <- as.matrix(do.call(`expand.grid`,as.data.frame(tmp)))

    ## Determine which combinations are such that \sum q_l =p
    #combinations <- combinations[which(apply(combinations,1,sum)==p),]
    #groups <- as.vector(combinations[1,])

    combinations <- rep(floor(p/k),k)
    if(sum(combinations)!=p){
        tmp.diff <- p-sum(combinations)
        combinations[1:tmp.diff] <- combinations[1:tmp.diff]+1
    }
    groups <- combinations

    return(groups)
}

## determine partition size without requiring q_l < n
partition.size.new <- function(p,k){
    lo <- floor(p/k)
    hi <- ceiling(p/k)

    constraints <- c(lo,hi)

    ## Each q_l is such that : lo <= q_l <= hi
    #tmp <- t(matrix(rep(constraints,k),nrow=k,ncol=2,byrow=TRUE))

    ## Determining different q_l values
    #combinations <- as.matrix(do.call(`expand.grid`,as.data.frame(tmp)))

    ## Determine which combinations are such that \sum q_l =p
    #combinations <- combinations[which(apply(combinations,1,sum)==p),]
    #groups <- as.vector(combinations[1,])

    combinations <- rep(floor(p/k),k)
    if(sum(combinations)!=p){
        tmp.diff <- p-sum(combinations)
        combinations[1:tmp.diff] <- combinations[1:tmp.diff]+1
    }
    groups <- combinations

    return(groups)
}


###############################################
## Function to randomly partition data index ##
###############################################
random.partition <- function(n,p,k){
    ## Size of each partition group
    if(p!=k){
        size.groups <- partition.size(n,p,k)
    } else {
        size.groups <- rep(1,k)
    }

    ## Name of each group
    names <- paste("group",seq(1,k),sep="")
    cut.by <- rep(names, times = size.groups)

    ## Get random index
    rand.index <- split(sample(1:p, p), cut.by)

    return(rand.index)
}

index.sort.partition <- function(n,k,beta.sort.ix){
    p <- length(beta.sort.ix)

    ## Size of each partition group
    size.groups <- partition.size(n,p,k)

    ## Name of each group
    names <- seq(1,k)
    index.sort <- rep(names, times = size.groups)

    names(index.sort) <- paste("X_",beta.sort.ix,sep="")
    return(index.sort)

}

## get sizes of fixed covariates partition
fixed.size <- function(fixed.covariates,k){
    whole.groups <- floor(length(fixed.covariates)/k)
    group.size <- rep(whole.groups,k)

    tmp.diff <- abs(sum(group.size)-length(fixed.covariates))

}

## Function to ensure certain covariates are in each group, and the rest are randomly placed
fixed.plus.random.partition <- function(fixed.covariates,n,p,k){
    ## Size of each partition group: minus fixed.covariates
    size.groups <- partition.size(n,p-length(fixed.covariates),k)

    ## Name of each group
    names <- paste("group",seq(1,k),sep="")
    cut.by <- rep(names, times = size.groups)

    ## Get all covariates that will be randomly sampled (i.e., remove fixed.covariates)
    all <- 1:p
    covariates.not.fixed <- all[!all%in%fixed.covariates]

    ## Get random index for covariates.not.fixed
    rand.index.covariates.not.fixed <- split(sample(covariates.not.fixed, length(covariates.not.fixed)), cut.by)

    ## Put fixed covariates back in
    fixed.cut.by <-rep(names,times=partition.size.new(length(fixed.covariates),k))
    fixed.covariates.partition <- split(fixed.covariates,fixed.cut.by)

    rand.index <- appendList(rand.index.covariates.not.fixed,fixed.covariates.partition)

    return(rand.index)
}

## function to append lists
appendList <- function (x, val){
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    for (v in names(val)) {
        x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
            appendList(x[[v]], val[[v]])
        else c(x[[v]], val[[v]])
    }
    x
}

## Ridge regresssion
ridge.regression <- function(XX,response){
    ##mydata <- data.frame(cbind(t(response),t(XX)))
    #### center data
    ##mydata <- apply(mydata,2,make.center)
    ##xnam.orig <- paste("X_",1:(nrow(XX)-1),sep="")
    ##xnam <- c("Fixed",xnam.orig)
    ##fmla <- as.formula(paste("response~0+",paste(xnam,collapse="+")))
    ##out <-lm.ridge(fmla,data=mydata,lambda=seq(0.00001,25,0.1))

    ## Center data
    X <- t(XX)
    X <- apply(X,2,make.center)

    yy <- response$yy
    delta <- response$delta

    if(is.null(delta)){
        #######################
        ## linear regression ##
        #######################
        y <- t(yy)
        ##y <- apply(y,2,make.center)
        family <- "gaussian"
        grouped <- FALSE
    } else {
        yy <- as.numeric(yy)
        delta <- as.numeric(delta)
        y <- cbind(time=yy,status=delta)
        colnames(y) <- c("time","status")
        family <- "cox"
        grouped <- TRUE
    }

    ## Find optimal lambda value
    cv.out <- cv.glmnet(X,y,standardize=FALSE,family=family,alpha=0,grouped=grouped)
    lambda.opt <- cv.out$lambda.min

    ## Ridge regression
    ridge.out <- glmnet(X,y,standardize=FALSE,family=family,alpha=0,lambda=lambda.opt)
    beta.values <- ridge.out$beta[-1,]  ## remove diet

    return(beta.values)
}

## We partition based on ascending values of ridge regression estimates
ridge.sort.partition <- function(k,XX,response){
    ## Get parameter estimates from ridge regression
    beta.values <- ridge.regression(XX,response)

    ## Sort beta.values in decreasing order
    beta.sort <- sort(abs(beta.values),decreasing=TRUE,index.return=TRUE)

    ## index of ordering
    beta.index <- beta.sort$ix

    ## Partition based on decreasing order
    rand.index <- sort.partition(beta.index,k)

    list(rand.index=rand.index,beta.index=beta.index)
}



## We partition based on sorting the ridge regression estimates
sort.partition <- function(beta.index,k){
    ## Name of each group
    names <- paste("group",seq(1,k),sep="")

    ## Get rand.index labels
    rand.index <- vector("list",length(names))
    names(rand.index) <- names

    ## Sort indices for each group
    tmp <-1
    size.groups <- NULL
    tmp.index <- NULL
    for(j in 1:k){
        indices <- seq(tmp,length(beta.index),by=k)
        size.groups <- c(size.groups,length(indices))
        tmp.index <- c(tmp.index,indices)
        tmp <- tmp + 1
    }

    ## Put groups together
    cut.by <- rep(names, times = size.groups)

    ## Get random index
    rand.index <- split(beta.index[tmp.index], cut.by)

    return(rand.index)
}



## Function to do designed partitioning
## index.group : variable designating which variable belongs to which group
designed.partition <- function(index.group,k){
    cluster <- index.group

    ## Name of each group
    names <- paste("group",seq(1,k),sep="")

    ## Get k-means partition of beta values
    rand.index <- vector("list",length(names))
    names(rand.index) <- names

    for(j in 1:k){
        cluster.tmp <- as.numeric(which(cluster==j))
        size.groups <- partition.size.new(length(cluster.tmp),k)
        cut.by <- rep(names, times = size.groups)
        rand.index.tmp <- split(sample(cluster.tmp,length(cluster.tmp)), cut.by)
        ##rand.index <- mapply(c,rand.index,rand.index.tmp,SIMPLIFY=FALSE)
        rand.index <- appendList(rand.index,rand.index.tmp)
    }

    ## Below won't work
    ##index.val <-1:length(index.group)
    ##mydf <- data.frame(index.val,index.group)
    ##rand.index <- vector("list",length(names))
    ##names(rand.index) <- names

    ##size.big.groups <- partition.size(n=ncol(XX),p=nrow(XX)-1,k=k)
    ##for(s in 1:k){
    ##  size.cluster <- partition.size.new(size.big.groups[s],k)
    ## tmp <- stratified(mydf,"index.group",size=size.cluster)
    ##  selected.group <- tmp[,"index.val"]
    ##  rand.index[[s]] <- selected.group

    ## Update mydf
    ## mydf <- mydf[-selected.group,]
    ##}

    return(rand.index)
}


columns.to.list <- function( df ) {
    ll<-apply(df,2,list)
    ll<-lapply(ll,unlist)
    return(ll)
}


## Function to do step AIC on group subset
step.selection <- function(factor.z,index,XX,response,type=c("AIC","BIC")[1],
                           direction=c("both","forward","backward")[3],real_data){

    if(real_data==FALSE){
        xnam.orig <- paste("X_",index,sep="")
    } else {
        xnam.orig <- rownames(XX[index+1,])
    }

    if(factor.z==TRUE){
        xnam <- c("factor(Fixed)",xnam.orig)
    } else {
        xnam <- c("Fixed",xnam.orig)
    }

    yy <- response$yy
    delta <- response$delta

    if(is.null(delta)){
        #######################
        ## linear regression ##
        #######################
        mydata <- data.frame(cbind(t(yy),t(XX)))
        fmla <- as.formula(paste("response~",paste(xnam,collapse="+")))
        fit <- lm(fmla,data=mydata)
    } else {
        #########################
        ## survival regression ##
        #########################
        X <- data.frame(t(XX))
        tmp.list <- columns.to.list(X)
        mydata <- list(time=as.numeric(yy),status=as.numeric(delta))
        mydata <- appendList(mydata,tmp.list)
        fmla <- as.formula(paste("Surv(time,status)~",paste(xnam,collapse="+")))
        fit <- coxph(fmla,data=mydata)
    }
    ##print(fmla)


    if(type=="AIC"){
        deg.free <- 2
    } else {
        deg.free <- log(length(yy))
    }

    if(factor.z==TRUE){
        step.reg <- step(fit,k=deg.free,direction=direction,
                         scope = list(lower = ~ factor(Fixed)),trace=FALSE,data=mydata)
    } else {
        step.reg <- step(fit,k=deg.free,direction=direction,
                         scope = list(lower = ~ Fixed),trace=FALSE,data=mydata)
    }

    ## XX selected
    results <- intersect(names(step.reg$coefficients),xnam.orig)

    ## Store results
    out <- data.frame(rep(0,nrow(XX)-1))
    rownames(out) <- rownames(XX)[-1]
    colnames(out) <- "results"
    out[xnam.orig,] <- as.numeric(!xnam.orig%in%results)
    return(out)
}




## Function to do step AIC on group subset
step.selection.inclusion <- function(factor.z,index,XX,response,type=c("AIC","BIC")[1],
                                     direction=c("both","forward","backward")[3],real_data){

    if(real_data==FALSE){
        xnam.orig <- paste("X_",index,sep="")
    } else {
        xnam.orig <- rownames(XX[index+1,])
    }

    if(factor.z==TRUE){
        xnam <- c("factor(Fixed)",xnam.orig)
    } else {
        xnam <- c("Fixed",xnam.orig)
    }

    yy <- response$yy
    delta <- response$delta

    if(is.null(delta)){
        #######################
        ## linear regression ##
        #######################
        mydata <- data.frame(cbind(t(yy),t(XX)))
        fmla <- as.formula(paste("response~",paste(xnam,collapse="+")))
        fit <- lm(fmla,data=mydata)
    } else {
        #########################
        ## survival regression ##
        #########################
        X <- data.frame(t(XX))
        tmp.list <- columns.to.list(X)
        mydata <- list(time=as.numeric(yy),status=as.numeric(delta))
        mydata <- appendList(mydata,tmp.list)
        fmla <- as.formula(paste("Surv(time,status)~",
                                 paste(xnam,collapse="+")))
        fit <- coxph(fmla,data=mydata)
    }
    ##print(fmla)
    ## go here for forward selection example: http://msekce.karlin.mff.cuni.cz/~pesta/NMFM404/ph.html



    if(type=="AIC"){
        deg.free <- 2
    } else {
        deg.free <- log(length(yy))
    }

    if(factor.z==TRUE){
        step.reg <- step(fit,k=deg.free,direction=direction,
                         scope = list(lower = ~ factor(Fixed)),trace=FALSE,data=mydata)
    } else {
        step.reg <- step(fit,k=deg.free,direction=direction,
                         scope = list(lower = ~ Fixed),trace=FALSE,data=mydata)
    }

    ## XX selected
    results <- intersect(names(step.reg$coefficients),xnam)

    ## Store results
    out <- data.frame(rep(0,nrow(XX)))
    rownames(out) <- rownames(XX)
    colnames(out) <- "results"
    out[xnam,] <- as.numeric(xnam%in%results)
    return(out)
}





