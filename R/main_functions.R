####################################
####################################
##
##
## Functions for the main method
##
####################################
####################################

#' Weigted lasso variable selection
#'
#' Performs variable selection with covariates multiplied by weights that direct which variables
#' are likely to be associated with the response.
#'
#' @param x (n by m) matrix of main covariates where m is the number of covariates and n is the sample size.
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable. This covariate should
#' always be selected. Can be NULL.
#' @param y (n by 1) a matrix corresponding to the response variable. If \code{regression.type} is "cox",
#' \code{y} contains the observed event times.
#' @param cox.delta (n by 1) a matrix that denotes censoring when \code{regression.type} is "cox" (1 denotes
#' survival event is observed, 0 denotes the survival event is censored). Can be NULL.
#' @param factor.z logical. If TRUE, the fixed variable z is a factor variable.
#' @param regression.type a character indicator that is either "linear" for linear regression
#' or "cox" for Cox proportional hazards regression. Default is "linear".
#' @param weight.type Character value denoting which weights to be used for the weighted lasso, where each covariate in \code{x}
#' is multiplied by a scalar weight. Options include
#' \itemize{
#'   \item{one:}{The scalar weight is one.}
#'   \item{corr.estimate:}{The scalar weight for covariate \eqn{x_j} is the Pearson correlation between
#'   \eqn{x_j} and \eqn{y}. }
#'   \item{corr.pvalue:}{The scalar weight for covariate \eqn{x_j} is the p-value of the coefficient of
#'   \eqn{x_j} in the regression of \eqn{y} on \eqn{x_j} }
#'   \item{corr.bh.pvalue:}{The scalar weight for covariate \eqn{x_j} is the Benjanmini-Hocbherg adjusted p-value
#'   from \code{corr.pvalue}.}
#'   \item{corr.qvalue:}{The scalar weight for covariate \eqn{x_j} is the q-value transform of the p-value
#'   from \code{corr.pvalue}.}
#'   \item{corr.tstat:}{The scalar weight for covariate \eqn{x_j} is the t-statistic associated with testing
#'   the significance of \eqn{x_j} in the regression of \eqn{y} on \eqn{x_j}.}
#'   \item{parcor.estimate:}{The scalar weight for covariate \eqn{x_j} is the partial correlation between
#'   \eqn{x_j} and \eqn{y} after adjustment for \eqn{z}. }
#'   \item{parcor.pvalue:}{The scalar weight for covariate \eqn{x_j} is the p-value of the coefficient of
#'   \eqn{x_j} in the regression of \eqn{y} on \eqn{z} and \eqn{x_j} }
#'   \item{parcor.bh.pvalue:}{The scalar weight for covariate \eqn{x_j} is the Benjanmini-Hocbherg adjusted p-value
#'   from \code{parcor.pvalue}.}
#'   \item{parcor.qvalue:}{The scalar weight for covariate \eqn{x_j} is the q-value transform of the p-value
#'   from \code{parcor.pvalue}.}
#'   \item{parcor.tstat:}{The scalar weight for covariate \eqn{x_j} is the t-statistic associated with testing
#'   the significance of \eqn{x_j} in the regression of \eqn{y} on \eqn{z} and \eqn{x_j}.}
#'   \item{exfrequency.random.partition.aic:}{The scalar weight for covariate \eqn{x_j} is an exclusion frequency.
#'   The exclusion frequency is obtained as follows: we first partition the covariates into \code{k.split}
#'   random groups, and we apply a stepwise linear/Cox regression of the response on each partition
#'   set of covariate. The final model is selected using an AIC criterion, and we track if \eqn{x_j} is
#'   excluded from the final model. We repeat
#'   this procedure \code{nboot} times and the exclusion frequency is the average number of times
#'   \eqn{x_j} is excluded.}
#'   \item{exfrequency.random.partition.bic:}{The scalar weight for covariate \eqn{x_j} is computed
#'   as in exfrequency.random.partition.aic, except that the final model within each stepwise regression
#'   is selected using a BIC criterion.}
#'   \item{exfrequency.kmeans.partition.aic:}{The scalar weight for covariate \eqn{x_j} is an exclusion frequency.
#'   The exclusion frequency is obtained as follows: we apply ridge regression of the response on all covariates and
#'   obtain ridge regression coefficients for each covariate.
#'   We then partitioned the covariates into \code{k.split}
#'   groups using a K-means criterion on the ridge regression coefficients, and we applied
#'   a stepwise linear/Cox regression of the response on each partition
#'   set of covariate. The final model is selected using an AIC criterion, and we track if \eqn{x_j} is
#'   excluded from the final model. We repeat
#'   this procedure \code{nboot} times and the exclusion frequency is the average number of times
#'   \eqn{x_j} is excluded.}
#'   \item{exfrequency.kmeans.partition.bic:}{The scalar weight for covariate \eqn{x_j} is computed
#'   as in exfrequency.kmeans.partition.aic, except that the final model within each stepwise regression
#'   is selected using a BIC criterion.}
#'   \item{exfrequency.kquartile.partition.aic:}{The scalar weight for covariate \eqn{x_j} is an exclusion frequency.
#'   The exclusion frequency is obtained as follows: we apply ridge regression of the response on all covariates and
#'   obtain ridge regression coefficients for each covariate.
#'   We then partitioned the covariates into \code{k.split}
#'   groups using k-quantiles of the ridge regression coefficients, and we applied a stepwise linear/Cox regression of the response on each partition
#'   set of covariate. The final model is selected using an AIC criterion, and we track if \eqn{x_j} is
#'   excluded from the final model. We repeat
#'   this procedure \code{nboot} times and the exclusion frequency is the average number of times
#'   \eqn{x_j} is excluded.}
#'   \item{exfrequency.kquartiles.partition.bic:}{The scalar weight for covariate \eqn{x_j} is computed
#'   as in exfrequency.kquartiles.partition.aic, except that the final model within each stepwise regression
#'   is selected using a BIC criterion.}
#'   \item{exfrequency.ksorted.partition.aic:}{The scalar weight for covariate \eqn{x_j} is an exclusion frequency.
#'   The exclusion frequency is obtained as follows: we apply ridge regression of the response on all covariates and
#'   obtain ridge regression coefficients for each covariate.
#'   We then partitioned the covariates into \code{k.split}
#'   groups by first ordering the ridge regression coefficients in descending order
#'   and splitting them into \code{k.split} groups. We then applied a stepwise linear/Cox regression of the response on each partition
#'   set of covariate. The final model is selected using an AIC criterion, and we track if \eqn{x_j} is
#'   excluded from the final model. We repeat
#'   this procedure \code{nboot} times and the exclusion frequency is the average number of times
#'   \eqn{x_j} is excluded.}
#'   \item{exfrequency.ksorted.partition.bic:}{The scalar weight for covariate \eqn{x_j} is computed
#'   as in exfrequency.ksorted.partition.aic, except that the final model within each stepwise regression
#'   is selected using a BIC criterion.}
#' }
#' @param weight_fn A user-defined function to be applied to the weights for the weighted lasso.
#' Default is an identify function.
#' @param ttest.pvalue logical indicator used when \code{weight.type} is "corr.pvalue","corr.bh.pvalue", "corr.qvalue",
#' "parcor.pvalue","parcor.bh.pvalue","parcor.qvalue". If TRUE, p-value for each covariate is computed from univariate
#' linear/cox regression of the response on each covariate. If FALSE, the
#' p-value is computed from correlation coefficients between the response and each covariate.
#' Default is FALSE.
#' @param q_opt_tuning_method character indicator used when \code{weight.type} is "corr.qvalue" or "parcor.qvalue".
#' Options are "bootstrap" or "smoother" to specify how the optimal tuning parameter is obtained when computing
#' q-values from Storey and Tibshirani (2003). Default is "smoother" (smoothing spline).
#' @param robust logical indicator used when \code{weight.type} is "corr.qvalue" or "parcor.qvalue".
#' If TRUE, q-values computed as in Storey and Tibshirani (2003)
#' are robust for small p-values.
#' @param pi0.known logical indicator used when \code{weight.type} is "corr.qvalue" or "parcor.qvalue".
#' If TRUE, when computing q-values, the estimate of the
#' true proportion of the null hypothesis is set to the value of pi0.val given by the user.
#' If FALSE, the estimate of the true proportion of the null hypothesis is
#' computed by bootstrap or smoothing spline as proposed in Storey and Tibshirani (2003). Default is FALSE.
#' @param pi0.val scalar used when \code{weight.type} is "corr.qvalue" or "parcor.qvalue".
#' A user supplied estimate of the true proportion of the null hypothesis. Used only when pi0.known is TRUE. Default is 0.9.
#' @param show.plots logical indicator. When \code{weight.type} is "corr.qvalue" or "parcor.qvalue",
#' \code{show.plots} refers to figures associated with q-value computations as proposed in Storey
#' and Tibshirani (2003). If \code{show.plots} is TRUE, we display the density histogram of original p-values,
#' density histogram of the q-values, scatter plot of \eqn{\hat\pi} versus \eqn{\lambda} in the
#' computation of q-values, and scatter plot of significant tests versus q-value cut-off. When
#' \code{penalty.type} is "penalized.loss", \code{show.plots} refers to plots associated with the
#' penalized loss criterion. If TRUE, a plot of the penalized loss criterion
#' versus steps in the LARS algorithm of Efron et al (2004) is displayed. Default of
#' \code{show.plots} is FALSE.
#' @param qval.alpha scalar value used when \code{weight.type} is "corr.qvalue" or "parcor.qvalue".
#' The choice of \code{qval.alpha} indicates the cut-off for q-values used to obtain the result \code{threshold.selection}
#' The result \code{threshold.selection} contains all covariates for which their q-value is less than \code{qval.alpha}.
#' @param alpha.bh scalar value used when \code{weight.type} is "corr.pvalue","corr.bh.pvalue",
#' "parcor.pvalue", "parcor.bh.pvalue".
#' The choice of \code{alpha.bh} indicates the cut-off for p-values
#' used to obtain the result in \code{threshold.selection}.
#' The result \code{threshold.selection} contains all covariates for which their
#' p-value is less than \code{alpha.bh}.
#' @param penalty.choice character that indicates the variable selection criterion. Options are "cv.mse" for
#' the K-fold cross-validated mean squared prediction error, "penalized.loss" for the penalized loss criterion which
#' requires specification of the penalization parameter \code{penalized.loss.delta},
#' "cv.penalized.loss" for the K-fold cross-validated criterion to determine delta in the penalized loss
#' criterion, and "deviance.criterion" for optimizing the
#' Cox proportional hazards deviance (only available when \code{regression.type} is "cox".) Defalt is "penalized.loss".
#' @param penalized.loss.delta scalar to indicate the choice of the penalization parameter delta in the
#' penalized loss criterion when \code{penalty.choice} is "penalized.loss".
#' @param est.MSE character that indicates how the mean squared error is estimated in the penalized loss
#' criterion when \code{penalty.choice} is "penalized.loss" or "cv.penalized.loss". Options are
#' "est.var" which means the MSE is sd(y) * sqrt(n/(n-1)) where n is the sample size, and
#' "step" which means we use the MSE from forward stepwise regression with AIC as the selection criterion. Default
#' is "est.var".
#' @param cv.folds scalar denoting the number of folds for cross-validation
#' when \code{penalty.choice} is "cv.mse" or "cv.penalized.loss". Default is 10.
#' @param mult.cv.folds scalar denoting the number of times we repeat the cross-validation procedures
#' of \code{penalty.choice} being "cv.mse" or "cv.penalized.loss". Default is 0.
#' @param nboot scalar denoting the number of bootstrap samples obtained for exclusion frequency weights when
#' \code{weight.type} is "exfrequency.random.partition.aic", "exfrequency.random.partition.bic",
#' "exfrequency.kmeans.partition.aic", "exfrequency.kmeans.partition.bic","exfrequency.kquartiles.partition.aic",
#' "exfrequency.kquartiles.partition.bic","exfrequency.ksorted.partition.aic","exfrequency.ksorted.partition.bic".
#' Default is 100.
#' @param k.split scalar that indicates the number of partitions used to compute the exclusion frequency weights
#' when \code{weight.type} is "exfrequency.random.partition.aic", "exfrequency.random.partition.bic",
#' "exfrequency.kmeans.partition.aic", "exfrequency.kmeans.partition.bic","exfrequency.kquartiles.partition.aic",
#' "exfrequency.kquartiles.partition.bic","exfrequency.ksorted.partition.aic","exfrequency.ksorted.partition.bic". Default is 4.
#' @param step.direction character that indicates the direction of stepwise regression used to compute the exclusion frequency weights
#' when \code{weight.type} is "exfrequency.random.partition.aic", "exfrequency.random.partition.bic",
#' "exfrequency.kmeans.partition.aic", "exfrequency.kmeans.partition.bic","exfrequency.kquartiles.partition.aic",
#' "exfrequency.kquartiles.partition.bic","exfrequency.ksorted.partition.aic","exfrequency.ksorted.partition.bic".
#'  One of "both", "forward" or "backward". Default is "backward".
#'
#' @references
#'
#' Efron, B., Hastie, T., Johnstone, I. AND Tibshirani, R. (2004). Least angle regression.
#' Annals of Statistics 32, 407–499.
#'
#' Garcia, T.P. and M¨uller, S. (2016). Cox regression with exclusion frequency-based weights to
#' identify neuroimaging markers relevant to Huntington’s disease onset. Annals of Applied Statistics, 10, 2130-2156.
#'
#' Garcia, T.P. and M¨uller, S. (2014). Influence of measures of significance-based weights in the weighted Lasso.
#' Journal of the Indian Society of Agricultural Statistics (Invited paper), 68, 131-144.
#'
#' Garcia, T.P., Mueller, S., Carroll, R.J., Dunn, T.N., Thomas, A.P., Adams, S.H., Pillai, S.D., and Walzem, R.L.
#' (2013). Structured variable selection with q-values. Biostatistics, DOI:10.1093/biostatistics/kxt012.
#'
#' Storey, J. D. and Tibshirani, R. (2003). Statistical significance for genomewide studies.
#' Proceedings of the National Academy of Sciences 100, 9440-9445.
#'
#' @return
#' \itemize{
#'    \item \strong{weights:} {weights used in the weighted Lasso. Weights computed depend on \code{weight.type} selected.}
#'    \item \strong{weighted.lasso.results:} {variable selection results from the LASSO when the covariates
#'    are multiplied by weights as specified by \code{weight.type}.}
#'    \item \strong{threshold.selection:} {variable selection results when weights are below a specified threshold.
#'    Results are reported only when \code{weight.type} are "corr.pvalue","corr.bh.pvalue",
#'    "corr.qvalue","parcor.pvalue","parcor.bh.pvalue","parcor.qvalue".}
#' }
#'
#' @importFrom stats kmeans
#'
#' @export
#'
#' @examples
#' x = matrix(rnorm(100*5, 0, 1),100,5)
#' z = matrix(rbinom(100, 1, 0.5),100,1)
#' y = matrix(z[,1] + 2*x[,1] - 2*x[,2] + rnorm(100, 0, 1), 100)
#'
#' dwl0 <- d2wlasso(x,z,y)
#' dwl1 <- d2wlasso(x,z=NULL,y,weight.type="corr.pvalue")
#' dwl2 <- d2wlasso(x,z,y,weight.type="parcor.qvalue")
#' dwl3 <- d2wlasso(x,z,y,weight.type="parcor.bh.pvalue")
#' dwl4 <- d2wlasso(x,z,y,weight.type="parcor.qvalue",mult.cv.folds=100)
#' dwl5 <- d2wlasso(x,z,y,weight.type="exfrequency.random.partition.aic")
#'
#' ## Cox model
#' x <- matrix(rnorm(100*5, 0, 1),100,5)
#' z <- matrix(rbinom(100, 1, 0.5),100,1)
#' y <- matrix(exp(z[,1] + 2*x[,1] - 2*x[,2] + rnorm(100, 0, 2)), 100)
#' cox.delta <- matrix(1,nrow=length(y),ncol=1)
#' dwl0.cox <- d2wlasso(x,z,y,cox.delta,regression.type="cox")
#' dwl1.cox <- d2wlasso(x,z=NULL,y,cox.delta,regression.type="cox",weight.type="corr.pvalue")
#' dwl2.cox <- d2wlasso(x,z,y,cox.delta,regression.type="cox",weight.type="parcor.qvalue")
#' dwl3.cox <- d2wlasso(x,z,y,cox.delta,regression.type="cox",weight.type="parcor.bh.pvalue")
#' dwl4.cox <- d2wlasso(x,z,y,cox.delta,regression.type="cox",weight.type="parcor.qvalue",mult.cv.folds=100)
#' dwl5.cox <- d2wlasso(x,z,y,cox.delta,regression.type="cox",weight.type="exfrequency.random.partition.aic")

d2wlasso <- function(x,z,y,cox.delta=NULL,
                     factor.z=TRUE,
                     regression.type=c("linear","cox")[1],
                     weight.type=c("one","corr.estimate","corr.pvalue","corr.bh.pvalue",
                                   "corr.tstat","corr.qvalue",
                                   "parcor.estimate","parcor.pvalue","parcor.bh.pvalue",
                                   "parcor.tstat","parcor.qvalue",
                                   "exfrequency.random.partition.aic",
                                   "exfrequency.random.partition.bic",
                                   "exfrequency.kmeans.partition.aic",
                                   "exfrequency.kmeans.partition.bic",
                                   "exfrequency.kquartiles.partition.aic",
                                   "exfrequency.kquartiles.partition.bic",
                                   "exfrequency.ksorted.partition.aic",
                                   "exfrequency.ksorted.partition.bic")[1],
                     weight_fn=function(x){x},
                     ttest.pvalue=TRUE,
                     q_opt_tuning_method=c("bootstrap","smoother")[2],
                     qval.alpha=0.15,
                     alpha.bh=0.05,
                     robust=TRUE,
                     show.plots=FALSE,
                     pi0.known=FALSE,
                     pi0.val=0.9,
                     penalty.choice=c("cv.mse","cv.penalized.loss","penalized.loss",
                                      "deviance.criterion")[3],
                     est.MSE=c("est.var","step")[1],
                     cv.folds=10,
                     mult.cv.folds=0,
                     penalized.loss.delta=2,
                     nboot=100,
                     k.split=4,
                     step.direction="backward"){

    ######################
    ## WARNING Messages ##
    ######################

    if(pi0.known==TRUE & is.null(pi0.val)){
        stop("User must provide pi0.val if pi0.known is TRUE.")
    }

    if(!is.matrix(x)){
        stop("x must be a n x m matrix.")
    }

    if(!is.matrix(z) | ncol(z)>1){
        stop("z must be a n x 1 matrix.")
    }

    if(!is.matrix(y) | ncol(y)>1){
        stop("y must be a n x 1 matrix.")
    }


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

    z.names <- colnames(z)
    x.names <- colnames(x)
    colnames(XX) <- colnames_use

    ## allocate names to y
    if(is.null(colnames(y))){
        colnames(y) <- paste0("y",1:ncol(y))
    }

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
    ## For some weights, we can apply a threshold to determine if a covariate should be included or not.
    threshold.selection <- NULL
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
        cor.out <- correlations(factor.z,x,z,response,partial=FALSE,ttest.pvalue=ttest.pvalue,
                                regression.type=regression.type)

        #########################
        ## Extract the weights ##
        #########################
        if(weight.type=="corr.tstat"){
            ## t-values with \beta_k in y= \beta_k * x_k
            weights <- cor.out$tvalues
        } else if(weight.type=="corr.pvalue"){
            ## p-values with \beta_k in y= \beta_k * x_k
            weights <- cor.out$pvalues
            benhoch.results <- ben.hoch.interest(weights,alpha=alpha.bh,padjust.method="none")
            threshold.selection <- benhoch.results$interest

            if(!is.null(z)){
                ## we ignore z
                threshold.selection <- rbind(z=0,threshold.selection)
            }

        } else if(weight.type=="corr.estimate"){
            ##  correlation between y and x_k
            weights <- cor.out$estimate
        } else if(weight.type=="corr.bh.pvalue"){
            ## Benjamini-Hochberg adjusted p-values from y=z + \beta_k * x_k
            benhoch.results <- ben.hoch.interest(cor.out$pvalues,alpha=alpha.bh,padjust.method="BH")
            weights <- benhoch.results$pval.adjust
            threshold.selection <- benhoch.results$interest

            if(!is.null(z)){
                ## we ignore z
                threshold.selection <- rbind(z=0,threshold.selection)
            }

        } else if(weight.type=="corr.qvalue"){

            ## Results for testing if a x_k has an effect on y, but NOT
            ##            accounting for z
            ## That is, we test : H_0 : \beta_{x_k}=0

            qvalues.results <- q.computations(cor.out, method=q_opt_tuning_method,
                                              show.plots=show.plots,robust=robust,
                                              pi0.known=pi0.known,pi0.val=pi0.val)

            weights <- qvalues.results$qval.mat
            threshold.selection <- q.interest(weights,alpha=qval.alpha,criteria="less")
            threshold.selection <- threshold.selection$interest

            if(!is.null(z)){
                ## we ignore z
                threshold.selection <- rbind(z=0,threshold.selection)
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
            parcor.out <- correlations(factor.z,x,z,response,partial=TRUE,ttest.pvalue=ttest.pvalue,
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
                benhoch.results <- ben.hoch.interest(weights,alpha=alpha.bh,padjust.method="none")
                threshold.selection <- rbind(z=1,benhoch.results$interest)
            } else if(weight.type=="parcor.estimate"){
                ## partial correlations
                weights <- parcor.out$estimate
            } else if(weight.type=="parcor.bh.pvalue"){
                ## Benjamini-Hochberg adjusted p-values from y=z + \beta_k * x_k
                benhoch.results <- ben.hoch.interest(parcor.out$pvalues,alpha=alpha.bh)
                weights <- benhoch.results$pval.adjust
                threshold.selection <- rbind(z=1,benhoch.results$interest)

            } else if(weight.type=="parcor.qvalue"){
                ## Results for testing if a microbe has an effect on phenotype, but AFTER
                ##            accounting for diet
                ## That is, we test : H_0 : \beta_{x_j|z}=0
                # compute q-value as used by JD Storey with some adjustments made
                qvalues.results <- q.computations(parcor.out,method=q_opt_tuning_method,
                                                             show.plots=show.plots,robust=robust,
                                                             pi0.known=pi0.known,pi0.val=pi0.val)
                weights <- qvalues.results$qval.mat
                threshold.selection <- q.interest(weights,alpha=qval.alpha,criteria="less")
                threshold.selection <- rbind(z=1,threshold.selection$interest)

            }
        } else {
            warning("Your choice of weight.type requires z to be non-NULL.")
        }
    } else if(weight.type=="exfrequency.random.partition.aic" |
              weight.type=="exfrequency.random.partition.bic" |
              weight.type=="exfrequency.kmeans.partition.aic" |
              weight.type=="exfrequency.kmeans.partition.bic" |
              weight.type=="exfrequency.kquartiles.partition.aic"|
              weight.type=="exfrequency.kquartiles.partition.bic"|
              weight.type=="exfrequency.ksorted.partition.aic"|
              weight.type=="exfrequency.ksorted.partition.bic"
              ){

        step.type <- NULL
        if(weight.type=="exfrequency.random.partition.aic" |
           weight.type=="exfrequency.kmeans.partition.aic" |
           weight.type=="exfrequency.kquartiles.partition.aic"|
           weight.type=="exfrequency.ksorted.partition.aic"){
            step.type <- "AIC"
        } else if(weight.type=="exfrequency.random.partition.bic" |
                  weight.type=="exfrequency.kmeans.partition.bic" |
                  weight.type=="exfrequency.kquartiles.partition.bic"|
                  weight.type=="exfrequency.ksorted.partition.bic"){
            step.type <- "BIC"
        }


        ###########################################
        ##  Set up storage for bootstrap results ##
        ###########################################
        tmp2.store <- as.data.frame(matrix(0, nrow = m, ncol = nboot,
                                           dimnames = list(x.names,
                                                           seq(1,nboot))))
        weight.boot <- tmp2.store



        #####################################
        ## Compute intermediary quantities ##
        #####################################
        ## Get parameter estimates from ridge regression
        beta.values <- ridge.regression(XX,response,x.names)

        if(weight.type=="exfrequency.kmeans.partition.aic" |
           weight.type=="exfrequency.kmeans.partition.bic"){
            ## K-means clustering
            kmeans.out <- stats::kmeans(beta.values,centers=k.split,iter.max=100)
            index.group.kmeans <- kmeans.out$cluster
        }

        if(weight.type=="exfrequency.kquartiles.partition.aic" |
           weight.type=="exfrequency.kquartiles.partition.bic"){

            ## K-quart clustering
            index.group.kquart <- cut(beta.values, breaks=quantile(beta.values,
                                        probs=seq(0,1, by=1/k.split)),include.lowest=TRUE)
            index.group.kquart <- as.numeric(factor(index.group.kquart, labels=1:k.split))
        }

        if(weight.type=="exfrequency.ksorted.partition.aic"|
           weight.type=="exfrequency.ksorted.partition.bic"){

            ## K-sort clustering
            beta.sort <- sort(abs(beta.values),decreasing=TRUE,index.return=TRUE)
            sort.beta.index <- beta.sort$ix

            ## index of ordering
            index.group.sort <- index.sort.partition(n=nrow(XX),k=k.split,sort.beta.index,x.names)
        }

        for(b in 1:nboot){
            ##print(b)
            if(weight.type=="exfrequency.random.partition.aic" |
               weight.type=="exfrequency.random.partition.bic"){
                ## Randomly partition the index
                rand.index <- random.partition(n=nrow(XX),p=m,k=k.split)
            }

            if(weight.type=="exfrequency.kmeans.partition.aic" |
               weight.type=="exfrequency.kmeans.partition.bic"){
                ## Partition the index using k-means
                kmeans.rand.index <- designed.partition(index.group.kmeans,k=k.split)
            }

            if(weight.type=="exfrequency.kquartiles.partition.aic" |
               weight.type=="exfrequency.kquartiles.partition.bic"){
                ## Partition the index using k-quartile
                kquart.rand.index <- designed.partition(index.group.kquart,k=k.split)
            }

            if(weight.type=="exfrequency.ksorted.partition.aic"|
               weight.type=="exfrequency.ksorted.partition.bic"){
                ## Partition the index using k-quartile
                sort.rand.index <- designed.partition(index.group.sort,k=k.split)

                ## Measures how often the largest beta from ridge regression is in the true,
                ##  non-zero beta coefficients
                #beta.index.sort[sort.beta.index[1]+1,j] <- as.numeric(sort.beta.index[1]%in%fixed.covariates)
            }

            ## Apply stepwise AIC to each group
            for(l in 1:k.split){
                ##print(l)
                if(weight.type=="exfrequency.random.partition.aic" |
                   weight.type=="exfrequency.random.partition.bic"){
                    ## Random partitioning
                    index <- as.numeric(unlist(rand.index[l]))

                } else if(weight.type=="exfrequency.kmeans.partition.aic" |
                          weight.type=="exfrequency.kmeans.partition.bic"){
                    index <- as.numeric(unlist(kmeans.rand.index[l]))

                } else if(weight.type=="exfrequency.kquartiles.partition.aic" |
                          weight.type=="exfrequency.kquartiles.partition.bic"){
                    index <- as.numeric(unlist(kquart.rand.index[l]))

                } else if(weight.type=="exfrequency.ksorted.partition.aic"|
                   weight.type=="exfrequency.ksorted.partition.bic"){
                    index <- as.numeric(unlist(sort.rand.index[l]))
                }

                if(length(index)!=0){

                    weight.boot[,b] <- weight.boot[,b] +
                                step.selection(factor.z,index,XX,x.names,z.names,
                                               response,type=step.type,
                                               direction=step.direction)
                }
            }
        }

        weights <- apply(weight.boot,1,sum)/nboot
        weights <- t(data.frame(weights))
    }

    ################################
    ## Change rownames of weights ##
    ################################
    rownames(weights) <- x.names

    ############################################################
    ##                                                        ##
    ##                                                        ##
    ##          COMPUTE THE WEIGHTED LASSO RESULTS            ##
    ##                                                        ##
    ##                                                        ##
    ############################################################

    if(mult.cv.folds > 0 & (penalty.choice=="cv.penalized.loss" | penalty.choice="cv.mse")){
        lasso.w <-0

        for(v in 1:mult.cv.folds){
            lasso.w.tmp <- weighted.lasso.computations(weights,weight_fn,response,
                                                       XX,z,z.names,
                                                       show.plots,
                                                       penalty.choice=penalty.choice,
                                                       est.MSE,
                                                       cv.folds,
                                                       penalized.loss.delta)

            lasso.w <- lasso.w + as.matrix(lasso.w.tmp$interest)
        }
        weighted.lasso.results <- lasso.w

    } else {
        lasso.w <- weighted.lasso.computations(weights,weight_fn,response,
                                               XX,z,z.names,
                                               show.plots,
                                               penalty.choice,
                                               est.MSE,
                                               cv.folds,
                                               penalized.loss.delta)

        weighted.lasso.results <- as.matrix(lasso.w$interest)
    }



    return(list(weights=weights,weighted.lasso.results=weighted.lasso.results,
                threshold.selection=threshold.selection))
}





###############################################
## Functions to get a sample weight function ##
###############################################

#' Weight function examples
#'
#' Returns possible weight functions that can be used in d2wlasso.
#'
#' @param weight_fn_type character indicating type of function to return.
#'
#' @return a function that can be applied to a scalar.
#
#' @export
#'
#' @examples
#' weight_fn <- get.weight.fn("sqrt") # return f(x)=sqrt(x)
#'
get.weight.fn <- function(weight_fn_type=c("identity","sqrt","inverse_abs","square")[1]){

    ## Weight functions
    if(weight_fn_type=="identity"){
        out <- function(x){
            return(x)
        }
    } else if(weight_fn_type=="sqrt"){
        out <- function(x){
            return(sqrt(x))
        }
    } else if(weight_fn_type=="inverse_abs"){
        out <- function(x){
            return(1/abs(x))
        }
    } else if(weight_fn_type=="square"){
        out <- function(x){
            return(x^2)
        }
    }
    return(out)
}


#########################################################
## Function for making data frames for storing results ##
#########################################################

#' Matrix to store results for covariates x and vector response y.
#'
#' Forms a matrix that will be used to store regression results of response y on covariates x.
#'
#' @param XX (n by m) matrix of main covariates where m is the number of covariates and n is the sample size.
#' @param yy (n by 1) a matrix corresponding to the response variable.
#'
#' @return (m by 1) matrix of zeroes
#'
#' @export
#'
#' @example
#'
#' xx = matrix(rnorm(100*5, 0, 1),100,5)
#' colnames(xx) <- paste0("X",1:ncol(xx))
#' yy = matrix(2*x[,1] - 2*x[,2] + rnorm(100, 0, 1), 100)
#' colnames(yy) <- "y"
#' store.xy(xx,yy)
store.xy <- function(XX,yy){
	out <- as.data.frame(matrix(0,nrow=ncol(XX),ncol=ncol(yy)))
	rownames(out) <- colnames(XX)
	colnames(out) <- colnames(yy)
	return(out)
}

#' Matrix to store results for covariates x.
#'
#' Forms a matrix that will be used to store regression results of response y on covariates x.
#'
#' @param XX (n by m) matrix of main covariates where m is the number of covariates and n is the sample size.
#'
#' @return (m by 1) matrix of zeroes
#'
#' @export
#'
#' @example
#'
#' xx = matrix(rnorm(100*5, 0, 1),100,5)
#' colnames(xx) <- paste0("X",1:ncol(xx))
#' store.x(xx)
store.x <- function(XX){
	out <- as.data.frame(rep(0,ncol(XX)))
	rownames(out) <- colnames(XX)
	return(out)
}


##########################################################
## Functions to get partial correlation and its p-value #
##########################################################

#' Pearson correlation, p-value and t-statistic associated with the regression between response y and a covariate x
#'
#' @param x (n by 1) matrix corresponding to a covariate vector where n is the sample size.
#' @param y (n by 1) a matrix corresponding to the response variable. If \code{regression.type} is "cox",
#' \code{y} contains the observed event times.
#' @param delta (n by 1) a matrix that denotes censoring when \code{regression.type} is "cox" (1 denotes
#' survival event is observed, 0 denotes the survival event is censored). Can be NULL.
#' @param regression.type a character indicator that is either "linear" for linear regression
#' or "cox" for Cox proportional hazards regression. Default is "linear".
#' @param method character indicating the type of correlation to compute. Default if "pearson"
#' @param alternative character indicating whether the p-value is computed using one-sided or two-sided
#' testing. Default is "two-sided".
#' @param ttest.pvalue logical indicator. If TRUE, p-value for each covariate is computed from univariate
#' linear/cox regression of the response on each covariate. If FALSE, the
#' p-value is computed from correlation coefficients between the response and each covariate.
#' Default is FALSE.
#'
#' @return
#' \itemize{
#'   \item \strong{p.value:}{p-value of the coefficient of x in the regression of y on x.}
#'   \item \strong{estimate:}{Pearson correlation between x and y.}
#'   \item \strong{t.stat:}{t-statistic with testing the significance of x in the regression
#'   of y on x.}
#' }
#'
#'
#'
#' @export
#' @examples
#' x = matrix(rnorm(100*5, 0, 1),100,1)
#' colnames(x) <- paste0("X",1:ncol(xx))
#' y = matrix(2*x[,1] - 2*x[,2] + rnorm(100, 0, 1), 100)
#' colnames(y) <- "y"
#' delta <- NULL
#'
#' output <- corr.pvalue(x,y,delta,regression.type="linear")
#'
#' @importFrom stats lm
#' @importFrom survival coxph Surv
corr.pvalue <- function(x,y,delta,method="pearson",
                        alternative="two.sided",ttest.pvalue=FALSE,regression.type){
    x <- as.numeric(x)
    y <- as.numeric(y)
    delta.use <- as.numeric(delta)

    if(regression.type=="linear"){
        out <- stats::cor.test(x,y,alternative=alternative,method=method,na.action=na.omit)
        estimate <- out$estimate
    } else {
        estimate <- 0
    }

    if(ttest.pvalue==FALSE & regression.type=="linear"){
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


#' Partial correlation, p-value and t-statistic associated with the regression between response y and a covariate x
#' after adjustment for a covariate z.
#'
#' @param x (n by 1) matrix corresponding to a covariate vector where n is the sample size.
#' @param y (n by 1) a matrix corresponding to the response variable. If \code{regression.type} is "cox",
#' \code{y} contains the observed event times.
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable. Can be NULL.
#' @param factor.z logical. If TRUE, the fixed variable z is a factor variable.
#' @param delta (n by 1) a matrix that denotes censoring when \code{regression.type} is "cox" (1 denotes
#' survival event is observed, 0 denotes the survival event is censored). Can be NULL.
#' @param regression.type a character indicator that is either "linear" for linear regression
#' or "cox" for Cox proportional hazards regression. Default is "linear".
#' @param method character indicating the type of correlation to compute. Default if "pearson"
#' @param alternative character indicating whether the p-value is computed using one-sided or two-sided
#' testing. Default is "two-sided".
#' @param ttest.pvalue logical indicator. If TRUE, p-value for each covariate is computed from univariate
#' linear/cox regression of the response on each covariate. If FALSE, the
#' p-value is computed from correlation coefficients between the response and each covariate.
#' Default is FALSE.
#'
#' @return
#' \itemize{
#'   \item \strong{p.value:}{p-value of the coefficient of x in the regression of y on x and z.}
#'   \item \strong{estimate:}{Partial correlation between x and y after adjustment for z.}
#'   \item \strong{t.stat:}{t-statistic with testing the significance of x in the regression
#'   of y on x and z.}
#' }
#'
#'
#'
#' @export
#' @examples
#' x = matrix(rnorm(100*5, 0, 1),100,1)
#' colnames(x) <- paste0("X",1:ncol(xx))
#' z = matrix(rbinom(100, 1, 0.5),100,1)
#' y = matrix(z[,1] + 2*x[,1] - 2*x[,2] + rnorm(100, 0, 1), 100)
#' colnames(y) <- "y"
#' delta <- NULL
#'
#' output <- parcorr.pvalue(factor.z=TRUE,x,y,z,delta,regression.type="linear")
#'
#' @importFrom stats lm
#' @importFrom survival coxph Surv
parcorr.pvalue <- function(factor.z,x,y,z,delta=NULL,method="pearson",alternative="two.sided",ttest.pvalue=FALSE,regression.type){
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
        out <- corr.pvalue(xres,yres,delta.use,method,alternative,ttest.pvalue=FALSE,regression.type=regression.type)
        estimate <- out$estimate
    } else {
        estimate <- 0
    }

    if(ttest.pvalue==FALSE & regression.type=="linear"){
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



#' Correlation, p-value and t-statistic associated with the regression between response y and a covariate x
#' after potential adjustment for a covariate z.
#'
#' @param x (n by m) matrix of main covariates where m is the number of covariates and n is the sample size.
#' @param response a list containing: yy and delta. yy is an (n by 1) matrix corresponding to the response variable. If \code{regression.type} is "cox",
#' \code{y} contains the observed event times. delta is an (n by 1) matrix that denotes censoring when \code{regression.type} is "cox" (1 denotes
#' survival event is observed, 0 denotes the survival event is censored). Can be NULL.
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable. Can be NULL.
#' @param factor.z logical. If TRUE, the fixed variable z is a factor variable.
#' @param regression.type a character indicator that is either "linear" for linear regression
#' or "cox" for Cox proportional hazards regression. Default is "linear".
#' @param partial logical indicator. If TRUE, the partial correlation is computed between the response
#' and each covariate in x after adjustment for z. If FALSE, the correlation is computed between the
#' response and each covariate in x.
#' @param ttest.pvalue logical indicator. If TRUE, p-value for each covariate is computed from univariate
#' linear/cox regression of the response on each covariate. If FALSE, the
#' p-value is computed from correlation coefficients between the response and each covariate.
#' Default is FALSE.
#'
#' @return
#' \itemize{
#'   \item \strong{pvalues:}{p-values for each x_k in the regression of y on x_k when \code{partial}=FALSE. Otherwise,
#'   p-values for each x_k in the regression of y on x_k and z when \code{partial}=TRUE.}
#'   \item \strong{estimate:}{Correlation between x_k and y when \code{partial}=FALSE.
#'   Otherwise, \code{estimate} is the partial correlation between x_k and y after adjustment for z
#'   when \code{partial} is TRUE.}
#'   \item \strong{tvalues:}{t-statistic with testing the significance of x_k in the regression
#'   of y on x_k and z when \code{partial} is TRUE. Otherwise it is the t-statistic when testing the
#'   significance of x_k in the regression of y on x_k when \code{partial} is FALSE.}
#' }
#'
#'
#'
#' @export
#' @examples
#' x = matrix(rnorm(100*5, 0, 1),100,1)
#' colnames(x) <- paste0("X",1:ncol(xx))
#' z = matrix(rbinom(100, 1, 0.5),100,1)
#' y = matrix(z[,1] + 2*x[,1] - 2*x[,2] + rnorm(100, 0, 1), 100)
#' colnames(y) <- "y"
#' delta <- NULL
#' response <- list(yy=y,delta=delta)
#'
#' output <- correlations(factor.z=TRUE,x,z,response, partial=TRUE,regression.type="linear")
#'
correlations <- function(factor.z,x,z,response,partial=FALSE,ttest.pvalue=FALSE,regression.type){

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
				                      z=data.z,delta=data.delta[,j],ttest.pvalue=ttest.pvalue,regression.type=regression.type)
			} else {
				tmp <- corr.pvalue(x=data.XX[,i],y=data.response[,j],delta=data.delta[,j],
				                   ttest.pvalue=ttest.pvalue,regression.type=regression.type)
			}
			correlation[i,j] <- tmp$estimate
			pvalues[i,j] <- tmp$p.value
			tvalues[i,j] <- tmp$t.stat
		}
	}
	list(estimate = correlation, pvalues = pvalues,tvalues=tvalues)
}




##################################################
# Function to get F-test statistics and p-values #
##################################################

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
    smooth.log.pi0 = FALSE,pi0.known=FALSE,pi0.val=0.9)
{
    if (is.null(p)) {
        ##qvalue.gui()
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
      # TG added pi0.known option
      if(pi0.known==TRUE){
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


qvalue.old <- function(p, alpha=NULL, lam=NULL, robust=F,pi0.known=FALSE,pi0.val=0.9)
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
    if(pi0.known==TRUE){
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

#' @import graphics
q.computations <- function(out, method=c("smoother","bootstrap")[2],
				show.plots=TRUE,robust=TRUE,
				pi0.known=FALSE,pi0.val=0.9){

	qval.mat <- matrix(0,nrow=nrow(out$pvalues),ncol=ncol(out$pvalues))
	qval.mat <- as.data.frame(qval.mat)
	rownames(qval.mat) <- rownames(out$pvalues)
	colnames(qval.mat) <- colnames(out$pvalues)


	for(i in 1:ncol(qval.mat)){
		pvalues <- out$pvalues[,i]
		estimate <- out$estimate[,i]

		#qobj <- qvalue.adj(pvalues,pi0.method=method,lambda=seq(0,0.95,by=0.01),robust=FALSE,pi0.known=pi0.known,pi0.val=pi0.val)
		#qval <- qobj$qvalues
		##if(q.old==FALSE){
		    qobj <- qvalue.adj(pvalues,pi0.method=method,lambda=seq(0,0.95,by=0.01),
		                       robust=robust,pi0.known=pi0.known,pi0.val=pi0.val)
		    qval <- qobj$qvalues
		##} else {
		##   qobj <- qvalue.old(pvalues,robust=robust,pi0.known=pi0.known,pi0.val=pi0.val)
		##    qval <- qobj$qvalue
		##}
		pi0 <- qobj$pi0
		qval.mat[,i] <- qval

		cnames <- colnames(out$pvalues)

		# Plots
		if(show.plots==TRUE){

			# Density histogram of p-values
			##postscript(paste(file,"_histpval_",cnames[i],".eps",sep=""))
		    hist(pvalues,main="",xlab="Density of observed p-values",ylab="",freq=FALSE,yaxt="n",ylim=c(0,5),cex.lab=2,cex.axis=2)
			abline(h=1,lty=1)
			abline(h=qobj$pi0,lty=2)
			axis(2,at = c(0,qobj$pi0,1,2,3,4),labels=c(0,round(qobj$pi0,2),1,2,3,4),las=2,cex.lab=2,cex.axis=2)
			box(lty=1,col="black") 	# for a box around plot
			##dev.off()

			# Density histogram of estimate
			##postscript(paste(file,"_histest_",cnames[i],".eps",sep=""))
			hist(estimate,main="",xlab="Density of observed statistic",ylab="",yaxt="n",freq=FALSE,ylim=c(0,5),cex.lab=2,cex.axis=2)
			axis(2,at = c(0,1,2,3,4),labels=c(0,1,2,3,4),las=2,cex.lab=2,cex.axis=2)
			##dev.off()

			# Plot of \hat \pi vs. \lambda
			##postscript(paste(file,"_pihat_vs_lambda_",cnames[i],".eps",sep=""))
			qplot2(qobj,rng=c(0,1),whichplot=1)
			##dev.off()

			# Plot of significant tests vs. q cut-off
			##postscript(paste(file,"_sigtest_vs_qval_",cnames[i],".eps",sep=""))
			qplot2(qobj,rng=c(0,1),whichplot=3)
			##dev.off()
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
#' @importFrom stats p.adjust
ben.hoch.interest <- function(pvalues,alpha=0.05,padjust.method="BH"){

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
    pval.adjust[,i] <- stats::p.adjust(pval,method=padjust.method)

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

#' Compute weighted lasso variable selection
#'
#' Performs variable selection with covariates multiplied by weights that direct which variables
#' are likely to be associated with the response.
#'
#' @param weights (m x 1) matrix that we use to multiply the m-covariates by.
#' @param weight_fn A user-defined function to be applied to the weights for the weighted lasso.
#' Default is an identify function.
#' @param yy (n by 1) a matrix corresponding to the response variable. If \code{regression.type} is "cox",
#' \code{yy} contains the observed event times.
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable. This covariate will
#' always be selected. Can be NULL.
#' @param z.names character denoting the column name of the z-covariate if z is not NULL. Can be NULL.
#' @param XX (n by K) matrix of main covariates where n is the sample size and K=m if z is NULL,
#' and K= m+1 otherwise. Here, m refers to the number of x-covariates.
#' @param data.delta (n by 1) a matrix that denotes censoring when \code{regression.type} is "cox" (1 denotes
#' survival event is observed, 0 denotes the survival event is censored). Can be NULL.
#' @param penalty.choice character that indicates the variable selection criterion. Options are "cv.mse" for
#' the K-fold cross-validated mean squared prediction error, "penalized.loss" for the penalized loss criterion which
#' requires specification of the penalization parameter \code{penalized.loss.delta},
#' "cv.penalized.loss" for the K-fold cross-validated criterion to determine delta in the penalized loss
#' criterion, and "deviance.criterion" for optimizing the
#' Cox proportional hazards deviance (only available when \code{regression.type} is "cox".) Defalt is "penalized.loss".
#' @param show.plots logical indicator. If TRUE and \code{penalty.type} is "penalized.lasso",
#' a plot of the penalized lasso criterion versus steps in the LARS algorithm of Efron et al (2004).
#' @param delta scalar to indicate the choice of the penalization parameter delta in the
#' penalized loss criterion when \code{penalty.choice} is "penalized.loss".
#' @param est.MSE character that indicates how the mean squared error is estimated in the penalized loss
#' criterion when \code{penalty.choice} is "penalized.loss" or "cv.penalized.loss". Options are
#' "est.var" which means the MSE is sd(y) * sqrt(n/(n-1)) where n is the sample size, and
#' "step" which means we use the MSE from forward stepwise regression with AIC as the selection criterion. Default
#' is "est.var".
#' @param cv.folds scalar denoting the number of folds for cross-validation
#' when \code{penalty.choice} is "cv.mse" or "cv.penalized.loss". Default is 10.
#'
#'
#' @references
#' Efron, B., Hastie, T., Johnstone, I. AND Tibshirani, R. (2004). Least angle regression.
#' Annals of Statistics 32, 407–499.
#'
#' Garcia, T.P. and M¨uller, S. (2016). Cox regression with exclusion frequency-based weights to
#' identify neuroimaging markers relevant to Huntington’s disease onset. Annals of Applied Statistics, 10, 2130-2156.
#'
#' Garcia, T.P. and M¨uller, S. (2014). Influence of measures of significance-based weights in the weighted Lasso.
#' Journal of the Indian Society of Agricultural Statistics (Invited paper), 68, 131-144.
#'
#' Garcia, T.P., Mueller, S., Carroll, R.J., Dunn, T.N., Thomas, A.P., Adams, S.H., Pillai, S.D., and Walzem, R.L.
#' (2013). Structured variable selection with q-values. Biostatistics, DOI:10.1093/biostatistics/kxt012.
#'
#' Storey, J. D. and Tibshirani, R. (2003). Statistical significance for genomewide studies.
#' Proceedings of the National Academy of Sciences 100, 9440-9445.
#'
#' @return
#' \itemize{
#'   \item \strong{sig.variables:}{m-dimensional vector where the kth entry is 1 if x_k is selected to be
#'   in the model, and 0 if not.}
#'   \item \strong{sign.of.variables:}{m-dimensional vector where the kth entry is 1 if the coefficient in
#'   front of x_k is positive, and -1 if not.}
#'   \item \strong{delta.out:}{the penalization parameter delta used in the penalized loss criterion.
#'   When \code{penalty.choice} is "penalized.loss", \code{delta.out} is the same as \code{delta}.
#'   When \code{penalty.choice} is "cv.penalized.loss", \code{delta.out} is the resulting delta-value
#'   obtained from the k-fold cross validation.}
#' }
#'
#' @export
#'
#' @import glmnet
#' @import lars
weighted.lasso <- function(weights,weight_fn=function(x){x},yy,XX,z,data.delta,z.names,
                           show.plots=FALSE,
                           penalty.choice,
                           est.MSE=c("est.var","step")[2],
                           cv.folds=10,
                           delta=2
                  ){

    #####################
    ## Set up the data ##
    #####################
    y <- yy
    X <- XX

    N <- length(y)
    y1 <- y
    delta.out <- delta

    # Standardized design matrix X
    X1 <- apply(X,2,make.std)

    ## Transform the weights
    weights.use <- weight_fn(weights)

    # Get rid of small weights
    weights.use <- sapply(weights.use,function(u){max(0.0001,abs(u))})

    ## Be sure to include z
    if(!is.null(z)){
        wts <- c(1e-5,weights.use)
    } else {
        wts <- c(1,weights.use)
    }

    X1 <- X1 %*% diag(1/wts)

    if(is.null(data.delta)){
        #############################
        ## Linear Regression LASSO ##
        #############################
        ## adjust use.Gram
        if(ncol(X1)>500){
            use.Gram <- FALSE
        } else {
            use.Gram <- TRUE
        }


        ## Run Lasso
        wLasso.out <-  lasso.procedure(y1,X1)


        ##order.variables <- as.numeric(wLasso.out$actions)

        if(penalty.choice=="cv.mse"){
            ## Computes the K-fold cross-validated mean squared prediction error for lars, lasso, or forward stagewise.
            ## good implementation, mode="step"
            wLasso.cv <- cv.lars(X1, y1, type = c("lasso"), K = cv.folds,
                                 trace = FALSE, normalize = FALSE, intercept= FALSE,
                                 plot.it=FALSE,mode="step",use.Gram=use.Gram)

            bestindex <- wLasso.cv$index[which.min(wLasso.cv$cv)]

            ## final best descriptive model
            s <- NULL
            predict.out <- predict(wLasso.out, X1,s=bestindex, type = "coefficients", mode="step")
            delta.out <- delta

        } else if(penalty.choice=="cv.penalized.loss"){
            cv.out <- cv.delta(y1,X1,z.names,K=cv.folds,est.MSE=est.MSE,show.plots)
            predict.out <- cv.out$predict.out
            delta.out <- cv.out$delta
        } else if(penalty.choice=="penalized.loss"){
            tmp.out <- penalized.loss.criterion(wLasso.out,y1,X1,z.names,
                                            delta,est.MSE,show.plots)
            predict.out <- tmp.out$predict.out
            delta.out <- delta
        }

    } else {

        ##########################
        ## Cox Regression LASSO ##
        ##########################
        if(is.null(z)){
            penalty <- c(0,rep(1,ncol(X1)-ncol(z)))
        } else {
            penalty <- rep(1,ncol(X1))
        }
        ## Run Lasso
        ytmp <- cbind(time=y1,status=data.delta)
        colnames(ytmp) <- c("time","status")


        if(penalty.choice=="cv.mse"){
            wLasso.cv <- cv.glmnet(X1, ytmp, standardize=FALSE,family="cox",alpha=1,penalty.factor=penalty)
            lambda.opt <- wLasso.cv$lambda.min

        } else if(penalty.choice=="cv.penalized.loss" | "penalized.loss"){
            ## not used
            stop("This method is not implemented for the Cox model.")
        } else if(penalty.choice=="deviance.criterion"){
            wLasso.out <-  glmnet(X1, ytmp, standardize=FALSE,family="cox",alpha=1,penalty.factor=penalty)

            ## BIC/Deviance criterion: deviance + k*log(n)
            deviance <- stats::deviance(wLasso.out)
            p.min <- which.min(deviance + wLasso.out$df*log(N))
            lambda.opt <- wLasso.out$lambda[p.min]
        }

        wLasso.fit <- glmnet(X1, ytmp, standardize=FALSE,family="cox",alpha=1,lambda=lambda.opt)
        ## final best descriptive model
        predict.out <- list(coefficients=wLasso.fit$beta)
    }

    ################################################
    ## Report variables where the coefficient > 0 ##
    ################################################
    ind <- which(as.logical(abs(predict.out$coefficients)>1e-10))

    sig.variables <- rep(0,ncol(XX))
    sig.variables[ind] <- 1

    ##########################################
    ## Report the sign of the coefficients  ##
    ##########################################

    sign.of.variables <- rep(0,ncol(XX))
    ind.pos <- which(as.logical(predict.out$coefficients >0))
    sign.of.variables[ind.pos] <- 1

    ind.neg <- which(as.logical(predict.out$coefficients <0))
    sign.of.variables[ind.neg] <- -1

    list(sig.variables=sig.variables,
         sign.of.variables=sign.of.variables,delta.out=delta.out)
}


#' Wrapper function to store results for  weighted lasso variable selection
#'
#' Performs variable selection with covariates multiplied by weights that direct which variables
#' are likely to be associated with the response.
#'
#' @param weights (m x 1) matrix that we use to multiply the m-covariates by.
#' @param weight_fn A user-defined function to be applied to the weights for the weighted lasso.
#' Default is an identify function.
#' @param response a list containing: yy and delta. yy is an (n by 1) matrix corresponding to the response variable. If \code{regression.type} is "cox",
#' \code{y} contains the observed event times. delta is an (n by 1) matrix that denotes censoring when \code{regression.type} is "cox" (1 denotes
#' survival event is observed, 0 denotes the survival event is censored). Can be NULL.
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable. This covariate will
#' always be selected. Can be NULL.
#' @param z.names character denoting the column name of the z-covariate if z is not NULL. Can be NULL.
#' @param XX (n by K) matrix of main covariates where n is the sample size and K=m if z is NULL,
#' and K= m+1 otherwise. Here, m refers to the number of x-covariates.
#' @param penalty.choice character that indicates the variable selection criterion. Options are "cv.mse" for
#' the K-fold cross-validated mean squared prediction error, "penalized.loss" for the penalized loss criterion which
#' requires specification of the penalization parameter \code{penalized.loss.delta},
#' "cv.penalized.loss" for the K-fold cross-validated criterion to determine delta in the penalized loss
#' criterion, and "deviance.criterion" for optimizing the
#' Cox proportional hazards deviance (only available when \code{regression.type} is "cox".) Defalt is "penalized.loss".
#' @param show.plots logical indicator. If TRUE and \code{penalty.type} is "penalized.loss", a plot of the penalized
#' loss criterion versus steps in the LARS algorithm of Efron et al (2004).
#' @param delta scalar to indicate the choice of the penalization parameter delta in the
#' penalized loss criterion when \code{penalty.choice} is "penalized.loss".
#' @param est.MSE character that indicates how the mean squared error is estimated in the penalized loss
#' criterion when \code{penalty.choice} is "penalized.loss" or "cv.penalized.loss". Options are
#' "est.var" which means the MSE is sd(y) * sqrt(n/(n-1)) where n is the sample size, and
#' "step" which means we use the MSE from forward stepwise regression with AIC as the selection criterion. Default
#' is "est.var".
#' @param cv.folds scalar denoting the number of folds for cross-validation
#' when \code{penalty.choice} is "cv.mse" or "cv.penalized.loss". Default is 10.
#'
#'
#' @references
#' Efron, B., Hastie, T., Johnstone, I. AND Tibshirani, R. (2004). Least angle regression.
#' Annals of Statistics 32, 407–499.
#'
#' Garcia, T.P. and M¨uller, S. (2016). Cox regression with exclusion frequency-based weights to
#' identify neuroimaging markers relevant to Huntington’s disease onset. Annals of Applied Statistics, 10, 2130-2156.
#'
#' Garcia, T.P. and M¨uller, S. (2014). Influence of measures of significance-based weights in the weighted Lasso.
#' Journal of the Indian Society of Agricultural Statistics (Invited paper), 68, 131-144.
#'
#' Garcia, T.P., Mueller, S., Carroll, R.J., Dunn, T.N., Thomas, A.P., Adams, S.H., Pillai, S.D., and Walzem, R.L.
#' (2013). Structured variable selection with q-values. Biostatistics, DOI:10.1093/biostatistics/kxt012.
#'
#' Storey, J. D. and Tibshirani, R. (2003). Statistical significance for genomewide studies.
#' Proceedings of the National Academy of Sciences 100, 9440-9445.
#'
#' @return
#' \itemize{
#'   \item \strong{interest:}{(m by 1) matrix where the kth entry is 1 if x_k is selected to be
#'   in the model, and 0 if not.}
#'   \item \strong{sign.interest:}{(m by 1) matrix where the kth entry is 1 if the coefficient in
#'   front of x_k is positive, and -1 if not.}
#'   \item \strong{delta.out:}{the penalization parameter delta used in the penalized loss criterion.
#'   When \code{penalty.choice} is "penalized.loss", \code{delta.out} is the same as \code{delta}.
#'   When \code{penalty.choice} is "cv.penalized.loss", \code{delta.out} is the resulting delta-value
#'   obtained from the k-fold cross validation.}
#' }
#'
#' @export
#'
weighted.lasso.computations <- function(weights,weight_fn=function(x){x},response,
                                        XX,z,z.names,
                                        show.plots,
                                        penalty.choice,
                                        est.MSE,
                                        cv.folds,
                                        delta){
    #print(response)
    data.response <- response$yy
    data.delta <- response$delta

    interest <- matrix(0,nrow=ncol(XX),ncol=ncol(data.response))
    interest <- as.data.frame(interest)
    rownames(interest) <- colnames(XX)
    colnames(interest) <- colnames(data.response)

    interest.sign <- interest		# matrix to store sign of lasso coefficients


    delta.out <- matrix(0,nrow=1,ncol=ncol(data.response))

    for(i in 1:ncol(interest)){
        #print(data.response)
        lasso.out <-  weighted.lasso(weights[,i],
                       weight_fn,
                       yy=data.response[,i],
                       XX,z,data.delta[,i],
                       z.names,
                       show.plots,
                       penalty.choice,
                       est.MSE,
                       cv.folds,delta)


        interest[,i] <- lasso.out$sig.variables
        interest.sign[,i] <- lasso.out$sign.of.variables

        delta.out[,i] <- lasso.out$delta.out

    }
    list(interest=interest,interest.sign=interest.sign,
         delta.out=delta.out)
}



########################################################################
## Cross-Validation function to select delta in penalized loss criterion ##
########################################################################


#' Cross-validation function to select regularization parameter delta in the penalized loss criterion.
#'
#' Performs variable selection with covariates multiplied by weights that direct which variables
#' are likely to be associated with the response.
#'
#' @param weights (m x 1) matrix that we use to multiply the m-covariates by.
#' @param weight_fn A user-defined function to be applied to the weights for the weighted lasso.
#' Default is an identify function.
#' @param response a list containing: yy and delta. yy is an (n by 1) matrix corresponding to the response variable. If \code{regression.type} is "cox",
#' \code{y} contains the observed event times. delta is an (n by 1) matrix that denotes censoring when \code{regression.type} is "cox" (1 denotes
#' survival event is observed, 0 denotes the survival event is censored). Can be NULL.
#' @param z (n by 1) matrix of additional fixed covariate affecting response variable. This covariate will
#' always be selected. Can be NULL.
#' @param z.names character denoting the column name of the z-covariate if z is not NULL. Can be NULL.
#' @param XX (n by K) matrix of main covariates where n is the sample size and K=m if z is NULL,
#' and K= m+1 otherwise. Here, m refers to the number of x-covariates.
#' @param penalty.choice character that indicates the variable selection criterion. Options are "cv.mse" for
#' the K-fold cross-validated mean squared prediction error, "penalized.loss" for the penalized loss criterion which
#' requires specification of the penalization parameter \code{penalized.loss.delta},
#' "cv.penalized.loss" for the K-fold cross-validated criterion to determine delta in the penalized loss
#' criterion, and "deviance.criterion" for optimizing the
#' Cox proportional hazards deviance (only available when \code{regression.type} is "cox".) Defalt is "penalized.loss".
#' @param show.plots logical indicator. If TRUE and \code{penalty.type} is "penalized.loss", a plot of the penalized
#' loss criterion versus steps in the LARS algorithm of Efron et al (2004).
#' @param delta scalar to indicate the choice of the penalization parameter delta in the
#' penalized loss criterion when \code{penalty.choice} is "penalized.loss".
#' @param est.MSE character that indicates how the mean squared error is estimated in the penalized loss
#' criterion when \code{penalty.choice} is "penalized.loss" or "cv.penalized.loss". Options are
#' "est.var" which means the MSE is sd(y) * sqrt(n/(n-1)) where n is the sample size, and
#' "step" which means we use the MSE from forward stepwise regression with AIC as the selection criterion. Default
#' is "est.var".
#' @param cv.folds scalar denoting the number of folds for cross-validation
#' when \code{penalty.choice} is "cv.mse" or "cv.penalized.loss". Default is 10.
#'
#'
#' @references
#' Efron, B., Hastie, T., Johnstone, I. AND Tibshirani, R. (2004). Least angle regression.
#' Annals of Statistics 32, 407–499.
#'
#' Garcia, T.P. and M¨uller, S. (2016). Cox regression with exclusion frequency-based weights to
#' identify neuroimaging markers relevant to Huntington’s disease onset. Annals of Applied Statistics, 10, 2130-2156.
#'
#' Garcia, T.P. and M¨uller, S. (2014). Influence of measures of significance-based weights in the weighted Lasso.
#' Journal of the Indian Society of Agricultural Statistics (Invited paper), 68, 131-144.
#'
#' Garcia, T.P., Mueller, S., Carroll, R.J., Dunn, T.N., Thomas, A.P., Adams, S.H., Pillai, S.D., and Walzem, R.L.
#' (2013). Structured variable selection with q-values. Biostatistics, DOI:10.1093/biostatistics/kxt012.
#'
#' Storey, J. D. and Tibshirani, R. (2003). Statistical significance for genomewide studies.
#' Proceedings of the National Academy of Sciences 100, 9440-9445.
#'
#' @return
#' \itemize{
#'   \item \strong{interest:}{(m by 1) matrix where the kth entry is 1 if x_k is selected to be
#'   in the model, and 0 if not.}
#'   \item \strong{sign.interest:}{(m by 1) matrix where the kth entry is 1 if the coefficient in
#'   front of x_k is positive, and -1 if not.}
#'   \item \strong{delta.out:}{the penalization parameter delta used in the penalized loss criterion.
#'   When \code{penalty.choice} is "penalized.loss", \code{delta.out} is the same as \code{delta}.
#'   When \code{penalty.choice} is "cv.penalized.loss", \code{delta.out} is the resulting delta-value
#'   obtained from the k-fold cross validation.}
#' }
#'
#' @export

cv.delta <- function(y1,X1,z.names,K=10,est.MSE=c("est.var","step")[1],show.plots){
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
            beta.omit <- penalized.loss.criterion(wLasso.out,y1[-omit],X1[-omit,,drop=FALSE],z.names,
                                            delta=delta.cv[d],est.MSE=est.MSE,
                                            show.plots)

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
    predict.out <- penalized.loss.criterion(wLasso.out,y1,X1,delta=delta,est.MSE=est.MSE,show.plots)$predict.out
    list(predict.out=predict.out,delta=delta)
}

#' @import lars
lasso.procedure <- function(y1,X1){

    ## adjust use.Gram
	if(ncol(X1)>500){
		use.Gram <- FALSE
	} else {
		use.Gram <- TRUE
	}


	# Run Lasso
	wLasso.out <-  lars::lars(X1, y1, type = c("lasso"),
                trace = FALSE, normalize = FALSE, intercept = FALSE,use.Gram=use.Gram)

	list(wLasso.out=wLasso.out)
}

penalized.loss.criterion <- function(wLasso.out,y1,X1,z.names,
                               delta,
                               est.MSE=c("est.var","step")[1],
                               show.plots){
    # Setup
    N <- length(y1)
    p <- ncol(X1)
    s <- length(wLasso.out$df)
    p.pos <- NULL
    RSS = NULL

    for (i in 1:s){
        RSS[i] = sum((y1-predict(wLasso.out, X1, s=i, type = c("fit"))$fit)**2)
        p.pre = predict(wLasso.out, X1, s=i, type = c("coefficients"))$coefficients
        p.pos = c(p.pos,length(p.pre[abs(p.pre)>0]))
    }

    ## Get estimated MSE
    if(est.MSE=="est.var") {
        MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
        MSE <- MSE^2
    } else {
        ## Get estimated MSE from best fit of forward stepwise regression with AIC as selection criterion
        X2 <- data.frame(X1)
        full.lm <- lm(y1~.,data=X2)
        start.fmla <- as.formula(paste("y1~-1+",z.names))
        start.lm <- lm(start.fmla,data=X2)

        lowest.step.forward <- step(lm(start.fmla, data=X2),
                                    list(lower=start.lm,upper=full.lm), direction='forward',trace=FALSE)

        MSE <- summary(lowest.step.forward)$sigma
        MSE <- summary(lowest.step.forward)$sigma^2

        if(MSE < 1e-5){
            MSE <-  sd(as.vector(y1)) * sqrt( N / (N-1) )
            MSE <- MSE^2
        }
    }

    p.min = which.min(RSS/MSE+delta*p.pos)


    ## final best descriptive model
    predict.out <- predict(wLasso.out, X1, s=p.min, type = c("coefficients"))

    if(show.plots==TRUE){
        ##postscript(paste(file,"_lasso1.eps",sep=""))
        plot(wLasso.out,cex.axis=1.5,cex.lab=1.5)
        ##dev.off()
        ##postscript(paste(file,"_lasso2.eps",sep=""))
        ##x11()
        ##plot(s,RSS+2*(s),type="l",cex.axis=1.5,cex.lab=1.5)
        ##abline(v=s.min,lty=2)
        ## Samuel's correction 11/9/2011
        par(mar=c(5, 4, 4, 2)+1)
        plot(1:s,RSS+2*(p.pos),type="l",cex.axis=1.5,cex.lab=1.5,ylab=substitute(M[n](that,p),
                                                                                 list(that=delta)), xlab="Steps")
        abline(v=p.min,lty=2)
        ##dev.off()
    }

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
        stop("k too low")
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

index.sort.partition <- function(n,k,beta.sort.ix,x.names){
    p <- length(beta.sort.ix)

    ## Size of each partition group
    size.groups <- partition.size(n,p,k)

    ## Name of each group
    names <- seq(1,k)
    index.sort <- rep(names, times = size.groups)

    names(index.sort) <- x.names[beta.sort.ix]
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
#' @importFrom glmnet glmnet cv.glmnet
ridge.regression <- function(XX,response,x.names){
    ## Center data
    X <- XX
    X <- apply(X,2,make.center)

    yy <- response$yy
    delta <- response$delta

    if(is.null(delta)){
        #######################
        ## linear regression ##
        #######################
        y <- yy
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

    ## Find optimal lambda value for ridge regression
    cv.out <- glmnet::cv.glmnet(X,y,standardize=FALSE,family=family,alpha=0,grouped=grouped)
    lambda.opt <- cv.out$lambda.min

    ## Ridge regression
    ridge.out <- glmnet::glmnet(X,y,standardize=FALSE,family=family,alpha=0,lambda=lambda.opt)
    beta.values <- ridge.out$beta[x.names,]  ## only take beta values from variables that are subject to selection

    return(beta.values)
}

## We partition based on ascending values of ridge regression estimates
ridge.sort.partition <- function(k,XX,response,x.names){
    ## Get parameter estimates from ridge regression
    beta.values <- ridge.regression(XX,response,x.names)

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
    names_use <- paste0("group",seq(1,k))

    ## Get k-means partition of beta values
    rand.index <- vector("list",length(names_use))
    names(rand.index) <- names_use

    for(j in 1:k){
        cluster.tmp <- as.numeric(which(cluster==j))
        size.groups <- partition.size.new(length(cluster.tmp),k)
        cut.by <- rep(names_use, times = size.groups)
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
#' @import survival
#' @import stats
step.selection <- function(factor.z,index,XX,x.names,z.names,
                           response,type=c("AIC","BIC")[1],
                           direction=c("both","forward","backward")[3]){

    xnam.orig <- x.names[index]

    if(factor.z==TRUE){
        xnam <- c(paste0("factor(",z.names,")"),xnam.orig)
    } else {
        xnam <- c(z.names,xnam.orig)
    }

    yy <- response$yy
    delta <- response$delta

    if(is.null(delta)){
        #######################
        ## linear regression ##
        #######################
        mydata <- data.frame(yy,XX)
        fmla <- as.formula(paste("yy~",paste(xnam,collapse="+")))
        fit <- lm(fmla,data=mydata)
    } else {
        #########################
        ## survival regression ##
        #########################
        X <- data.frame(XX)
        tmp.list <- columns.to.list(X)
        mydata <- list(time=as.numeric(yy),status=as.numeric(delta))
        mydata <- appendList(mydata,tmp.list)
        fmla <- as.formula(paste("Surv(time,status)~",paste(xnam,collapse="+")))
        fit <- survival::coxph(fmla,data=mydata)
    }
    ##print(fmla)


    if(type=="AIC"){
        deg.free <- 2
    } else {
        deg.free <- log(length(yy))
    }

    if(factor.z==TRUE){
        lower.fmla <- as.formula(paste0("~ factor(",z.names,")"))
    } else {
        lower.fmla <- as.formula(paste0("~",z.names))
    }

    step.reg <- stats::step(fit,k=deg.free,direction=direction,
                            scope = list(lower = lower.fmla),trace=FALSE,data=mydata)

    ## XX selected
    results <- intersect(names(step.reg$coefficients),xnam.orig)

    ## Store results
    out <- data.frame(rep(0,length(x.names)))
    rownames(out) <- x.names
    colnames(out) <- "results"
    out[xnam.orig,] <- as.numeric(!xnam.orig%in%results)
    return(out)
}




## Function to do step AIC on group subset
## NOT USED
step.selection.inclusion <- function(factor.z,index,XX,x.names,z.names,
                                     response,type=c("AIC","BIC")[1],
                                     direction=c("both","forward","backward")[3]){

    xnam.orig <- x.names[index]

    if(factor.z==TRUE){
        xnam <- c(paste0("factor(",z.names,")"),xnam.orig)
    } else {
        xnam <- c(z.names,xnam.orig)
    }

    yy <- response$yy
    delta <- response$delta

    if(is.null(delta)){
        #######################
        ## linear regression ##
        #######################
        mydata <- data.frame(yy,XX)
        fmla <- as.formula(paste("yy~",paste(xnam,collapse="+")))
        fit <- lm(fmla,data=mydata)
    } else {
        #########################
        ## survival regression ##
        #########################
        X <- data.frame(XX)
        tmp.list <- columns.to.list(X)
        mydata <- list(time=as.numeric(yy),status=as.numeric(delta))
        mydata <- appendList(mydata,tmp.list)
        fmla <- as.formula(paste("Surv(time,status)~",paste(xnam,collapse="+")))
        fit <- survival::coxph(fmla,data=mydata)
    }
    ##print(fmla)


    if(type=="AIC"){
        deg.free <- 2
    } else {
        deg.free <- log(length(yy))
    }

    if(factor.z==TRUE){
        lower.fmla <- as.formula(paste0("~ factor(",z.names,")"))
    } else {
        lower.fmla <- as.formula(paste0("~",z.names))
    }

    step.reg <- stats::step(fit,k=deg.free,direction=direction,
                            scope = list(lower = lower.fmla),trace=FALSE,data=mydata)

    ## XX selected
    results <- intersect(names(step.reg$coefficients),xnam.orig)

    ## Store results
    out <- data.frame(rep(0,length(x.names)))
    rownames(out) <- x.names
    colnames(out) <- "results"

    out[xnam,] <- as.numeric(xnam%in%results)
    return(out)
}





