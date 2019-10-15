##############################################################
## R code for "Structured Variable Selection with q-values" ##
##############################################################

###############
## Libraries ##
###############

library(xtable)		# to create LaTeX tables
library(lars)		# for LASSO approach
library(plotrix)		# for computing standard errors of mean in simulation study
library(survival)   # for survival analysis
library(glmnet)     # for ridge regression


#########################################################
## Function for making data frames for storing results ##
#########################################################

store.micropheno <- function(microbes,phenotypes){
	out <- as.data.frame(matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes)))
	rownames(out) <- rownames(microbes)
	colnames(out) <- rownames(phenotypes)
	return(out)
}

store.micro <- function(microbes){
	out <- as.data.frame(rep(0,nrow(microbes)))
	rownames(out) <- rownames(microbes)
	return(out)
}


##########################################################
## Functions to get partial correlation and its p-value #
##########################################################

# When pearson correlation is 0 (because std. deviation is 0), I am setting
# p-value  to 1.

corr.pvalue <- function(x,y,delta,method="pearson",alternative="two.sided",ttest=FALSE,reg.type){
    ##print(reg.type)
    x <- as.numeric(x)
    y <- as.numeric(y)
    delta.use <- as.numeric(delta)

    if(reg.type=="linear"){
        out <- cor.test(x,y,alternative=alternative,method=method,na.action=na.omit)
        estimate <- out$estimate
    } else {
        estimate <- 0
    }

    if(ttest==FALSE & reg.type=="linear"){
        p.value <- out$p.value
        t.stat=NULL
    } else {
        y1 <- y
        x1 <- x
        delta1 <- delta.use

        if(reg.type=="linear"){
            #######################
            ## linear regression ##
            #######################
            summary.out <- summary(lm(y1~  x1))
            p.value <- summary.out$coefficients["x1","Pr(>|t|)"]
            t.stat <- summary.out$coefficients["x1","t value"]
            ##summary.out <- summary(lm(y~x))
            ##p.value <- summary.out$coefficients["x","Pr(>|t|)"]
        } else {
            #########################
            ## survival regression ##
            #########################
            survival.data <- list(time=y1,status=delta1,x1=x1)
            summary.out <- summary(coxph(Surv(time,status)~ x1,data=survival.data))
            p.value <- summary.out$coefficients["x1","Pr(>|z|)"]
            t.stat  <- summary.out$coefficients["x1","z"]
        }

    }
    list(p.value=p.value,estimate=estimate,t.stat=t.stat)
}


parcorr.pvalue <- function(factor.z,x,y,z,delta=NULL,method="pearson",alternative="two.sided",ttest=FALSE,reg.type){
    x <- as.numeric(x)
    y <- as.numeric(y)
    z <- as.numeric(z)
    delta.use <- as.numeric(delta)

    if(reg.type=="linear"){
        if(factor.z==TRUE){
            xres <- residuals(lm(x~factor(z)))
            yres <- residuals(lm(y~factor(z)))
        } else {
            xres <- residuals(lm(x~z))
            yres <- residuals(lm(y~z))
        }
        out <- corr.pvalue(xres,yres,delta.use,method,alternative,ttest=FALSE,reg.type=reg.type)
        estimate <- out$estimate
    } else {
        estimate <- 0
    }

    if(ttest==FALSE & reg.type=="linear"){
        p.value <- out$p.value
        t.stat <- NULL
    } else {
        y1 <- y
        x1 <- x
        delta1 <- delta.use

        if(reg.type=="linear"){
            #######################
            ## linear regression ##
            #######################

            if(factor.z==TRUE){
                summary.out <- summary(lm(y1 ~ factor(z) +  x1))  ## may have to change due to nature of z!!
            } else {
                summary.out <- summary(lm(y1 ~ z +  x1))  ## may have to change due to nature of z!!
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
                summary.out <- summary(coxph(Surv(time,status)~ factor(z) + x1,data=survival.data))
            } else{
                summary.out <- summary(coxph(Surv(time,status)~ z + x1,data=survival.data))
            }
            p.value <- summary.out$coefficients["x1","Pr(>|z|)"]
            t.stat  <- summary.out$coefficients["x1","z"]
        }
    }


    list(p.value=p.value,estimate=estimate,t.stat=t.stat)
}

correlations <- function(factor.z,microbes,phenotypes,partial=FALSE,ttest=FALSE,format.data=TRUE,reg.type){

    ## Formatting data
    data.response <- phenotypes$yy
    data.delta <- phenotypes$delta
    #if(format.data==TRUE){
    #    data.phenotypes <- data.response[-c(which(rownames(phenotypes)=="Fixed"),which(rownames(phenotypes)=="Cohort")),]
    #} else {
        data.phenotypes <- data.response
    #}
	data.microbes <- microbes[-which(rownames(microbes)=="Fixed"),]
	diet <- microbes["Fixed",]

	# Setting up matrices to store Pearson correlations and p-values
	correlation <- store.micropheno(data.microbes,data.phenotypes)
	pvalues <- store.micropheno(data.microbes,data.phenotypes)
	tvalues <- store.micropheno(data.microbes,data.phenotypes)

	# Computing pearson correlations and p-values
	for (i in 1: nrow(data.microbes)) {
		for(j in 1:nrow(data.phenotypes)) {
			if(partial==TRUE){
				tmp <- parcorr.pvalue(factor.z=factor.z,x=data.microbes[i,],y=data.phenotypes[j,],z=diet,delta=data.delta[j,],ttest=ttest,reg.type=reg.type)
			} else {
				tmp <- corr.pvalue(x=data.microbes[i,],y=data.phenotypes[j,],delta=data.delta[j,],ttest=ttest,reg.type=reg.type)
			}
			correlation[i,j] <- tmp$estimate
			pvalues[i,j] <- tmp$p.value
			tvalues[i,j] <- tmp$t.stat
		}
	}
	list(estimate = correlation, pvalues = pvalues,tvalues=tvalues)
}

# function to compute partial correlations using pcor.R

correlations.pcor <- function(microbes,phenotypes,partial="FALSE",ttest="FALSE",format.data=TRUE){
    ## Formatting data
    if(format.data==TRUE){
        data.phenotypes <- phenotypes[-c(which(rownames(phenotypes)=="Diet"),which(rownames(phenotypes)=="Cohort")),]
    } else {
        data.phenotypes <- phenotypes
    }

    data.microbes <- microbes[-which(rownames(microbes)=="Diet"),]
    diet <- microbes["Diet",]

    ## Setting up matrices to store Pearson correlations and p-values
    correlation <- store.micropheno(data.microbes,data.phenotypes)
    pvalues <- store.micropheno(data.microbes,data.phenotypes)

    ## Computing pearson correlations and p-values
    for (i in 1: nrow(data.microbes)) {
        for(j in 1:nrow(data.phenotypes)) {
            tmp <- pcor.test(as.numeric(data.microbes[i,]),
                             as.numeric(data.phenotypes[j,]),as.numeric(diet))

            correlation[i,j] <- tmp$estimate
            pvalues[i,j] <- tmp$p.value
        }
    }
    list(estimate = correlation, pvalues = pvalues)
}



##################################################
# Function to get F-test statistics and p-values #
##################################################

fstat.pvalue <- function(factor.z,microbe,diet){
    if(factor.z==TRUE){
        diet <- factor(as.numeric(diet))
    } else {
        diet <- as.numeric(diet)
    }
    microbe <- as.numeric(microbe)

    # Assuming variances in both diet groups are the same
    #microbe.lm <- lm(microbe ~ diet,na.action=na.omit)
    #microbe.anova <- anova(microbe.lm)
    #stat <- microbe.anova$F[1]
    #p.value <- microbe.anova$Pr[1]

    # Assuming variances in both diet groups are different
    microbe.ttest <- t.test(microbe ~ diet)
    stat <- microbe.ttest$statistic
    p.value <- microbe.ttest$p.value

    #if(is.na(Fstat)){
    #	Fstat <- NA
    #	p.value <- 1
    #}
    list(Fstat = stat, p.value = p.value)
}

ftests <- function(factor.z,XX){
    # Formatting data
    diet <- XX["Fixed",]
    data.XX <- XX[-which(row.names(XX)=="Fixed"),]

    # Setting up matrices to store F-statistics and p-values
    fstat <- store.micro(data.XX)
    pvalues <- store.micro(data.XX)

    # Computing pearson correlations and p-values
    for (i in 1: nrow(data.XX)) {
        tmp <- fstat.pvalue(factor.z,data.XX[i,],diet)
        fstat[i,] <- tmp$Fstat
        pvalues[i,] <- tmp$p.value
    }

    list(estimate = fstat, pvalues = pvalues)
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
        rng[2]), type = "l", xlab = "q-value cut-off", ylab = "Number of significant microbes",cex.lab=1.5,cex.axis=1.5)
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
		    #print("bye")
		    qobj <- qvalue.adj(pvalues,pi0.method=method,lambda=seq(0,0.95,by=0.01),robust=robust,pi0.true=pi0.true,pi0.val=pi0.val)
		    #qobj <- qvalue.adj(pvalues,pi0.method=method)
		    qval <- qobj$qvalues
		} else {
		    #print("hi there")
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
    pval.adjust[,i] <- p.adjust(pval,method="BH")

    ind <- which(pval.adjust[,i] <= alpha)

    if(length(ind)>0){
      interest[ind,i] <- 1
    }
  }
  list(interest=interest,pval.adjust=pval.adjust)
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
        print(ytmp)
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
    ind <- which(abs(predict.out$coefficients)>1e-10)
    ##ind <- which(predict.out$coefficients!=0)
    sig.variables <- rep(0,nrow(XX))
    sig.variables[ind] <- 1
    sign.of.variables <- rep(0,nrow(XX))
    ind.pos <- which(predict.out$coefficients >0)
    sign.of.variables[ind.pos] <- 1
    ind.neg <- which(predict.out$coefficients <0)
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
                               format.data = TRUE,
                               diet.wt=100,thresh.q=FALSE,corr.g=FALSE,delta=2,std.y="TRUE",
                               est.MSE=c("TRUE","est.var","step")[2],
                               cv.criterion=FALSE,vfold=10){
    #print(response)
    data.response <- response$yy
    data.delta <- response$delta
    #print(data.response)
    #if(format.data==TRUE){
    #    data.response <- data.response[-c(1,2),]
    #} else {
        data.response <- data.response
    #}
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

cv.delta <- function(y1,X1,K=10){
	# sequence of delta values
	delta.cv <- seq(0.75,2,by=0.1)

	# Randomly partition the data
	all.folds <- cv.folds(length(y1), K)

	# Matrix to store residuals
	residmat <- matrix(0, length(delta.cv), K)

	for(j in seq(K)){
		# data set to omit
		omit <- all.folds[[j]]

		# Run Lasso with after removing omit data
		wLasso.out <- lasso.procedure(y1[-omit],X1[-omit,,drop=FALSE])$wLasso.out

		for(d in 1:length(delta.cv)){

			# Find best-fitting model for specified delta
			beta.omit <- lasso.delta.choice(wLasso.out,y1[-omit],X1[-omit,,drop=FALSE],delta=delta.cv[d])

			# Find final fit with data omitted
			fit <- X1[omit,,drop=FALSE] %*% beta.omit$predict.out$coefficients

			# Store residual
			residmat[d,j] <-  apply((y1[omit] - fit)^2, 2, sum)
		}
	}

	cv <- apply(residmat,1,mean)

	# Check which delta's lead to min(cv)
	delta.ind <- which(cv==min(cv))
	delta.opt <- delta.cv[delta.ind]
      delta <- mean(delta.opt)		## takes average of delta values

      wLasso.out <- lasso.procedure(y1,X1)$wLasso.out
	predict.out <- lasso.delta.choice(wLasso.out,y1,X1,delta=delta)$predict.out
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


lasso.delta.choice <- function(wLasso.out,y1,X1,delta){
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

	MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
	MSE <- MSE^2

	p.min = which.min(RSS/MSE+delta*p.pos)


	## final best descriptive model
	predict.out <- predict(wLasso.out, X1, s=p.min, type = c("coefficients"))

	list(wLasso.out=wLasso.out,predict.out=predict.out)
}
