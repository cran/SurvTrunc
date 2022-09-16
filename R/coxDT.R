#'Fit Cox Proportional Hazards Regression Model Under Independent Double Truncation
#'
#'Fits a Cox proportional hazards regression model when the survival time is subject to both left and right truncation. Assumes that the truncation times are independent of survival times, and no censoring is present in the data.
#'@importFrom survival coxph Surv
#'@importFrom stats model.frame model.matrix model.response na.omit pnorm quantile sd
#'@param formula a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function. NOTE: \code{coxDT} does not handle censoring.
#'@param L vector of left truncation times
#'@param R vector right truncation times
#'@param data mandatory data.frame matrix needed to interpret variables named in the \code{formula}
#'@param subset an optional vector specifying a subset of observations to be used in the fitting process. All observations are included by default.
#'@param time.var default = FALSE. If TRUE, specifies that time varying covariates are included in the model.
#'@param subject a vector of subject identification numbers. Only needed if time.var=TRUE (see example).
#'@param B.SE.np number of iterations for bootstrapped standard error (default = 200)
#'@param CI.boot requests bootstrap confidence intervals (default==FALSE)
#'@param B.CI.boot number of iterations for bootstrapped confidence intervals (default = 2000)
#'@param pvalue.boot requests bootstrap p-values (default==FALSE)
#'@param B.pvalue.boot number of iterations for numerator (estimate) of bootstrapped test statistic (default = 500)
#'@param B.pvalue.se.boot number of iterations for denominator (standard error) of bootstrapped test statistic (default = 100)
#'@param trunc.weight Truncates weights at a prespecified level (default=100). Trade off is a slight increase in bias for reduction in variance.
#'@param print.weights requests the output of nonparametric selection probabilities (default==FALSE)
#'@param error convergence criterion for nonparametric selection probabilities (default = 10e-6)
#'@param n.iter maximum number of iterations for computation of nonparamteric selection probabilities (default = 1000)
#'@details Fits a Cox proportional hazards model in the presence of left and right truncation
#'by weighting each subject in the score equation of the Cox model by the probability that they
#'are observed in the sample. These selection probabilities are computed nonparametrically.
#'The estimation procedure here is performed using coxph {survival} and inserting these estimated
#'selection probabilities in the 'weights' option. This method assumes that the survival and
#'truncation times are independent. Furthermore, this method does not accommodate censoring.
#'Note: If only left truncation is present, set R=infinity.
#'If only right truncation is present, set L = -infinity.
#'
#'@return
#'\item{results.beta}{Displays the estimate, standard error, lower and upper 95\% Wald confidence limits,
#'Wald test statistic and corresponding p-value for each regression coefficient}
#'\item{CI}{Method used for computation of confidence interval: Normal approximation (default) or bootstrap}
#'\item{p.value}{Method used for computation of p-values: Normal approximation (default) or bootstrap}
#'\item{weights}{If print.weights=TRUE, displays the weights used in the Cox model}
#'@references Rennert L and Xie SX (2018). Cox regression model with doubly truncated data. Biometrics, 74(2), 725-733. http://dx.doi.org/10.1111/biom.12809.
#'@export
#'@examples
#'###### Example: AIDS data set #####
#'coxDT(Surv(Induction.time)~Adult,L.time,R.time,data=AIDS,B.SE.np=2)
#'
#'# WARNING: To save computation time, number of bootstrap resamples for standard error set to 2.
#'# Note: The minimum recommendation is 200, which is the default setting.
#'@examples
#'##### Including time-dependent covariates #####
#'# Accomodating time-dependent covariates in the model is similar to the accomodation in coxph
#'
#'# The data set may look like the following:
#'
#'# subject start stop event treatment test.score L.time R.time
#'# 1       T.10  T.11  1    X.1       Z.1        L.1    R.1
#'# 2       T.20  T.21  0    X.2       Z.21       L.2    R.2
#'# 2       T.21  T.22  1    X.2       Z.22       L.2    R.2
#'#...
#'
#'# Here the variable 'treatment' and the trunction times 'L.time' and 'R.time' stay the same
#'# from line to line. The variable 'test.score' will vary line to line. In this example,
#'# subject 1 has only one recorded measurement for test.score, and therefore only has one row
#'# of observations. Subject 2 has two recorded measurements for test score, and therefore has
#'# two rows of observations. In this example, it is assumed that the test score for subject 2
#'# is fixed at Z.21 between (T.20,T.21] and fixed at Z.22 between (T.21,T.22]. Notice that
#'# the event indicator is 0 in the first row of observations corresponding to subject 2,
#'# since they have not yet experienced the event. The status variable changes to 1 in the
#'# row where the event occurs.
#'
#'# Note: Start time cannot preceed left truncation time and must be strictly less than stop time.
#'
#'# example
#'
#'
#'test.data <- data.frame(
#'list(subject.id = c(1, 2, 2, 3, 4, 4, 5, 6, 7, 8, 8, 9, 10),
#'     start      = c(3, 5, 7, 2, 1, 2, 6, 5, 6, 6, 7, 2, 17),
#'     stop       = c(4, 7, 8, 5, 2, 6, 9, 8, 7, 7, 9, 6, 21),
#'     event      = c(1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1),
#'     treatment  = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
#'     test.score = c(5, 6, 7, 4, 6, 9, 3, 4, 4, 7, 6, 4, 12),
#'     L.time     = c(2, 4, 4, 2, 1, 1, 4, 5, 4, 3, 3, 1, 10),
#'     R.time     = c(6, 9, 9, 6, 7, 7, 9, 9, 8, 8, 9, 8, 24)))
#'
#'  coxDT(Surv(start,stop,event)~treatment+test.score,L.time,R.time,data=test.data,
#'  time.var=TRUE,subject=subject.id,B.SE.np=2)




coxDT = function(formula,L,R,data,subset,time.var=FALSE,subject=NULL,B.SE.np=200,CI.boot=FALSE,B.CI.boot=2000,pvalue.boot=FALSE,
                 B.pvalue.boot=500,B.pvalue.se.boot=100,trunc.weight=100,print.weights=FALSE,error=10^-6,n.iter=1000)
{
  if(missing(data)==TRUE) stop("Must specify data fame in 'data =' statement");

  set.seed(1312018)
  data=data[subset,]

  # removing missing observations
  nrows.data=dim(data)[1];
  data=na.omit(data);
  nrows.data.omit=dim(data)[1];


  # extracting outcomes and covariates
  mf = model.frame(formula=formula,data=data)
  X=model.matrix(attr(mf,"terms"),data=mf)[,-1]
  p=1; n=length(X);
  if(length(dim(X))>0) {p=dim(X)[2]; n=dim(X)[1]}  # number of predictors and observations

  Y<-model.response(mf)
  if(length(dim(Y))>0){
    if(any(Y[,2]==0)) stop("Censored observations can not be handled by this method!")
  }

  Y=as.numeric(model.response(mf))[1:n];

  # extracting truncation times
  L=deparse(substitute(L)); R=deparse(substitute(R));
  formula.temp=paste(L,R,sep="~")
  mf.temp=model.frame(formula=formula.temp,data=data)
  obs.data=sapply(rownames(mf),rownames(mf.temp),FUN=function(x,y) which(x==y))
  L=mf.temp[obs.data,1]; R=mf.temp[obs.data,2];



  # weight estimation
  P.obs.y.np=cdfDT(Y,L,R,error,n.iter,display=FALSE)$P.K # estimating selection probabilities for each subject
  weights.np=1/P.obs.y.np                   # weights
  ### truncating weights ###
  weights.np[which(weights.np>trunc.weight)]=trunc.weight
  ### # # # # # # # # #  ###
  P.obs.NP=n*(sum(1/P.obs.y.np))^(-1)       # estimating probability that random subject observed


  # computing estimates of nonparametric weighted estimator
  data.new=data.frame(data,weights.np)
  beta.np=coxph(formula,data=data.new,weights=weights.np)$coefficients

  # computing bootstrapped standard errors

  # first, we import the vector of subject id's for bootstrapping data with time-varying coefficients
  if(time.var==TRUE)
  {
  subject=deparse(substitute(subject))
  formula.temp2=paste(subject,subject,sep="~");
  subjects=model.response(model.frame(formula.temp2,data=data));
  n.subject=length(unique(subjects));
  }

  B=B.SE.np
  if(CI.boot==TRUE) B=max(B.SE.np,B.CI.boot)

  beta.boot.np=matrix(0,nrow=B,ncol=p)
  for(b in 1:B) {
    repeat{ # creating repeat loop in case Shen algorithm fails
      if(time.var==FALSE) {temp.sample=sort(sample(n,replace=TRUE))};
      if(time.var==TRUE) {
        temp.obs=sort(sample(n.subject,replace=TRUE))
        temp.sample=unlist(sapply(temp.obs,subjects,FUN=function(x,y) which(x==y)))
      }
      Y.temp=Y[temp.sample]; L.temp=L[temp.sample]; R.temp=R[temp.sample];
      if(p==1) X.temp=X[temp.sample];
      if(p>=2) X.temp=X[temp.sample,];

      out.temp=cdfDT(Y.temp,L.temp,R.temp,error,n.iter,display=FALSE)

      if(out.temp$max.iter_reached==0) {break}
    } # ending repeat loop

    # non-parametric weights (Shen) for cox regression
    P.obs.NP.temp=out.temp$P.K

    weights.np.temp=1/P.obs.NP.temp
    ### truncating weights ###
    weights.np.temp[which(weights.np.temp>trunc.weight)]=trunc.weight
    ### # # # # # # # # #  ###

    # updating data set to include bootstrapped observations
    data.temp=data.frame(data[temp.sample,],weights.np.temp)

    # computing estimates of nonparametric weighted estimator
    beta.boot.np[b,]=coxph(formula,data=data.temp,weights=weights.np.temp)$coefficients
  }
  # standard error
  se.beta.np=apply(beta.boot.np,2,sd);

  # If bootstrap not requested, return normal confidence intervals
  if(CI.boot==FALSE) {
    CI.lower=beta.np-1.96*se.beta.np; CI.upper=beta.np+1.96*se.beta.np;
    CI.beta.np=round(cbind(CI.lower,CI.upper),3)
  }
  # If bootstrap not requested, return p-values based on normality assumption
  if(pvalue.boot==FALSE)  {
    Test.statistic=(beta.np/se.beta.np)^2;
    p.value=round(2*(1-pnorm(abs(beta.np/se.beta.np))),4)
  }

  # computing 95% confidence intervals
  if(CI.boot==TRUE)
  {
    CI.beta.np=matrix(0,nrow=p,ncol=2)
    for(k in 1:p) CI.beta.np[k,]=round(c(quantile(beta.boot.np[,k],seq(0,1,0.025))[2],quantile(beta.boot.np[,k],seq(0,1,0.025))[40]),3)
  }



  if(pvalue.boot==TRUE)
  {
    B1=B.pvalue.boot;
    B2=B.pvalue.se.boot;
    beta.boot.np1=matrix(0,nrow=B1,ncol=p); beta.boot.np2=matrix(0,nrow=B2,ncol=p); beta.boot.np.sd1=matrix(0,nrow=B1,ncol=p)
    for(b1 in 1:B1) {
      repeat{ # creating repeat loop in case Shen algorithm fails
        if(time.var==FALSE) {temp.sample1=sort(sample(n,replace=TRUE))};
        if(time.var==TRUE) {
          temp.obs1=sort(sample(n.subject,replace=TRUE))
          temp.sample1=unlist(sapply(temp.obs1,subjects,FUN=function(x,y) which(x==y)))
        }
        Y.temp1=Y[temp.sample1]; L.temp1=L[temp.sample1]; R.temp1=R[temp.sample1];
        if(p==1) X.temp1=X[temp.sample1];
        if(p>=2) X.temp1=X[temp.sample1,];

        out.temp1=cdfDT(Y.temp1,L.temp1,R.temp1,error,n.iter,display=FALSE)


        if(out.temp1$max.iter_reached==0) {break}
      } # ending repeat loop

      P.obs.NP.temp1=out.temp1$P.K
      weights.np.temp1=1/P.obs.NP.temp1
      ### truncating weights ###
      weights.np.temp1[which(weights.np.temp1>trunc.weight)]=trunc.weight
      ### # # # # # # # # #  ###

      # updating data set to include bootstrapped observations
      data.temp1=data.frame(data[temp.sample1,],weights.np.temp1)
      # computing estimates of nonparametric weighted estimator
      beta.boot.np1[b1,]=coxph(formula,data=data.temp1,weights=weights.np.temp1)$coefficients



      # The loop below is to compute the standard error of each bootstrap estimate
      for(b2 in 1:B2) {
        repeat{ # creating repeat loop in case Shen algorithm fails
          if(time.var==FALSE) {temp.sample2=sort(sample(temp.sample1,replace=TRUE))};
          if(time.var==TRUE) {
            temp.obs2=sort(sample(temp.obs1,replace=TRUE))
            temp.sample2=unlist(sapply(temp.obs2,subjects,FUN=function(x,y) which(x==y)))
          }
          Y.temp2=Y[temp.sample2]; L.temp2=L[temp.sample2]; R.temp2=R[temp.sample2];
          if(p==1) X.temp2=X[temp.sample2];
          if(p>=2) X.temp2=X[temp.sample2,];

          out.temp2=cdfDT(Y.temp2,L.temp2,R.temp2,error,n.iter,display=FALSE)

          if(out.temp2$max.iter_reached==0) {break}
        } # ending repeat loop

        P.obs.NP.temp2=out.temp2$P.K
        weights.np.temp2=1/P.obs.NP.temp2
        ### truncating weights ###
        weights.np.temp2[which(weights.np.temp2>trunc.weight)]=trunc.weight
        ### # # # # # # # # #  ###


        # updating data set to include bootstrapped observations
        data.temp2=data.frame(data[temp.sample2,],weights.np.temp2)
        # computing estimates of nonparametric weighted estimator
        beta.boot.np2[b2,]=coxph(formula,data=data.temp2,weights=weights.np.temp2)$coefficients
      }
      beta.boot.np.sd1[b1,]=apply(beta.boot.np2,2,"sd")
    }



    # computing test statistics, bootstrapped test statistics, and p-values
    test_statistic.beta.np=numeric(p)
    test_statistic.beta.np.boot=matrix(0,nrow=B1,ncol=p)
    p.value=numeric(p)
    for(k in 1:p) {
      test_statistic.beta.np[k]=(beta.np[k]/se.beta.np[k])^2
      test_statistic.beta.np.boot[,k]=((beta.boot.np1[,k]-beta.np[k])/beta.boot.np.sd1[,k])^2
      p.value[k]=round(length(which(test_statistic.beta.np.boot[,k]>test_statistic.beta.np[k]))/B1,4)
    }
    Test.statistic=test_statistic.beta.np
  }
  # rounding output
  beta.np=round(beta.np,4);
  se.beta.np=round(se.beta.np,4);
  Test.statistic=round(Test.statistic,2)
  p.value[which(p.value==0)]="<.0005"


  results.beta=noquote(cbind(beta.np,se.beta.np,CI.beta.np,Test.statistic,p.value))
  rownames(results.beta)=colnames(X); colnames(results.beta)=c("Estimate","SE","CI.lower","CI.upper","Wald statistic","p-value")

  weights="print option not requested";
  if(print.weights==TRUE) weights=weights.np;

  cat("number of observations read:", nrows.data, "\n")
  cat("number of observations used:", nrows.data.omit, "\n\n")
  if((CI.boot==TRUE)&(pvalue.boot==FALSE)) return(list(results.beta=results.beta,CI="Bootstrap",p.value="Normal approximation",weights=weights));
  if((CI.boot==FALSE)&(pvalue.boot==TRUE)) return(list(results.beta=results.beta,CI="Normal approximation",p.value="Bootstrap",weights=weights));
  if((CI.boot==TRUE)&(pvalue.boot==TRUE))  return(list(results.beta=results.beta,CI="Bootstrap",p.value="Bootstrap",weights=weights));
  if((CI.boot==FALSE)&(pvalue.boot==FALSE))  return(list(results.beta=results.beta,CI="Normal approximation",p.value="Normal approximation",weights=weights));


}

