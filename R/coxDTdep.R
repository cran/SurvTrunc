#'Fit Cox Proportional Hazards Regression Model Under Dependent Double Truncation
#'
#'Fits a Cox proportional hazards regression model under dependent double truncation. That is, when the survival time is subject to left truncation and/or right truncation and the survival times are dependent on at least one of the truncation times.
#'@importFrom survival coxph Surv
#'@importFrom stats model.frame model.matrix model.response na.omit pnorm quantile sd qnorm as.formula
#'@param formula a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv}. function. NOTE: \code{coxDTdep} does not handle censoring.
#'@param data mandatory data.frame matrix needed to interpret variables named in the \code{formula}
#'@param L vector of left truncation times. If only right truncation is present, set L = -infinity.
#'@param R vector right truncation times. If only left truncation is present, set R=infinity.
#'@param error tolerance for convergence, default is 1e-4.
#'@param n.iter maximun number of iterations for EM algorithm, default is 1000.
#'@param n.boot number of bootstraps for computing standard errors, default is 100.
#'@param CI.level a numeric value between 0.5 and 1 representing the confidence level for two-sided confidence intervals, default is 0.95

#'@details Fits a Cox proportional hazards model in the presence of left, right, or double truncation when the survival times are
#'dependent on at least one of the truncation times. An EM algorithm is employed to obtain point estimates for the regression coefficients.
#'The standard errors are calculated using the bootstrap method. This method assumes no censoring is present in the data.
#'Note: If only left truncation is present, set R=infinity. If only right truncation is present, set L = -infinity.
#'
#'@return
#'Displays the estimate, standard error, lower and upper bounds of confidence interval, Wald test statistic and p-value for each regression coefficient
#'@references Rennert L and Xie SX (2022). Cox regression model under dependent truncation. Biometrics. 78(2), 460-473.  https://doi.org/10.1111/biom.13451
#'@export
#'@examples
#'###### Example: AIDS data set #####
#'\dontrun{coxDTdep(Surv(Induction.time)~Adult, L=AIDS$L.time, R=AIDS$R.time, data=AIDS, n.boot=2)}
#'
#'# WARNING: To save computation time, number of bootstrap resamples for standard error set to 2.
#'# Note: The minimum recommendation is 100, which is the default setting.


coxDTdep<-function(formula, L, R, data,
                   error=.0001, n.iter=1000, n.boot=100, CI.level=0.95){

  if(missing(data)==TRUE) stop("Must specify data fame in 'data =' statement")
  if(CI.level>=1|CI.level<=0.5) stop("CI.level has to be within (0.5, 1)!")

  # extracting outcomes and covariates
  mf<-model.frame(formula=formula, data=data)
  z<-model.matrix(attr(mf,"terms"),data=mf)[,-1]
  z.names<-colnames(model.matrix(attr(mf,"terms"),data=mf))[-1]
  p<-length(z.names)
  y<-model.response(mf)
  if(length(dim(y))>0){
    if(any(y[,2]==0)) stop("Censored observations can not be handled by this method!")
    y<-y[, 1] # no censoring
  }
  y<-as.numeric(y)
  n.obs<-length(y)

  ### Applying EM algorithm ###


  beta.em<-fun.EM(formula=formula,data=data, u=L,v=R,error=error,n.iter=n.iter)$beta.hat

  ### Applying weighted estimator that accounts for double truncation, but not dependence (Rennert and Xie, 2018)
  p.observed.NP<-weights.NP(y=y,u=L,v=R,error=error,S=n.iter)$P.K
  weights.DTdep<-1/p.observed.NP
  data$weights.DTdep<-weights.DTdep
  beta.w.np<-coxph(formula, weights=weights.DTdep, data=data)$coefficients




  ########## BOOTSTRAP STANDARD ERRORS #################

  # n.boot: number of bootstrap replications: this value can be lowered to speed up processing, but analysis based on B < 100 is not recommended
  set.seed(10042020)

  beta.boot.em<-numeric(n.boot)
  beta.boot.w.np<-numeric(n.boot) # change two covariates to one
  system.time(
    for(b in 1:n.boot) {
      repeat{ # creating repeat loop in case DTDA program crashes
        boot.sample<-sort(sample(n.obs,replace=TRUE))
        data.boot<-data[boot.sample, ]
        Y.boot<-y[boot.sample]
        U.boot<-L[boot.sample]
        V.boot<-R[boot.sample]  # Ignoring A.boot<-A[boot.sample]

        ### Applying EM algorithm ###
        out.EM.boot<-fun.EM(formula,u=U.boot,v=V.boot, data=data.boot, error=error,n.iter=n.iter)

        if(out.EM.boot$max.iter_reached==0)  {break}
      } # ending repeat loop

      ### em algorithm ###
      beta.boot.em[b]<-out.EM.boot$beta.hat


    }
  )


  # standard errors
  se.beta.em<-sd(beta.boot.em)

  ### confidence intervals based on bootstrapped standard errors ###
  CI.q<--qnorm((1-CI.level)/2)
  # EM
  CI.lower.em<-beta.em-CI.q*se.beta.em
  CI.upper.em<-beta.em+CI.q*se.beta.em
  Test.statistic<-(beta.em/se.beta.em)^2;
  p.value<-2*(1-pnorm(abs(beta.em/se.beta.em)))


  ############### RESULTS #######################
  table.Results<-c(beta.em,se.beta.em, CI.lower.em, CI.upper.em, Test.statistic, p.value)
  table.Results<-matrix(table.Results, nrow=p)
  rownames(table.Results)<-z.names
  colnames(table.Results)<-c("Estimate","SE","CI.lower","CI.upper","Wald statistic","p-value")

  return(table.Results)
}
