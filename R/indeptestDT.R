#'Testing quasi-independence between survival and truncation times
#'
#'This function tests for quasi-independence between the survival and truncation times.
#'The survival and truncation times must be quasi-independent to use coxDT and cdfDT.
#'@importFrom stats na.omit pchisq cor
#'@param y vector of event times
#'@param l vector of left truncation times
#'@param r vector of right truncation times
#'
#'@details Testing for quasi-independence between the survival and truncation times using the
#'conditional Kendall's tau introduced by Martin and Betensky (2005). More details are given in their paper.
#'
#'@return
#'\item{tau}{Conditional Kendall's tau for survival time and left truncation time and survival time and right truncation time}
#'\item{X2}{Chi-squared test statistic to test null hypothesis that survival and truncation times are quasi-independent. Default degrees of freedom (DF) is 2. If left and right truncation times perfectly correlated, DF = 1}
#'\item{p}{p-value for null hypothesis that survival and truncation times are quasi-independent}
#'
#'@references Martin and Betensky (2005). Testing Quasi-Independence of Failure and Truncation Times via Conditional Kendall's Tau. JASA. 100(470):484-492.
#'@export
#'@examples
#'# Generating independent survival and truncation times
#' set.seed(123)
#' y=rnorm(30); l=min(y)-abs(rnorm(30)); r=max(y)+abs(rnorm(30))
#'
#'indeptestDT(y,l,r)
#'
#'# Null hypothesis not rejected ==> not enough evidence to reject quasi-independence assumption
indeptestDT=function(y,l,r)
{
  # removing missing observations
  temp.data=data.frame(y,l,r);
  nrows.data=dim(temp.data)[1];
  temp.data=na.omit(temp.data);
  nrows.data.omit=dim(temp.data)[1];

  x=temp.data[,1]; u=temp.data[,2]; v=temp.data[,3];

  n=length(x);
  O=function(a,b,c) I(max(b[1],b[2])<=min(a[1],a[2]))*I(max(a[1],a[2])<=min(c[1],c[2]));
  M=0; # number of comparible pairs
  temp.i=matrix(0,nrow=n,ncol=2); # vector for tau statistic
  temp.i.uv=numeric(n); # vector for tau statistic for truncation times (note this is not used as part of tau-statistic here)
  for(i in 1:(n-1)) {
    temp.j=matrix(0,nrow=(n-i),ncol=2); temp.j.uv=numeric(n-i)
    for(j in (i+1):n) {
      # Tau statistic
      temp.j[(j-i),1]=sign((x[i]-x[j])*(u[i]-u[j]))*O(c(x[i],x[j]),c(u[i],u[j]),c(v[i],v[j]));
      temp.j[(j-i),2]=sign((x[i]-x[j])*(v[i]-v[j]))*O(c(x[i],x[j]),c(u[i],u[j]),c(v[i],v[j]));
      temp.j.uv[(j-i)]=sign((u[i]-u[j])*(v[i]-v[j]))*O(c(x[i],x[j]),c(u[i],u[j]),c(v[i],v[j]));
      if(O(c(x[i],x[j]),c(u[i],u[j]),c(v[i],v[j]))==1) M=M+1;
    }
    # Estimators for tau
    temp.i[i,1]=sum(temp.j[,1]); temp.i[i,2]=sum(temp.j[,2]);
    # estimator for (additional) tau for truncation times
    temp.i.uv[i]=sum(temp.j.uv);
  }

  # Variance estimator
  var.v=numeric(n); var.w=numeric(n); var.vw=numeric(n);
  for(i in 1:n) {
    # temporary vectors for variance estimator
    temp.var.v=numeric(n); temp.var.w=numeric(n); temp.var.gamma.11=numeric(n); temp.var.gamma.12=numeric(n); temp.var.gamma.22=numeric(n);
    for(j in 1:n) {
      # (1.1 terms)
      temp.var.v[j]=sign((x[i]-x[j])*(u[i]-u[j]))*O(c(x[i],x[j]),c(u[i],u[j]),c(v[i],v[j]));
      temp.var.gamma.11[j]=temp.var.v[j]^2

      # (2.2 terms)
      temp.var.w[j]=sign((x[i]-x[j])*(v[i]-v[j]))*O(c(x[i],x[j]),c(u[i],u[j]),c(v[i],v[j]));
      temp.var.gamma.22[j]=temp.var.w[j]^2

      # (1.2 terms)
      temp.var.gamma.12[j]=temp.var.v[j]*temp.var.w[j]
    }


    # estimators for variance
    var.v[i]=sum(temp.var.v[-i])^2-sum(temp.var.gamma.11[-i]);
    var.w[i]=sum(temp.var.w[-i])^2-sum(temp.var.gamma.22[-i]);
    var.vw[i]=sum(temp.var.v[-i])*sum(temp.var.w[-i])-sum(temp.var.gamma.12[-i])


  }
  tau=matrix(apply(temp.i,2,"sum"),ncol=2)/M;
  colnames(tau)=c("tau.L","tau.R")
  tau.uv=sum(temp.i.uv)/M;
  mu=M/choose(n,2);
  U.c=tau*mu;
  #tau.LX=U.c[1,1]/mu; tau.XR=U.c[1,2]/mu;
  var.11=sum(var.v)/(2*n*choose(n-1,2))
  var.12=sum(var.vw)/(2*n*choose(n-1,2))
  var.22=sum(var.w)/(2*n*choose(n-1,2))

  if(cor(u,v)==1){
    X2=0.25*n*U.c[,1]^2/var.11;
    p=1-pchisq(X2,1);
    p=round(1-pchisq(X2,2),3);
    p[which(p==0)]="<.0005";
    cat("Left and right truncation times perfectly correlated","\n");
    cat("DF(chi-squared test statistic) = 1","\n")};

  if(cor(u,v)<1) {
    X2=0.25*n*U.c%*%solve(matrix(c(var.11,var.12,var.12,var.22),nrow=2,ncol=2))%*%t(U.c);
    p=round(1-pchisq(X2,2),3);
    p[which(p==0)]="<.0005";
  cat("DF(chi-squared test statistic) = 2","\n")};


  cat("number of observations read:", nrows.data, "\n")
  cat("number of observations used:", nrows.data.omit, "\n\n")
  return(list(tau=round(tau,3),test.statistic=round(X2,3),p.value=noquote(p)))
}
