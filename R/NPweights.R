#'Distribution function estimation under double truncation
#'
#'This function computes the NPMLE of the event time distribution and truncation time distribution, when the event times
#'are subject to double truncation.
#'@importFrom grDevices dev.new
#'@importFrom graphics lines par persp plot
#'@importFrom stats na.omit
#'@param y vector of event times
#'@param l vector of left truncation times
#'@param r vector of right truncation times
#'@param error prespecified error for convergence (default = 1e-6)
#'@param n.iter maximum number of iterations
#'@param boot Logical. Default=FALSE. If TRUE, the simple bootstrap method is applied to estimate the standard error
#'and pointwise confidence intervals of the event time distribution
#'@param B.boot Numeric value for number of bootstrap resamples. Default is 200.
#'@param joint Logical. Default=FALSE. If TRUE, computes joint and marginal distributions of the truncation times
#'@param plot.cdf Logical. Default is FALSE. If TRUE, the estimated cumulative distribution and survival functions
#'of the event times are plotted. If boot=TRUE, confidence intervals are also plotted.
#'@param plot.joint Logical. Default is FALSE. If TRUE, the estimated marginal distribution functions of the truncation times,
#'and the joint distribution of the truncation times, are plotted. Note: Plot will only be displayed if both plot.joint=TRUE
#'and joint=TRUE.
#'@param display Logical. Default is TRUE. If FALSE, output will not be displayed upon execution of function.
#'@details Estimates the distribution function of the survival time in the presence
#'of left and right truncation. Also estimates the joint cumulative distribution function
#'and marginal cumulative distribution functions of the left and right truncation times.
#'The computation is performed using the algorithm introduced in Shen (2010).
#'This is an iterative algorithm that converges to the NPMLE after a number of iterations.
#'Note that the survival, left, and right truncation times must be the same length.
#'If either of these vectors have missing observations, the entire observation will be excluded.
#'@return \item{time}{Unique event times of the event time vector y}
#'\item{n.event}{Number of events that occurred at each timepoint}
#'\item{F}{Estimated cumulative distribution function of Y at each distinct value of y}
#'\item{Survival}{Estimated survival function of Y at each distinct value of y (equal to 1-F)}
#'\item{sigma.F}{Bootstrapped standard error of F at each distinct value of y (displayed if boot=TRUE)}
#'\item{CI.lower.F}{Estimated lower limits of the Wald confidence intervals for F (displayed if boot=TRUE).}
#'\item{CI.upper.F}{Estimated upper limits of the Wald confidence intervals for F (displayed if boot=TRUE).}
#'\item{P.K}{Probability of the observed vector y falling within a random truncation interval [L,R]}
#'\item{Joint.LR}{Estimated joint distribution function of (l,r)}
#'\item{Marginal.L}{Estimated marginal cumulative distribution function of L at each observed l}
#'\item{Marginal.R}{Estimated marginal cumulative distribution function of R at each observed r}
#'\item{n.iterations}{Number of iterations needed for convergence}
#'\item{max.iter_reached}{0 indicates convergence, 1 indicates that number of iterations exceeded n.iter}
#'@references Shen P-S (2010). Nonparametric analysis of doubly truncated data. Ann Inst Stat Math 62(5):835-853
#'@export
#'@examples
#'#AIDS data set:
#'out=cdfDT(AIDS$Induction.time,AIDS$L.time,AIDS$R.time,plot.cdf=TRUE)
#'out


cdfDT=function(y,l,r,error=1e-6,n.iter=10000,boot=FALSE,B.boot=200,joint=FALSE,plot.cdf=FALSE,plot.joint=FALSE,display=TRUE)
{
  if(joint==FALSE) plot.joint=FALSE;

  # removing missing observations
  temp.data=data.frame(y,l,r);
  nrows.data=dim(temp.data)[1];
  temp.data=na.omit(temp.data);
  nrows.data.omit=dim(temp.data)[1];


  y=temp.data[,1]; u=temp.data[,2]; v=temp.data[,3];
  n=length(y);

  fun.U=function(y,u) I(y>=u)*1;
  fun.V=function(y,v) I(y<=v)*1;

  fun.DT=function(y,u,v)
  {
    n=length(y);
    temp.U=sapply(y,u,FUN="fun.U")
    temp.V=sapply(y,v,FUN="fun.V")

    J=temp.U*temp.V

    K=matrix(0,nrow=n.iter+1,ncol=n); F=matrix(0,nrow=n.iter+1,ncol=n);
    f=matrix(0,nrow=n.iter+1,ncol=n); k=matrix(0,nrow=n.iter+1,ncol=n);

    # S0
    F.0=apply(J,2,"sum")/n



    #S1
    k[1,]=(sum(1/F.0)^-1)/F.0
    K[1,]=apply(k[1,]*J,2,"sum")
    #S2
    f[1,]=(sum(1/K[1,])^-1)/K[1,]
    F[1,]=apply(f[1,]*t(J),2,"sum")

    #S1.2
    k[2,]=(sum(1/F[1,])^-1)/F[1,]
    K[2,]=apply(k[2,]*J,2,"sum")
    #S2.2
    f[2,]=(sum(1/K[2,])^-1)/K[2,]
    F[2,]=apply(f[2,]*t(J),2,"sum")


    s=2;
    while(sum(abs(f[s,]-f[s-1,]))>error)
    {
      s=s+1;

      #S1
      k[s,]=(sum(1/F[s-1,])^-1)/F[s-1,]
      K[s,]=apply(k[s,]*J,2,"sum")

      #S2
      f[s,]=(sum(1/K[s,])^-1)/K[s,]
      F[s,]=apply(f[s,]*t(J),2,"sum")


      if(s>n.iter) break;
    }
    P.K=K[s,]; P.F=F[s,]

    n.unique.y=length(unique(y))
    distF=numeric(n.unique.y)
    for(i in 1:n.unique.y) distF[i]=round(sum(1/P.K)^-1*sum(I(y<=sort(unique(y))[i])/P.K),4);

    f=round(f[s,],4); k=round(k[s,],4);


    max.iter_reached=0; if(s>=n.iter) max.iter_reached=1;
    return(list(f=f,k=k,P.K=P.K,P.F=P.F,distF=distF,n.iterations=s,max.iter_reached=max.iter_reached))
  }

  out=fun.DT(y,u,v);
  P.K=out$P.K; P.F=out$P.F; distF=out$distF; f=out$f; k=out$k; n.iterations=out$n.iterations; max.iter_reached=out$max.iter_reached;
  ###################################################################################
  # NPMLE of truncation distribution (Shen: Semiparametric Analysis of Doubly Truncated Data p. 3183)
  if(joint==TRUE)
  {
    # NPMLE of joint truncation distribution
    unique.u=sort(unique(u)); unique.v=sort(unique(v));
    Joint.UV=matrix(0,nrow=length(unique.u),ncol=length(unique.v));

    for(a in 1:length(unique.u)) {
      for(b in 1:length(unique.v)) {
        Joint.UV[a,b]=(sum(1/(P.F)))^-1*sum(I(u<=unique.u[a])*I(v<=unique.v[b])/(P.F))
      }
    }
    # NPMLE of marginal truncation distributions
    Q.U=Joint.UV[,length(unique.v)]; R.V=Joint.UV[length(unique.u),]
    for(a in 1:length(unique.u)) {
      for(b in 1:length(unique.v)) {
        Joint.UV[a,b]=(sum(1/P.F))^-1*sum(I(u<=unique.u[a])*I(v<=unique.v[b])/P.F)
      }
    }
    Q.U=round(Joint.UV[,length(unique.v)],4); R.V=round(Joint.UV[length(unique.u),],4)
    Joint.UV=round(Joint.UV,4)
  }

  ####################################################################################


  if(boot==TRUE)
  {

    temp.data=data.frame(y,u,v);
    temp.data=temp.data[order(temp.data$y),]
    y.sort=temp.data$y; u.sort=temp.data$u; v.sort=temp.data$v;
    y.unique=sort(unique(y));
    n.unique.y=length(unique(y));

    F.boot=matrix(-1,nrow=B.boot,ncol=n.unique.y)
    for(b in 1:B.boot) {
      repeat{ # creating repeat loop in case program does not converge
        temp.sample=sort(sample(n,replace=TRUE))
        # Creating new data set based off of bootstrapped samples
        y.boot=y.sort[temp.sample]; u.boot=u.sort[temp.sample]; v.boot=v.sort[temp.sample];
        y.boot.unique=(unique(y.boot))


        x1=which(is.element(y.unique,y.boot.unique)==FALSE);
        x2=which(is.element(y.unique,y.boot.unique)==TRUE);
        x3=x1[which((x1>min(x2))&(x1<max(x2)))];
        x3min=x1[which(x1<min(x2))]; x3max=x1[which(x1>max(x2))];

        out.boot=fun.DT(y.boot,u.boot,v.boot)

        if(out.boot$max.iter_reached==0) {break}
      } # ending repeat loop

      F.boot[b,x2]=out.boot$distF

      if(length(x1)>0)
      {
        if(length(x3min)>0) F.boot[b,x3min]=0;
        if(length(x3max)>0) F.boot[b,x3max]=1;
        if(length(x3)>0)
        {
          while(min(F.boot[b,x3])<0)
          {
            F.boot[b,x3]=F.boot[b,x3-1]
          }
        }
      }

    }
    sigma=apply(F.boot,2,"sd")

    CI.lower.F=distF-1.96*sigma; CI.upper.F=distF+1.96*sigma;
  }

  f=f[which(duplicated(y)==FALSE)];
  f=f[order(unique(y))];


  if(display==TRUE)
  {
  if(max.iter_reached==0)
  {
  cat("number of iterations", n.iterations, "\n")
  summary <- cbind(event.time = sort(unique(y)), n.event = table(sort(y)), F = distF, Survival = 1-distF)
  colnames(summary) <- c("time", "n.event",  "cumulative.df", "survival")
  rownames(summary) <- rep("", times = length(unique(y)))
  print(summary, digits = 4, justify = "left")
  cat("number of observations read:", nrows.data, "\n")
  cat("number of observations used:", nrows.data.omit, "\n")
  }
  if(max.iter_reached==1) print("Maximum number of iterations reached. Program did not converge")
  }

  if(plot.cdf==TRUE)
  {
    dev.new()
    par(mfrow=c(1,2))
    plot(distF~sort(unique(y)),ylim=c(0,1),xlab="event time",ylab="",main="Cumulative distribution function")
    lines(distF~sort(unique(y)))
    if(boot==TRUE)
    {
    lines(CI.lower.F~sort(unique(y)),lty=2)
    lines(CI.upper.F~sort(unique(y)),lty=2)
    }
    plot((1-distF)~sort(unique(y)),ylim=c(0,1),xlab="event time",ylab="",main="Survival function")
    lines((1-distF)~sort(unique(y)))
    if(boot==TRUE)
    {
    lines((1-CI.lower.F)~sort(unique(y)),lty=2)
    lines((1-CI.upper.F)~sort(unique(y)),lty=2)
    }
  }

  if(plot.joint==TRUE)
  {
    dev.new()
    par(mfrow=c(1,2))
    plot(Q.U~sort(unique(u)),ylim=c(0,1),xlab="left truncation time",ylab="",main="Marginal cdf (left)")
    lines(Q.U~sort(unique(u)))
    plot(R.V~sort(unique(v)),ylim=c(0,1),xlab="right truncation time",ylab="",main="Marginal cdf (right)")
    lines(R.V~sort(unique(v)))

    dev.new()
    persp(sort(unique(u)),sort(unique(v)),Joint.UV,theta=30,expand=0.75,col="lightblue",main="Joint truncation distribution",
          xlab="left truncation time",ylab="right truncation time",zlab="")
  }


  if(boot==TRUE)
  {

    if(joint==TRUE) return(invisible(list(time=round(sort(unique(y)),4),n.event = table(sort(y)), F = distF, Survival = 1-distF,sigma.F=sigma,CI.lower.F=CI.lower.F,CI.upper.F=CI.upper.F,P.K=P.K,Joint.LR=Joint.UV,Marginal.L=Q.U,Marginal.R=R.V,n.iterations=n.iterations,max.iter_reached=max.iter_reached)));
    if(joint==FALSE) return(invisible(list(time=round(sort(unique(y)),4),n.event = table(sort(y)), F = distF, Survival = 1-distF,F=distF,sigma.F=sigma,CI.lower.F=CI.lower.F,CI.upper.F=CI.upper.F,P.K=P.K,n.iterations=n.iterations,max.iter_reached=max.iter_reached)));
  }

  if(boot==FALSE)
  {
    if(joint==TRUE) return(invisible(list(time=round(sort(unique(y)),4),n.event = table(sort(y)), F = distF, Survival = 1-distF,P.K=P.K,Joint.LR=Joint.UV,Marginal.L=Q.U,Marginal.R=R.V,n.iterations=n.iterations,max.iter_reached=max.iter_reached)));
    if(joint==FALSE) return(invisible(list(time=round(sort(unique(y)),4),n.event = table(sort(y)), F = distF, Survival = 1-distF,P.K=P.K,n.iterations=n.iterations,max.iter_reached=max.iter_reached)));
  }



}

