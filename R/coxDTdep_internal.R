#library(survival);

# Necessary Functions
fun.geq=function(a,b) I(a>=b)*1;
fun.leq=function(a,b) I(a<=b)*1;
fun.eq=function(a,b) I(a==b)*1;
fun.length.which=function(a,b) length(which(a==b));

########## EM ALGORITHM ###########

# function to compute baseline hazard
lambda.0=function(beta,y,z,t) return(1/apply(sapply(y,t,FUN="fun.geq")%*%as.matrix(exp(z%*%beta),nrow=length(y)),1,sum))
# function to compute weight matrices (detailed in M step of paper)
fun.W=function(beta,lambda,y.unique,y,z,u,v) {
  n=length(y);
  temp.1=sapply(y.unique,y,FUN="fun.eq") # n x d matrix I(T_i = t_j)
  temp.2=sapply(y.unique,u,FUN="fun.leq")+sapply(y.unique,v,FUN="fun.geq") # computing n x d matrix I(t_j < U_i) + I(t_j > V_i)
  temp.3=exp(z%*%beta)%*%t(lambda);
  temp.4=exp(-exp(z%*%beta)%*%t(cumsum(lambda)));
  temp.5=matrix(rep(exp(-exp(z%*%beta)*sapply(y.unique,u,FUN="fun.leq")%*%lambda),length(y.unique)),nrow=n,ncol=length(y.unique));
  temp.6=matrix(rep(exp(-exp(z%*%beta)*sapply(y.unique,v,FUN="fun.leq")%*%lambda),length(y.unique)),nrow=n,ncol=length(y.unique));
  W=temp.1+temp.2*temp.3*temp.4/(temp.5-temp.6)
  return(W)}

# EM algorithm
#==============================================================
# WAS fun.EM=function(formula,u,v,error,n.iter) {
fun.EM=function(formula,u,v,error,n.iter, data) {
#==============================================================


  #==============================================================
  # WAS mf = model.frame(formula=formula);
  raw.data = data
  mf = model.frame(formula=formula, data=data)
  #==============================================================
  z=model.matrix(attr(mf,"terms"),data=mf)[,-1]
  y=model.response(mf);

  data=data.frame(y,u,v,z); n=dim(data)[1];
  newdata=data[order(y),]; # Ordering data set by survival time
  y=as.numeric(newdata$y)[1:n]; u=newdata$u; v=newdata$v; z=as.matrix(newdata[,4:dim(data)[2]]); y.unique=unique(y)
  d=length(y.unique); num.cov=dim(newdata)[2]-3; # number of unique observations, and covariates

  ### Beginning EM Algorithm
  beta.EM=matrix(0,nrow=n.iter,ncol=num.cov)
  lambda.EM=matrix(0,nrow=n.iter,ncol=d)

  # Step 0 - coeffecients from standard cox model
  #==============================================================
  # WAS beta.EM[1,]=coxph(formula)$coefficients
  beta.EM[1,]=coxph(formula=formula, data=raw.data)$coefficients
  #==============================================================
  lambda.EM[1,]=lambda.0(beta.EM[1,],y,z,y.unique)

  # creating weights
  W=round(fun.W(beta.EM[1,],lambda.EM[1,],y.unique,y,z,u,v),2);
  w=as.vector(t(W));
  w.plus.j=apply(W,2,"sum")
  w.i.plus=apply(W,1,"sum")

  ##### Creating new formula to apply the coxph function with weight vector of length n*d
  y.new=rep(y.unique,length(y))
  status.obs.new=rep(1,length(y.new))
  z.temp=matrix(0,nrow=n*d,ncol=num.cov)
  for(i in 1:num.cov) z.temp[,i]=rep(z[,i],each=d)
  znam <- paste0("z.new", 1:num.cov)
  colnames(z.temp)=znam;

  new.data=data.frame(y.new,status.obs.new,z.temp)
  new.formula=as.formula(paste("Surv(y.new,status.obs.new) ~ ", paste(znam, collapse= "+")))

  ### step 2
  beta.EM[2,]=coxph(new.formula,data=new.data,weights=w,subset=which(w>0))$coefficients

  temp=t(W)%*%exp(z%*%beta.EM[2,]) # column j of W times exp(z*beta)
  lambda.EM[2,]=sapply(1:d, function(j) w.plus.j[j]/sum(temp[j:d]))

  ### Step K:
  k=2;
  while(max(abs(beta.EM[k,]-beta.EM[k-1,]))>error)
  {
    if(k>=n.iter) break;
    k=k+1;
    W=round(fun.W(beta.EM[k-1,],lambda.EM[k-1,],y.unique,y,z,u,v),2)
    w=as.vector(t(W));
    w.plus.j=apply(W,2,"sum")
    w.i.plus=apply(W,1,"sum")
    beta.EM[k,]=coxph(new.formula,data=new.data,weights=w,subset=which(w>0))$coefficients
    temp=t(W)%*%exp(z%*%beta.EM[k,]) # column j of W times exp(z*beta)
    lambda.EM[k,]=sapply(1:d, function(j) w.plus.j[j]/sum(temp[j:d]))
    if(k>n.iter) break;
    #print(k)
  }
  beta.hat.EM=beta.EM[k,]
  lambda.hat.EM=lambda.EM[k,]

  max.iter_reached=0; if(k>=n.iter) max.iter_reached=1;

  if(k<n.iter) return(list(beta.hat=beta.hat.EM,lambda.hat=lambda.hat.EM,n.iterations=k,max.iter_reached=max.iter_reached))
  if(k>=n.iter) return(list(max.iter_reached=max.iter_reached))
}



############ Weighted Estimator: Rennert and Xie (2018)
# Function to compute NP weights
fun.U=function(y,u) I(y>=u)*1;
fun.V=function(y,v) I(y<=v)*1;
fun.temp=function(x) sum(x)^-1
fun.cumulative=function(x) {
  temp=upper.tri(matrix(0,nrow=length(x),ncol=length(x)),diag=TRUE)*1
  temp2=numeric(length(x))
  result=apply(x*temp,2,"sum")
  return(result)
}

weights.NP=function(y,u,v,joint=FALSE,error=1e-6,S=1000)
{
  n=length(y)
  temp.U=sapply(y,u,FUN="fun.U")
  temp.V=sapply(y,v,FUN="fun.V")

  J=temp.U*temp.V

  K=matrix(0,nrow=S+1,ncol=n); F=matrix(0,nrow=S+1,ncol=n);
  f=matrix(0,nrow=S+1,ncol=n); k=matrix(0,nrow=S+1,ncol=n);

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


    if(s>S) break;
  }


  distF=round(fun.cumulative(f[s,]),5)
  P.K=K[s,]
  ###################################################################################
  # NPMLE of K (Shen: Semiparametric Analysis of Doubly Truncated Data p. 3183)
  #x=seq(0,ceiling(max(y)),0.25);
  #G.np=t(sapply(x,u,v,FUN="func"))%*%k[s,];

  if(joint==TRUE)
  {
    # NPMLE of joint CDF of U,V
    max.u.pop=ceiling(max(y)) # or max.u.pop=ceiling(max(u));
    max.v.pop=ceiling(max(v)) # or max.v.pop=ceiling(max(v));
    range.u=seq(0,max.u.pop,0.25); range.v=seq(0,max.v.pop,0.5);
    Joint.UV=matrix(0,nrow=length(range.u),ncol=length(range.v));

    F.v=numeric(length(v)); F.u_m=numeric(length(u));
    f.final=f[s,];

    for(i in 1:n) {
      F.v[i]=sum(f.final*I(y<=v[i]));
      F.u_m[i]=sum(f.final*I(y<u[i]));
    }

    for(a in 1:length(range.u)) {
      for(b in 1:length(range.v)) {
        Joint.UV[a,b]=(sum(1/(F.v-F.u_m)))^-1*sum(I(u<=range.u[a])*I(v<=range.v[b])/(F.v-F.u_m))
      }
    }
    Q.U=Joint.UV[,length(range.v)]; R.V=Joint.UV[length(range.u),]

    # creating separate matrix to compute Joint.UV at the values of the survival times y
    Joint.UV.y=matrix(0,nrow=n,ncol=n+1);
    for(a in 1:n) {
      for(b in 1:n) {
        Joint.UV.y[a,b]=(sum(1/(F.v-F.u_m)))^-1*sum(I(u<=y[a])*I(v<=y[b])/(F.v-F.u_m))
      }
      Joint.UV.y[a,n+1]=(sum(1/(F.v-F.u_m)))^-1*sum(I(u<=y[a])*I(v<=1000)/(F.v-F.u_m))
    }

    Q.U.y=Joint.UV.y[,n+1]; R.V.y=Joint.UV.y[n,1:n]
  }

  ####################################################################################
  f=round(f[s,],5); k=round(k[s,],5);

  max.iter_reached=0;
  if(s<S) {
    if(joint==TRUE) return(list(f=f,F=distF,k=k,P.K=P.K,Joint.UV=Joint.UV,Q.U=Q.U,R.V=R.V,Q.U.y=Q.U.y,R.V.y=R.V.y,n.iterations=s,max.iter_reached=max.iter_reached));
    if(joint==FALSE) return(list(f=f,F=distF,k=k,P.K=P.K,n.iterations=s,max.iter_reached=max.iter_reached));
  }
  if(s>=S) {max.iter_reached=1; return(list(max.iter_reached=max.iter_reached))
  }
}
