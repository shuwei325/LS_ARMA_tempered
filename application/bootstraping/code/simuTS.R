library(stabledist)

#asymmetric stable
rS.kawai<-function(n=1,alpha=1.5,a){
  gamma<-(-a*gamma(-alpha)*cos(pi*alpha/2))^(1/alpha)
  sample<-rstable(n=n,alpha=alpha,beta=1,gamma=gamma,delta=0,pm=1)
  return(sample)
}

#one-sided tempered stable
rTS.kawai2.one.side<-function(m=1,alpha=1.5,a=1,lambda=1,c=0.14,Delta=0.01){
  #Approximative sampling with stable proposal distribution
  #In this case, according to Rroji(2015)'s notation: 
  #C_+=a
  #alpha=alpha
  #lambda+=lambda
  sample.ts<-numeric(0)
  acceptance<-numeric(0)
  for(i in 1:m){
    U<-1
    cond<-0   
    k<-0
    while(U>cond){
      U<-runif(n=1,min=0,max=1)
      V<-rS.kawai(n=1,alpha=alpha,a=Delta*a)
      cond<- exp(-lambda*(V+c))
      k=k+1
    }
    sample.ts[i]<- V-Delta*gamma(1-alpha)*a*lambda^(alpha-1)
    acceptance[i]<-k
  }
  return(list(sample=sample.ts,acceptance=acceptance))
}

#tempered stable
rTS.kawai2<-function(m=1,alpha.p=1.5,alpha.n=1.5,C.p=1,C.n=1,lambda.p=1,lambda.n=1,c=0.14,Delta=0.01){
  X.p <- rTS.kawai2.one.side(m=m,alpha=alpha.p,a=C.p,lambda=lambda.p,c=c,Delta=Delta)
  X.n <- rTS.kawai2.one.side(m=m,alpha=alpha.n,a=C.n,lambda=lambda.n,c=c,Delta=Delta)
  X   <- X.p$sample-X.n$sample
  return(list(sample=X,acceptance.p=X.p$acceptance,acceptance.n=X.n$acceptance))
}

####standardized TS distribution Simulation
#tempered stable
rstdTS.kawai2<-function(m=1,alpha = 1.5,lambda.p = 1,lambda.n = 1,c=0.14,Delta=1){
     C = 1/gamma(2-alpha)/(lambda.p^(alpha-2)+lambda.n^(alpha-2))
     amostra<-rTS.kawai2(m=m,alpha.p=alpha,alpha.n=alpha,
                         C.p=C,C.n=C,lambda.p=lambda.p,lambda.n=lambda.n,
                         c=c,Delta=Delta)
    return(amostra)
}


