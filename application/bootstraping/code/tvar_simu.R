
tvARsim<-function(M=10,N=500,phia0,phia1,phib0,phib1,phic0,phic1,phid0,phid1,
                  alpha = 0.5,lambda.p = 1,lambda.n = 1,
                  sigma0=1,sigma1=0.1){
  #---Generate processes----
  
  if(alpha<=1){c=0}else{c=4}
  SS = rstdTS.kawai2(m=(N+M),alpha = alpha,lambda.p = lambda.p,lambda.n = lambda.n,c=c,Delta=1)$sample
  LSAR = rep(NA, length = N+M)
  
  TT.<-N
  t.<-1:TT.
  u.<-t./TT.
  phia.<- phia0+phia1*u.
  phib.<- phib0+phib1*u.
  phic.<- phic0+phic1*u.
  phid.<- phid0+phid1*u.
  
  sigma.<-sigma0+sigma1*u.
  
  phia<-c(rep(phia0,M),phia.)
  phib<-c(rep(phib0,M),phib.)
  phic<-c(rep(phic0,M),phic.)
  phid<-c(rep(phid0,M),phid.)
  sigma<-c(rep(sigma0,M),sigma.)
  
  LSAR[1:4]=sigma[1:4]*SS[1:4]
  
  for (i in 5:length(LSAR)) {
      LSAR[i] = -phia[i]*LSAR[i-1]-phib[i]*LSAR[i-2]-phic[i]*LSAR[i-3]-phid[i]*LSAR[i-4]+sigma[i]*SS[i]
  }
  return(LSAR=LSAR[(M+1):length(LSAR)])
}


