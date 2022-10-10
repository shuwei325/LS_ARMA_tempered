library(stabledist)
tvARMAsim<-function(M,N,phia0,phia1,phib0,phib1,alpha=2,beta=0,gamma=1,delta=0,sigma0=1,sigma1=0.1){
  #---Generate processes----
  SS = matrix(rstable(M*(N+1),alpha,beta,gamma,delta,pm=1),ncol=M)
  tvARMA = matrix(NA, nrow = N+1, ncol = M)
  tvARMA[1,]=SS[1,]
  TT.<-N+1
  t.<-1:TT.
  u.<-t./TT.
  theta1<- phia0+phia1*u.
  theta2<- phib0+phib1*u.
  sigma<-sigma0+sigma1*u.
  
  for (j in 1:M){
    for (i in 2:(N+1)) {
      tvARMA[i,j] = -theta1[i]*tvARMA[i-1,j] + sigma[i]*(theta2[i]*SS[i-1,j]+SS[i,j])
    }
  }
  return(tvARMA=tvARMA[2:(N+1),])
}


