
tvARMAsim<-function(N=1000,phia0,phia1,phib0,phib1,sigma=1,alpha=1.5,lambda.p=1,lambda.n=1,c=3,Delta=1){
  #---Generate processes----
  amostra = rstdTS.kawai2(m=(N+1),alpha = alpha,lambda.p = lambda.p,lambda.n = lambda.n,c=c,Delta=1)
  SS = amostra$sample
  tvARMA = matrix(NA, nrow = N+1, ncol = 2)
  tvARMA[,1]=SS
  TT.<-N+1
  t.<-1:TT.
  u.<-t./TT.
  theta1<- phia0+phia1*u.
  theta2<- phib0+phib1*u.
  tvARMA[1,2]<-tvARMA[1,1]
    for (i in 2:(N+1)) {
      tvARMA[i,2] = -theta1[i]*tvARMA[i-1,2] + sigma*(theta2[i]*SS[i-1]+SS[i])
    }
  return(tvARMA=tvARMA[2:(N+1),2])
}


