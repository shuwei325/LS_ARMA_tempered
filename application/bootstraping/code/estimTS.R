stdCTS.charact.function<-function(t,alpha=1.5,lambda.p=1,lambda.n=1){
  if( alpha <= 0 | alpha >= 2 ) stop('alpha not in (0, 2)')
  if (alpha ==1 ){
    charact.function<-exp( 1/(lambda.p^(-1)+(lambda.n^(-1)))*
                             (  (lambda.p-1i*t)*log(1-1i*t/lambda.p)+ (lambda.n+1i*t)*log(1+1i*t/lambda.n) ) )
  }else{
    charact.function<-exp( ((lambda.p-t*1i)^alpha  - lambda.p^alpha +(lambda.n+t*1i)^alpha - lambda.n^alpha)/(alpha*(alpha-1)*(lambda.p^(alpha-2)+lambda.n^(alpha-2)))+
                             t*1i*(lambda.p^(alpha-1)-lambda.n^(alpha-1))/((alpha-1)*(lambda.p^(alpha-2)+lambda.n^(alpha-2)))  )
  }
  return(charact.function)
}

characteristic_function_to_density <- function(
  n,   # Number of points, ideally a power of 2
  a, b, # Evaluate the density on [a,b[
  alpha=1.5,
  lambda.p=1,
  lambda.n=1
) {
  i <- 0:(n-1)            # Indices
  dx <- (b-a)/n           # Step size, for the density
  x <- a + i * dx         # Grid, for the density
  dt <- 2*pi / ( n * dx ) # Step size, frequency space
  c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
  d <-  n/2 * dt          # (center the interval on zero)
  t <- c + i * dt         # Grid, frequency space
  phi_t <- stdCTS.charact.function(t,alpha=alpha,lambda.p=lambda.p,lambda.n=lambda.n)
  X <- exp( -(0+1i) * i * dt * a ) * phi_t
  Y <- fft(X)
  density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
  data.frame(
    i = i,
    t = t,
    characteristic_function = phi_t,
    x = x,
    density = Re(density)
  )
}

loglik.stdCTS<-function(omega,amostra,n=2^8){
  alpha = omega[1]
  lambda.p = omega[2]
  lambda.n = omega[3]
  
  a <- min(amostra)-0.1
  b <- max(amostra)+0.1
  
  d <- characteristic_function_to_density(
    n=n,
    a=a, b=b,
    alpha=alpha,
    lambda.p=lambda.p,
    lambda.n=lambda.n
  )
  
  loglik<-0
  for(i in 1:length(amostra)){
    xout = amostra[i]
    prox.x=d$x[abs(d$x-xout) %in% sort(abs(d$x-xout), partial=1:2)[1:2]] #two x closer to the sample
    index<-c(which(d$x==prox.x[1]),which(d$x==prox.x[2]))                #two index closer to the sample
    solution = approx(d$x[index],d$density[index],xout=xout)
    if(solution$y<=0|is.na(solution$y)){
      lik.l<-0.0001
    }else{
      lik.l<-solution$y
    }
    loglik<-log(lik.l)+loglik
  }  
  return(loglik)
}

neg.loglik.stdCTS<-function(omega,amostra,n=2^8){
  alpha = omega[1]
  lambda.p = omega[2]
  lambda.n = omega[3]
  
  a <- min(amostra)-0.1
  b <- max(amostra)+0.1
  
  d <- characteristic_function_to_density(
    n=n,
    a=a, b=b,
    alpha=alpha,
    lambda.p=lambda.p,
    lambda.n=lambda.n
  )
  
  loglik<-0
  for(i in 1:length(amostra)){
    xout = amostra[i]
    prox.x=d$x[abs(d$x-xout) %in% sort(abs(d$x-xout), partial=1:2)[1:2]] #two x closer to the sample
    index<-c(which(d$x==prox.x[1]),which(d$x==prox.x[2]))                #two index closer to the sample
    solution = approx(d$x[index],d$density[index],xout=xout)
    if(solution$y<=0|is.na(solution$y)){
      lik.l<-0.0001
    }else{
      lik.l<-solution$y
    }
    loglik<-log(lik.l)+loglik
  }
  return(-loglik)
}

stdCTS.estimate<-function(initial.sim,y,n=2^8){
  mod<-optim(par=initial.sim,fn=neg.loglik.stdCTS,amostra=y,n=n,lower=c(0.01,0.01,0.01),upper=c(1.99,Inf,Inf))
  return(estimates=mod$par)
}

stdCTS.estimate.bfgs<-function(initial.sim,y,n=2^8){
  mod<-optim(par=initial.sim,fn=neg.loglik.stdCTS,amostra=y,n=n,method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(1.99,Inf,Inf))
  return(estimates=mod$par)
}

stdCTS.estimate.nlminb<-function(initial.sim,y,n=2^8){
  mod<-nlminb(start=initial.sim,neg.loglik.stdCTS,amostra=y,n=n,lower=c(0.01,0.01,0.01),upper=c(1.99,10,10))
  return(estimates=mod$par)
}

