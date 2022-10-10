source("./code/tvar_simu.R")
source("./code/simuTS.R")
source("./code/estimTS.R")
source("./code/ARMA_root_check.R")


library(stabledist)
library(LSTS)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(tcltk)

R=1000

#define cores
no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

progress <- function(n) cat(sprintf("task %d is complete\n", n))
opts <- list(progress=progress)

system.time(
  results<-foreach(i = 1:N, .combine='rbind',
                   .packages=c("stabledist","LSTS"),.options.snow=opts) %dopar% {

  n.iterac<-30
  TT.=1007
  
  Phia0=-1.5985
  Phia1= 0.3305
  Phib0= 0.9135
  Phib1= 0.0207
  Phic0= -0.0585
  Phic1= -0.7153
  Phid0= -0.1316
  Phid1= 0.5454
  Sigma0= 0.0077
  Sigma1= 0.0152
  Alpha=0.5017
  Lambda.p=0.6702
  Lambda.n=0.6798
  
  dadosreais=tvARsim(M=10,N=TT.,phia0=Phia0,phia1=Phia1,phib0=Phib0,phib1=Phib1,
                     phic0=Phic0,phic1=Phic1,phid0=Phid0,phid1=Phid1,
                    alpha = Alpha,lambda.p = Lambda.p,lambda.n = Lambda.n,
                    sigma0=Sigma0,sigma1=Sigma1)
  TT.=length(dadosreais)
  
  
  
  phia0.= 1.5
  phia1.= -0.3
  phib0.= -0.9
  phib1.= -0.04
  phic0.= 0.27
  phic1.= 0.7
  phid0.= 0.10
  phid1.= -0.54
  sigma0.= 0.007 
  sigma1.= 0.003
  
  initial.sim= c(phia0.,phia1.,phib0.,phib1.,phic0.,phic1.,phid0.,phid1.,sigma0.,sigma1.)
  initial.sim.TS=c(0.5,0.7,0.7)
 
  fit<-LS.whittle(dadosreais, start=initial.sim,
                  order = c(p = 4, q = 0),ar.order=c(1,1,1,1),sd.order = c(1))
  
  coef.whittle<-fit$coef*c(-rep(1,8),rep(1,2))
 
  fit$residuals
  
  mod.residuals<-stdCTS.estimate.nlminb(initial.sim=initial.sim.TS,y=fit$residuals,n=2^10)
  
  return(c(mod.residuals,coef.whittle))
  
}
)

stopCluster(cl)

    
#save(results,file="MonteCarloSimulation_wind_ar4_tempered1.R")

