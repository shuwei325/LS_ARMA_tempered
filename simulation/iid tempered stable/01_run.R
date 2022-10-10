source("00_packages.R")
source("./code/simuTS.R")
source("./code/estimTS.R")

R=1000

no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

progress <- function(n) cat(sprintf("task %d is complete\n", n))
opts <- list(progress=progress)

system.time(
  results<-foreach(i = 1:R, .combine='rbind',
                   .packages=c("MixedTS","stabledist"),.options.snow=opts) %dopar% {
                     
                     #Especificar parametros               
                     N=500 #1000 or 1500
                     
                     Alpha    = 0.9
                     Lambda_p = 0.8
                     Lambda_n = 1.3
                     
                     #Especificar valores iniciais      
                     omega = c(0.8,1,1)
                     
                     #1. simulate
                     if(Alpha<=1){c=0}else{c=5}
                     amostra<-rstdTS.kawai2(m=N,alpha = Alpha,lambda.p = Lambda_p,lambda.n = Lambda_n,c=c,Delta=1)
                     amostra<-amostra$sample
                     
                     #estimates1.1<-stdCTS.estimate.bfgs(initial.sim=omega,y=amostra,n=2^8)
                     estimates1.2<-stdCTS.estimate.nlminb(initial.sim=omega,y=amostra,n=2^8)
                     
                     return(estimates1.2)
                   }    
)

stopCluster(cl)

####################################
####################################

#save(results,file="./results/TSD_simulation_0_1_T1500.Rdata")



