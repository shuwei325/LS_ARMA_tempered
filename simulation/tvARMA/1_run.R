source("00_packages.R")
source("./code/simuTS.R")
source("./code/estimTS.R")
source("./code/tvarma_simu.R")



R=1000

no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

progress <- function(n) cat(sprintf("task %d is complete\n", n))
opts <- list(progress=progress)

system.time(
  results<-foreach(i = 1:R, .combine='rbind',
                   .packages=c("stabledist","LSTS"),.options.snow=opts) %dopar% {
                     
                     N=500 #1000, 1500, 2000 or 3000
                     phia0= 0.3
                     phia1= -0.3
                     phib0= -0.5
                     phib1= 0.4
                     sigma=1.2
                     
                     Alpha=0.7
                     lambda_p=1
                     lambda_n=0.5
                     
                     initial.sim.whittle = c(0.2,0,-0.3,0.5,1.5)
                     initial.sim.TS<-c(0.2,1.2,1.2)
                     
                     if(Alpha<=1){c=0}else{c=4}
                     y_sim=tvARMAsim(N=N,phia0=phia0,phia1=phia1,phib0=phib0,phib1=phib1,sigma=sigma,alpha=Alpha,lambda.p=lambda_p,lambda.n=lambda_n,c=c,Delta=1)
                     
                     
                     TT.=length(y_sim)
                     
                     
                     N. = trunc(TT.^{0.8})
                     S. = trunc(0.2*N.)
                     M. = trunc((TT. - N.)/S. + 1)
                     
                     try(mod.whittle<-LS.whittle(y_sim, start=initial.sim.whittle, order = c(p = 1, q = 1),
                                                 ar.order = c(1), ma.order=c(1) ,sd.order = 0),silent=T)
                     coef.whittle<-c(mod.whittle$coef*c(-1,-1,1,1,1))
                    
                     omega<-c(initial.sim.TS,coef.whittle)
                     
                     TT.<-length(y_sim)
                     t.<-1:TT.
                     u.<-t./TT.
                     
                     phia0<-omega[4]
                     phia1<-omega[5]
                     phib0<-omega[6]
                     phib1<-omega[7]
                     phia<- phia0+phia1*u.
                     phib<- phib0+phib1*u.
                     sigma<-omega[8]
                     
                     e.<-numeric()
                     e.[1]<-y_sim[1]
                     for(i in 2:TT.){
                       e.[i] = y_sim[i] + phia[i]*y_sim[i-1] - phib[i]*e.[i-1]
                     }
                     
                     z<- e./sigma
                     
                     mod.residuals<-stdCTS.estimate.nlminb(initial.sim=initial.sim.TS,y=z,n=2^8)
                     return(c(mod.residuals,coef.whittle))
     }
)

stopCluster(cl)

####################################3
####################################3


#save(results,file="./results/tvARMA_simulation_0_7_T1500.Rdata")

