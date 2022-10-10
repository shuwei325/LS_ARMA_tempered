load("MonteCarloSimulation_wind_ar4_tempered.Rdata")


mean.sd<-function(y){
  tabla<-as.matrix(c(mean(y),sd(y)))
  rownames(tabla)<-c("mean","sd")
  return(tabla)
}

tabla<-apply(results,2,mean.sd)
round(tabla,5)

library(xtable)
xtable(t(tabla),digits=4)
