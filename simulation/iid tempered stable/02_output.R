library(car)
library(moments)
library(xtable)
library(latex2exp)


load("./results/TSD_simulation_0_4_T500.Rdata")
resultsA2T1<-results
load("./results/TSD_simulation_0_4_T1000.Rdata")
resultsA2T2<-results
load("./results/TSD_simulation_0_4_T1500.Rdata")
resultsA2T3<-results

row.names(resultsA2T1)<-NULL
row.names(resultsA2T2)<-NULL
row.names(resultsA2T3)<-NULL

###############Summary
mean.sd<-function(y){
  tabla<-as.matrix(c(mean(y),sd(y)))
  rownames(tabla)<-c("mean","sd")
  return(tabla)
}

kurtosis.skewness<-function(y){
  tabla<-as.matrix(c(kurtosis(y),skewness(y)))
  return(tabla)
}

tablaA2T1<-apply(resultsA2T1,2,mean.sd)
tablaA2T2<-apply(resultsA2T2,2,mean.sd)
tablaA2T3<-apply(resultsA2T3,2,mean.sd)

parametrosA2<-c(0.4,0.5,1.5)
parametrosA2

tablaA2.descriptiva<-round(rbind(tablaA2T1,tablaA2T2,tablaA2T3),4)
tablaA2.descriptiva


tabla.stdCTS<-round(rbind(tablaA2.descriptiva),4)
colnames(tabla.stdCTS)<-c("alpha","lambda_p","lambda_n")
tabla.stdCTS


tablaA2T1.m<-apply(resultsA2T1,2,kurtosis.skewness)
tablaA2T2.m<-apply(resultsA2T2,2,kurtosis.skewness)
tablaA2T3.m<-apply(resultsA2T3,2,kurtosis.skewness)
rownames(tablaA2T1.m)<-c("kurtosis","skewness")
rownames(tablaA2T2.m)<-c("kurtosis","skewness")
rownames(tablaA2T3.m)<-c("kurtosis","skewness")

momments.A2<-rbind(tablaA2T1,tablaA2T1.m,tablaA2T2,tablaA2T2.m,tablaA2T3,tablaA2T3.m)

xtable(momments.A2,digits=4)


##############################fig


pdf("./output/stdCTS_alpha0_4.pdf")

par(mfrow=c(1,1),cex=1, lwd=1)

qqPlot(resultsA2T1[,1],ylab=TeX("$\\alpha$"))
qqPlot(resultsA2T1[,2],ylab=TeX("$\\lambda_+$"))
qqPlot(resultsA2T1[,3],ylab=TeX("$\\lambda_-$"))

qqPlot(resultsA2T2[,1],ylab=TeX("$\\alpha$"))
qqPlot(resultsA2T2[,2],ylab=TeX("$\\lambda_+$"))
qqPlot(resultsA2T2[,3],ylab=TeX("$\\lambda_-$"))

qqPlot(resultsA2T3[,1],ylab=TeX("$\\alpha$"))
qqPlot(resultsA2T3[,2],ylab=TeX("$\\lambda_+$"))
qqPlot(resultsA2T3[,3],ylab=TeX("$\\lambda_-$"))


par(mfrow=c(1,1),cex=1.5, lwd=5)
hist(resultsA2T1[,1],freq=FALSE,xlab=TeX("$\\alpha$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[1],col=3)
curve(dnorm(x,mean(resultsA2T1[,1]),sd(resultsA2T1[,1])),add=T,col=2)
hist(resultsA2T1[,2],freq=FALSE,xlab=TeX("$\\lambda_+$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[2],col=3)
curve(dnorm(x,mean(resultsA2T1[,2]),sd(resultsA2T1[,2])),add=T,col=2)
hist(resultsA2T1[,3],freq=FALSE,xlab=TeX("$\\lambda_-$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[3],col=3)
curve(dnorm(x,mean(resultsA2T1[,3]),sd(resultsA2T1[,3])),add=T,col=2)

hist(resultsA2T2[,1],freq=FALSE,xlab=TeX("$\\alpha$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[1],col=3)
curve(dnorm(x,mean(resultsA2T2[,1]),sd(resultsA2T2[,1])),add=T,col=2)
hist(resultsA2T2[,2],freq=FALSE,xlab=TeX("$\\lambda_+$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[2],col=3)
curve(dnorm(x,mean(resultsA2T2[,2]),sd(resultsA2T2[,2])),add=T,col=2)
hist(resultsA2T2[,3],freq=FALSE,xlab=TeX("$\\lambda_-$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[3],col=3)
curve(dnorm(x,mean(resultsA2T2[,3]),sd(resultsA2T2[,3])),add=T,col=2)


hist(resultsA2T3[,1],freq=FALSE,xlab=TeX("$\\alpha$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[1],col=3)
curve(dnorm(x,mean(resultsA2T3[,1]),sd(resultsA2T3[,1])),add=T,col=2)
hist(resultsA2T3[,2],freq=FALSE,xlab=TeX("$\\lambda_+$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[2],col=3)
curve(dnorm(x,mean(resultsA2T3[,2]),sd(resultsA2T3[,2])),add=T,col=2)
hist(resultsA2T3[,3],freq=FALSE,xlab=TeX("$\\lambda_-$"),main="", cex.lab = 1.5)
abline(v=parametrosA2[3],col=3)
curve(dnorm(x,mean(resultsA2T3[,3]),sd(resultsA2T3[,3])),add=T,col=2)


plot(density(resultsA2T1[,1],adjust = 1.5),lty=1,col=2,ylim=c(0,5),xlab=TeX("$\\alpha$"),main="", cex.lab = 1.5)
points(density(resultsA2T2[,1],adjust = 1.5),type="l",lty=5,col=3)
points(density(resultsA2T3[,1],adjust = 1.5),type="l",lty=4,col=4)
abline(v=parametrosA2[1],col=6)
legend("topright",c("T=500","T=1000","T=1500"),lty=c(1,5,4),col=c(2,3,4))

plot(density(resultsA2T1[,2],adjust = 1.5),lty=1,col=2,ylim=c(0,8),xlab=TeX("$\\lambda_+$"),main="", cex.lab = 1.5)
points(density(resultsA2T2[,2],adjust = 1.5),type="l",lty=5,col=3)
points(density(resultsA2T3[,2],adjust = 1.5),type="l",lty=4,col=4)
abline(v=parametrosA2[2],col=6)
legend("topright",c("T=500","T=1000","T=1500"),lty=c(1,5,4),col=c(2,3,4))

plot(density(resultsA2T1[,3],adjust = 1.5),lty=1,col=2,ylim=c(0,7),xlab=TeX("$\\lambda_-$"),main="", cex.lab = 1.5)
points(density(resultsA2T2[,3],adjust = 1.5),type="l",lty=5,col=3)
points(density(resultsA2T3[,3],adjust = 1.5),type="l",lty=4,col=4)
abline(v=parametrosA2[3],col=6)
legend("topright",c("T=500","T=1000","T=1500"),lty=c(1,5,4),col=c(2,3,4))

dev.off()


# fig1 --------------------------------------------------------------------

jpeg(filename = "./output/fig1.jpeg",width=3000,height=1000,res=350)
par(mfrow=c(1,3), lwd=2.5)
plot(density(resultsA2T1[,1],adjust = 1.5),lty=1,col=2,ylim=c(0,5),xlab=TeX("$\\alpha$"),main="", cex.lab = 1.5)
points(density(resultsA2T2[,1],adjust = 1.5),type="l",lty=5,col=3)
points(density(resultsA2T3[,1],adjust = 1.5),type="l",lty=4,col=4)
abline(v=parametrosA2[1],col=6)
legend("topright",c("T=500","T=1000","T=1500"),lty=c(1,5,4),col=c(2,3,4))

plot(density(resultsA2T1[,2],adjust = 1.5),lty=1,col=2,ylim=c(0,8),xlab=TeX("$\\lambda_+$"),main="", cex.lab = 1.5)
points(density(resultsA2T2[,2],adjust = 1.5),type="l",lty=5,col=3)
points(density(resultsA2T3[,2],adjust = 1.5),type="l",lty=4,col=4)
abline(v=parametrosA2[2],col=6)
legend("topright",c("T=500","T=1000","T=1500"),lty=c(1,5,4),col=c(2,3,4))

plot(density(resultsA2T1[,3],adjust = 1.5),lty=1,col=2,ylim=c(0,7),xlab=TeX("$\\lambda_-$"),main="", cex.lab = 1.5)
points(density(resultsA2T2[,3],adjust = 1.5),type="l",lty=5,col=3)
points(density(resultsA2T3[,3],adjust = 1.5),type="l",lty=4,col=4)
abline(v=parametrosA2[3],col=6)
legend("topright",c("T=500","T=1000","T=1500"),lty=c(1,5,4),col=c(2,3,4))
dev.off()



# fig 2 -------------------------------------------------------------------

jpeg(filename = "./output/fig2.jpeg",width=3000,height=2000,res=350)
par(mfrow=c(2,3), lwd=2)
hist(resultsA2T1[,1],ylim=c(0,3),freq=FALSE,xlab=TeX("$\\alpha$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[1],col=3)
curve(dnorm(x,mean(resultsA2T1[,1]),sd(resultsA2T1[,1])),add=T,col=2)
hist(resultsA2T1[,2],freq=FALSE,xlab=TeX("$\\lambda_+$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[2],col=3)
curve(dnorm(x,mean(resultsA2T1[,2]),sd(resultsA2T1[,2])),add=T,col=2)
hist(resultsA2T1[,3],ylim=c(0,3.5),freq=FALSE,xlab=TeX("$\\lambda_-$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[3],col=3)
curve(dnorm(x,mean(resultsA2T1[,3]),sd(resultsA2T1[,3])),add=T,col=2)

qqPlot(resultsA2T1[,1],ylab=TeX("$\\alpha$"))
qqPlot(resultsA2T1[,2],ylab=TeX("$\\lambda_+$"))
qqPlot(resultsA2T1[,3],ylab=TeX("$\\lambda_-$"))

dev.off()


jpeg(filename = "./output/fig3.jpeg",width=3000,height=2000,res=350)
par(mfrow=c(2,3), lwd=2)
hist(resultsA2T2[,1],freq=FALSE,xlab=TeX("$\\alpha$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[1],col=3)
curve(dnorm(x,mean(resultsA2T2[,1]),sd(resultsA2T2[,1])),add=T,col=2)
hist(resultsA2T2[,2],freq=FALSE,xlab=TeX("$\\lambda_+$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[2],col=3)
curve(dnorm(x,mean(resultsA2T2[,2]),sd(resultsA2T2[,2])),add=T,col=2)
hist(resultsA2T2[,3],freq=FALSE,xlab=TeX("$\\lambda_-$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[3],col=3)
curve(dnorm(x,mean(resultsA2T2[,3]),sd(resultsA2T2[,3])),add=T,col=2)

qqPlot(resultsA2T2[,1],ylab=TeX("$\\alpha$"))
qqPlot(resultsA2T2[,2],ylab=TeX("$\\lambda_+$"))
qqPlot(resultsA2T2[,3],ylab=TeX("$\\lambda_-$"))

dev.off()



jpeg(filename = "./output/fig4.jpeg",width=3000,height=2000,res=350)
par(mfrow=c(2,3), lwd=2)


hist(resultsA2T3[,1],ylim=c(0,5),freq=FALSE,xlab=TeX("$\\alpha$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[1],col=3)
curve(dnorm(x,mean(resultsA2T3[,1]),sd(resultsA2T3[,1])),add=T,col=2)
hist(resultsA2T3[,2],freq=FALSE,xlab=TeX("$\\lambda_+$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[2],col=3)
curve(dnorm(x,mean(resultsA2T3[,2]),sd(resultsA2T3[,2])),add=T,col=2)
hist(resultsA2T3[,3],ylim=c(0,6),freq=FALSE,xlab=TeX("$\\lambda_-$"),main="", cex.lab = 1.5,col=0)
abline(v=parametrosA2[3],col=3)
curve(dnorm(x,mean(resultsA2T3[,3]),sd(resultsA2T3[,3])),add=T,col=2)


qqPlot(resultsA2T3[,1],ylab=TeX("$\\alpha$"))
qqPlot(resultsA2T3[,2],ylab=TeX("$\\lambda_+$"))
qqPlot(resultsA2T3[,3],ylab=TeX("$\\lambda_-$"))

dev.off()
