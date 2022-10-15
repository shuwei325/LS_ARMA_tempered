source("./00_packages.R")

tabla<-read.table("./data/TS.CF.OFFSHORE.30yr.date.txt",header=TRUE)
tabla1<- tabla %>% 
  dplyr::filter(Year==2015|Year==2014|Year==2013|Year==2012)

#Select 2015:

tabla.DE<-tabla1[,c(1:7,10)]
DE.list<-list()
for(i in 1:4){
  anos <- 2012:2015
  DE.list[[i]]<- tabla.DE %>% 
    dplyr::filter(Year==anos[i])
}
par(mfrow=c(1,2))
for(i in 1:4){
  #y<-DE.list[[i]]$DE
  y<-diff(DE.list[[i]]$DE)
}
i<-4  #2015
z.2015<-DE.list[[i]]$DE
y.2015<-diff(z.2015)

z =z.2015[3985:4992]
y = diff(z)

pdf("./output/wind.pdf")

par(mfrow=c(2,1))
ts.plot(z,ylab=expression(paste("wind power (",y[t],")")))
points(3985:4992,z.2015[3985:4992],col=2,type="l")
ts.plot(y,ylab=expression(paste(Delta,"y"[t])))
points(3985:4992,y.2015[3985:4992],col=2,type="l")


par(mfrow=c(2,1))
acf(y,lag.max = 30,main="")
pacf(y,lag.max = 30,main="")

dev.off()

# Analysis by blocks of phi and sigma parameters
T. = length(y)
N = 200
S = 100
M = trunc((T. - N)/S + 1)

par(mfrow=c(3,2))
for(j in 1:M){
  x = y[(1 + S*(j-1)):(N + S*(j-1))]
  acf(x, lag = 30, ylim = c(-0.5,0.5), bty = "n", las = 1, main = "")
  pacf(x, lag = 10, ylim = c(-1,1), bty = "n", las = 1, main = "")
}

###tvAR(4)
table = c()
for(j in 1:M){
  x = y[(1 + S*(j-1)):(N + S*(j-1))]
  table = rbind(table,nlminb(start = c(0,0,0,0, 10), N = N, objective = LS.whittle.loglik,
                             series = x, order = c(p = 4, q = 0))$par)
}
u = (N/2 + S*(1:M-1))/T.
table = as.data.frame(cbind(u, table))
colnames(table) = c("u", "theta1","theta2","theta3","theta4", "sigma")
spar = 0.6
par(mfrow=c(2,3))
plot(smooth.spline(theta1, spar = spar)$y ~ u, data = table, pch = 20,
     xlim = c(0,1), ylim=c(0,3),las = 1, bty = "n", ylab = expression(alpha[1](u)),type="b")
plot(smooth.spline(theta2, spar = spar)$y ~ u, data = table, pch = 20,
     xlim = c(0,1), ylim=c(-2,0),las = 1, bty = "n", ylab = expression(alpha[1](u)),type="b")
plot(smooth.spline(theta3, spar = spar)$y ~ u, data = table, pch = 20,
     xlim = c(0,1), ylim=c(0,2),las = 1, bty = "n", ylab = expression(alpha[1](u)),type="b")
plot(smooth.spline(theta4, spar = spar)$y ~ u, data = table, pch = 20,
     xlim = c(0,1), ylim=c(-1,1),las = 1, bty = "n", ylab = expression(alpha[1](u)),type="b")
plot(smooth.spline(sigma, spar = spar)$y ~ u, data = table, pch = 20,
     xlim = c(0,1),ylim=c(0,0.05), las = 1, bty = "n", ylab = expression(gamma(u)),type="b")

theta1 = smooth.spline(table$theta1, spar = spar)$y
fit.1a = nls(theta1 ~ a0+a1*u, start = list(a0 = 0,a1 = 0))
theta2 = smooth.spline(table$theta2, spar = spar)$y
fit.1b = nls(theta2 ~ a0+a1*u, start = list(a0 = 0,a1 = 0))
theta3 = smooth.spline(table$theta3, spar = spar)$y
fit.1c = nls(theta3 ~ a0+a1*u, start = list(a0 = 0,a1 = 0))
theta4 = smooth.spline(table$theta4, spar = spar)$y
fit.1d = nls(theta4 ~ a0+a1*u, start = list(a0 = 0,a1 = 0))
sigma = smooth.spline(table$sigma, spar = spar)$y
fit.2 = nls(sigma ~ b0+b1*u, start = list(b0 = 0.65, b1 = 0))


fit<-LS.whittle(y, start=c(coef(fit.1a),coef(fit.1b),coef(fit.1c),coef(fit.1d),coef(fit.2)),
                order = c(p = 4, q = 0),ar.order=c(1,1,1,1),sd.order = c(1))
LS.summary(fit)


par(mfrow=c(2,2))
ts.plot(fit$residuals)
qqPlot(fit$residuals)
acf(fit$residuals)
pacf(fit$residuals)


jpeg(filename = "fig10.jpeg",width=2000,height=2000,res=350)
ts.diag(fit$residuals,lag=15)
dev.off()

Box.test(fit$residuals,type="Ljung",lag=15)


skewness(fit$residuals)
kurtosis(fit$residuals)
shapiro.test(fit$residuals)
jarque.bera.test(fit$residuals)


source("./code/estimTS.R")

initial.sim.TS<-c(0.8,1,1)
mod.residuals<-stdCTS.estimate.nlminb(initial.sim=initial.sim.TS,y=fit$residuals,n=2^10)
mod.residuals


pdf("./output/whittle_ar4_residual.pdf")

boxplot(fit$residuals)
hist(fit$residuals,freq=FALSE,ylim=c(0,1),breaks=60,main="",xlab="residual")
curve(dstable(x,alpha=1.42,beta=0,gamma=0.40,delta=0,pm=1),add=T,col=2,lwd=3)
curve(dnorm(x,mean(fit$residuals),sd(fit$residuals)),add=T,col=3,lwd=3)

d<-characteristic_function_to_density(n=2^10,a=-10, b=10, alpha=mod.residuals[1],lambda.p=mod.residuals[2], lambda.n=mod.residuals[3])
points(d$x,d$density,col=4,type="l",lwd=3)
legend("topleft",c("stable","normal","stdCTS"),col=c(2,3,4),lty=c(1),lwd=3)


dev.off()



jpeg(filename = "fig9.jpeg",width=2000,height=2000,res=350)

hist(fit$residuals,freq=FALSE,ylim=c(0,1),breaks=60,main="",xlab="residual")
curve(dstable(x,alpha=1.42,beta=0,gamma=0.40,delta=0,pm=1),add=T,col=2,lwd=2.5,lty=5)
curve(dnorm(x,mean(fit$residuals),sd(fit$residuals)),add=T,col=3,lwd=2.5,lty=3)

d<-characteristic_function_to_density(n=2^10,a=-10, b=10, alpha=mod.residuals[1],lambda.p=mod.residuals[2], lambda.n=mod.residuals[3])
points(d$x,d$density,col=4,type="l",lwd=2.5)
legend("topleft",c("stable","normal","stdCTS"),col=c(2,3,4),lty=c(5,3,1),lwd=2.5)


dev.off()

