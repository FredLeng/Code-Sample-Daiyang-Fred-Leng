
install.packages("SemiPar")

library(SemiPar)
data(lidar)
attach(lidar)

y<-logratio
x<-range

plot(x,y, main="LIDAR")

#Knots
K=(700-400)/12.5  ## K=24

l<-1:24
knots<-400+12.5*l  ## K=24
length(knots)

###Model Matrix
n<-length(x)
p<-length(knots)+2   #### p=26, K=24
X<-matrix(1, nrow=n, ncol=p)
X[,2]<-x
for (j in 1:(p-2)){
  X[,(j+2)]<-(x-knots[j])*(x>knots[j])	
}

#Fitting and RSS and CV
loglam=seq(-3,12,0.01)
lam=exp(loglam)
nlam=length(lam)

D<-diag(c(0,0, rep(1, p-2)))
yhat=rep(0,length(lam))

CV<-rep(0, nlam)
RSS=rep(0, nlam)
for (k in 1:nlam){
  D<-diag(c(0,0, rep(1, p-2)))
  yhat<-X%*%solve(t(X)%*%X+lam[k]^2*D)%*%t(X)%*%y
  L<-diag(X%*%solve(t(X)%*%X+lam[k]^2*D)%*%t(X))
  RSS[k]=sum((y-yhat)^2)
  CV[k]=0
  CVtemp=0
  for (i in 1:n){
  CVtemp<-((y[i]-yhat[i])/(1-L[i]))^2
  CV[k]=CV[k]+CVtemp}
}
#Plot for RSS and CV
plot(loglam,RSS, ylab="RSS and CV",pch=0.1, cex=0.1, xlim=c(-2,12), main="RSS and CV for LIDAR")
lines(loglam,CV,lty=2)

abline(v=loglam[which.min(CV)])

legend(0,3, legend=c("RSS", "CV"), lty=1:2, cex=0.8)

#Data with regression fit

plot(x,y, main="LIDAR with Spline Fit")

loglambdamin=loglam[which.min(CV)]

lambda=exp(loglambdamin)

yhatfinal<-X%*%solve(t(X)%*%X+lambda^2*D)%*%t(X)%*%y

lines(x,yhatfinal)
