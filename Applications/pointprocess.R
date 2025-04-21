rm(list=ls())
library(sp)
library(gstat)
library(fields)
library(classInt)
library(spatstat)
library(DiceKriging)
library(DiceDesign)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp") #adjust for your openmp setting
Sys.setenv("PKG_LIBS"="-fopenmp") #adjust for your openmp setting
library(Rcpp)
library(RcppArmadillo)
library(batchmeans)
sourceCpp("pointprocess.cpp")

### prparing 1A2A experiment data
RSVB2 <- read.table("rsvb2.txt", header = TRUE, sep = "\t")
assign("RSVB2", RSVB2, envir = .GlobalEnv)  
R <- 0 
rho <- function(x){
  result <- 0 + (x>100)*1 + ((x>r1) & (x<=100))*(1 + 1/( t3*(x-r2) )^2) + 
    ((x>R) & (x<= r1)) * (t1 - ( sqrt(t1)*(x-t2)/(t2-R) )^2)
  return(result)}  
X <- as.matrix( RSVB2 )
r.X <- rdist(X)
plot(X)
center <- c(2000,2000); radius <- (1350)
W <- disc(radius, center)



##--------------------------------------##
##### Standard AVM#####
##--------------------------------------##
# calculating initial likelihood 
hat = c(3e-4,1.35,12,0.25)  
lambda = hat[1]
t1 <- hat[2]; t2 <- hat[3]; t3 <- hat[4]   # theta parameter (for interaction function)

# solution for r1 r2 solving two equations #
a <- (t2-R)^2 / (t1*t3^2)
b <- t3^2 * (t1 - 1)
d <- 27*a*b^2
c <- ( d + sqrt((d+2)^2-4) + 2 )^(1/3)
deltar <- 1/(3*b) * ( c/2^(1/3) + 2^(1/3)/c + 1 )

r1 <- a / deltar^1.5 + t2
r2 <- r1 - sqrt(deltar)

rho.X <- rho(r.X)                            
diag(rho.X) <- NA
rho.X.sum <- apply(log(rho.X),2,sum,na.rm=T)
l.h.X <- sum( pmin(rho.X.sum,1.2) ) + (dim(X)[1])*log(lambda)   # likelihood for the initial data
l.h.X

### Make Design points ###
hat                                            # initial guess  
l.h.X                                            # likelihood for the initial data
n <- 30000                       # inner sampler
COV <- diag(c(3e-10,0.01,0.1,0.01))

ptm<-proc.time()
avm<-DMH(center,radius,X,r.X,l.h.X,40000,matrix(hat,1),COV,n)
proc.time()-ptm
save(avm,file="standard_avm40000.RData")

##--------------------------------------##
##### DA-AVM with Function emulation #####
##--------------------------------------##
#short run of dmh: 1500 runs for 100 particles, 2000 runs for 200 particles, 3000 runs for 400 particles
ptm <- proc.time() 
dmh <- DMH(center, radius, X, r.X, l.h.X, 2000, matrix(hat,1), COV, n) 
proc.time() - ptm

# calculating initial likelihood 
parameter = dmh[1000:2000,]
hat = apply(parameter[1:1000,],2,mean)
hat



#getting 200 particles
p <- 4                                          # dimension of parameter 
num.point <- 200                                # suppose we have 'num.point data'
th <- matrix(0,num.point,p)
for(i in 1:4){ th[,i] = unique(parameter[,i])[1:num.point] }
length(unique(th[,1]))


### Calculate importance sampling estimate  (sampling and calculating) ###
N <- 2000 # number of importance sampling estimate
n <- 30000 # number of Birth Death MCMC run 
ptm <- proc.time()
Sample = pResponse(center,(radius),X,r.X,hat,n,th,N,20)
IStime = proc.time()-ptm 

Max = apply(Sample,2,max) 
maxx = matrix(rep(Max,N),N,byrow=T)
y = Max + log( apply( exp((Sample)-maxx), 2, mean ) )

# unnormalized likelihood 
lhX <- c()
for(i in 1:200){ lhX[i] <- evalunNormX(th[i,], r.X) }
# full likelihood
y = lhX - y 

### Final run of MCMC by using krigging estimate ###
Niter <- 40000
theta <- matrix(hat,1)
matrix(y[1:num.point],num.point,1)
# calculating initial value of normalize const #
ptm <- proc.time()

m <- km(~ ., design = th[1:num.point,], response = matrix(y[1:num.point],num.point,1),covtype = "matern3_2")
GPtime = proc.time()-ptm
GPtime
x.point <- data.frame( t(theta[1,]) )

# Kriging
pred.m <- predict(m,x.point, "UK")
lhXZ <- ypred <- pred.m$mean
lhXZ

# starting cov matrix for proposal
COV <- diag(c(3e-10,0.01,0.1,0.001))                                 

# GLS bets coefficients 
beta.hat.gls <- c(coef(m)$trend1,coef(m)$trend2,coef(m)$trend3,coef(m)$trend4,coef(m)$trend5)



#delayed-acceptance MCMC Run using GP emulator
ptm<-proc.time()
parameterGP=DelayedAcceptanceMCMC(40000,theta,COV,lhXZ,l.h.X,beta.hat.gls,c(coef(m)$range,coef(m)$sd2),th,y,r.X,center, radius, X,r.X,n,1)
MCMC2time=proc.time()-ptm ;MCMC2time
save(parameterGP,file="GPdelayed200_40000.RData")



##--------------------------------------##
##### DA-AVM with Subsampling #####
##--------------------------------------##
# (a) divide into 16 parts with angle_width pi/2
### ideal subsample
numshard <- 16
angle_step <- 360 / numshard  # Rotate by 22.5 degrees per shard
angle_width <- 90            


dx <- X[,1] - center[1]
dy <- X[,2] - center[2]


angles <- atan2(dy, dx) * 180 / pi
angles[angles < 0] <- angles[angles < 0] + 360 

Data2 <- vector("list", numshard)

for (i in 1:numshard) {
  center_angle <- (i - 1) * angle_step
  start_angle <- (center_angle - angle_width/2) %% 360
  end_angle <- (center_angle + angle_width/2) %% 360
  
  if (start_angle < end_angle) {
    idx <- which(angles >= start_angle & angles <= end_angle)
  } else {
    
    idx <- which(angles >= start_angle | angles <= end_angle)
  }
  
  Data2[[i]] <- X[idx, ]
}
l.h.X.sub.list <- vector("list", 16)
r.X.sub.list<-vector("list",16)
for (i in 1:16) {
  X.sub <- Data2[[i]]
  r.X.sub <- rdist(X.sub)
  rho.X.sub <- rho(r.X.sub)
  diag(rho.X.sub) <- NA
  rho.X.sub.sum <- apply(log(rho.X.sub), 2, sum, na.rm = TRUE)
  l.h.X.sub <- sum(pmin(rho.X.sub.sum, 1.2)) + (dim(X.sub)[1]) * log(lambda)
  l.h.X.sub.list[[i]] <- l.h.X.sub
  r.X.sub.list[[i]]<- r.X.sub
}

l.h.X.sub.list<-unlist(l.h.X.sub.list)


###(b) divide into 32 parts with angle width pi/4
numshard <- 32
angle_step <- 360 / numshard  # Rotate by 11.25 degrees per shard
angle_width <- 45             

# 중심 기준 좌표 차이
dx <- X[,1] - center[1]
dy <- X[,2] - center[2]

angles <- atan2(dy, dx) * 180 / pi
angles[angles < 0] <- angles[angles < 0] + 360

Data2 <- vector("list", numshard)

for (i in 1:numshard) {
  center_angle <- (i - 1) * angle_step
  start_angle <- (center_angle - angle_width/2) %% 360
  end_angle <- (center_angle + angle_width/2) %% 360
  
  if (start_angle < end_angle) {
    idx <- which(angles >= start_angle & angles <= end_angle)
  } else {
    idx <- which(angles >= start_angle | angles <= end_angle)
  }
  
  Data2[[i]] <- X[idx, ]
}

l.h.X.sub.list <- vector("list", 32)
r.X.sub.list<-vector("list",32)
for (i in 1:32) {
  X.sub <- Data2[[i]]
  r.X.sub <- rdist(X.sub)
  rho.X.sub <- rho(r.X.sub)
  diag(rho.X.sub) <- NA
  rho.X.sub.sum <- apply(log(rho.X.sub), 2, sum, na.rm = TRUE)
  l.h.X.sub <- sum(pmin(rho.X.sub.sum, 1.2)) + (dim(X.sub)[1]) * log(lambda)
  l.h.X.sub.list[[i]] <- l.h.X.sub
  r.X.sub.list[[i]]<- r.X.sub
}

l.h.X.sub.list<-unlist(l.h.X.sub.list)

###(c) divide into 64 parts with angle width pi/8
numshard <- 64
angle_step <- 360 / numshard  # Rotate by 5.625 degrees per shard
angle_width <- 22.5           

dx <- X[,1] - center[1]
dy <- X[,2] - center[2]

angles <- atan2(dy, dx) * 180 / pi
angles[angles < 0] <- angles[angles < 0] + 360

Data2 <- vector("list", numshard)

for (i in 1:numshard) {
  center_angle <- (i - 1) * angle_step
  start_angle <- (center_angle - angle_width/2) %% 360
  end_angle <- (center_angle + angle_width/2) %% 360
  
  if (start_angle < end_angle) {
    idx <- which(angles >= start_angle & angles <= end_angle)
  } else {
    idx <- which(angles >= start_angle | angles <= end_angle)
  }
  
  Data2[[i]] <- X[idx, ]
}

l.h.X.sub.list <- vector("list", 64)
r.X.sub.list<-vector("list",64)
for (i in 1:64) {
  X.sub <- Data2[[i]]
  r.X.sub <- rdist(X.sub)
  rho.X.sub <- rho(r.X.sub)
  diag(rho.X.sub) <- NA
  rho.X.sub.sum <- apply(log(rho.X.sub), 2, sum, na.rm = TRUE)
  l.h.X.sub <- sum(pmin(rho.X.sub.sum, 1.2)) + (dim(X.sub)[1]) * log(lambda)
  l.h.X.sub.list[[i]] <- l.h.X.sub
  r.X.sub.list[[i]]<- r.X.sub
}

l.h.X.sub.list<-unlist(l.h.X.sub.list)
###(a) sampling
n=30000
pt<-proc.time()
parameter_ideal_random<-DADMHIdealRandom(center,radius,X,r.X,l.h.X,Data2,r.X.sub.list,l.h.X.sub.list,40000,matrix(hat,1),COV,n,numshard,2)
end<-proc.time()
print(end-pt)
save(parameter_ideal_random,file="subsample4delayed40000idealrandom.RData")


###(b) sampling
n=30000
pt<-proc.time()
parameter_ideal_random2<-DADMHIdealRandom(center,radius,X,r.X,l.h.X,Data2,r.X.sub.list,l.h.X.sub.list,40000,matrix(hat,1),COV,n,numshard,4)
end<-proc.time()
print(end-pt)
parameter_ideal_random2$theta
save(parameter_ideal_random2,file="subsample8delayed40000idealrandom.RData")
###(c) sampling
n=30000
pt<-proc.time()
parameter_ideal_random3<-DADMHIdealRandom(center,radius,X,r.X,l.h.X,Data2,r.X.sub.list,l.h.X.sub.list,40000,matrix(hat,1),COV,n,numshard,8)
end<-proc.time()
print(end-pt)
save(parameter_ideal_random3,file="subsample16delayed40000idealrandom.RData")


###Result###
burn_in=10000
avm_after_burnin<-avm[(burn_in+1):nrow(dmh),]
gp_after_burnin<-parameterGP$theta[(burn_in):nrow(parameterGP$theta),]
sub4_after_burnin<-parameter_ideal_random$theta[(burn_in):nrow(parameter_ideal_random$theta),]
###Traceplot
par(mfrow=c(2,2))
ts.plot(avm_after_burnin[,1],main='AVM',ylab='lambda')
abline(h=mean(avm_after_burnin[,1]),col='red',lty=2)
ts.plot(avm_after_burnin[,2],main='AVM',ylab='theta1')
abline(h=mean(avm_after_burnin[,2]),col='red',lty=2)
ts.plot(avm_after_burnin[,3],main='AVM',ylab='theta2')
abline(h=mean(avm_after_burnin[,3]),col='red',lty=2)
ts.plot(avm_after_burnin[,4],main='AVM',ylab='theta3')
abline(h=mean(avm_after_burnin[,4]),col='red',lty=2)
ts.plot(gp_after_burnin[,1],main='GP',ylab='lambda')
abline(h=mean(gp_after_burnin[,1]),col='red',lty=2)
ts.plot(gp_after_burnin[,2],main='GP',ylab='theta1')
abline(h=mean(gp_after_burnin[,2]),col='red',lty=2)
ts.plot(gp_after_burnin[,3],main='GP',ylab='theta2')
abline(h=mean(gp_after_burnin[,3]),col='red',lty=2)
ts.plot(gp_after_burnin[,4],main='GP',ylab='theta3')
abline(h=mean(gp_after_burnin[,4]),col='red',lty=2)
ts.plot(sub4_after_burnin[,1],main='Sub',ylab='lambda')
abline(h=mean(sub4_after_burnin[,1]),col='red',lty=2)
ts.plot(sub4_after_burnin[,2],main='Sub',ylab='theta1')
abline(h=mean(sub4_after_burnin[,2]),col='red',lty=2)
ts.plot(sub4_after_burnin[,3],main='Sub',ylab='theta2')
abline(h=mean(sub4_after_burnin[,3]),col='red',lty=2)
ts.plot(sub4_after_burnin[,4],main='Sub',ylab='theta3')
abline(h=mean(sub4_after_burnin[,4]),col='red',lty=2)

## Posterior mean
colMeans(avm_after_burnin)
colMeans(gp_after_burnin)
colMeans(sub4_after_burnin)

##a 95% HPD interval
hdi(avm_after_burnin)
hdi(gp_after_burnin)
hdi(sub4_after_burnin)


## Acceptance rate
length(unique(avm[,1]))/40000
parameterGP$GP_acceptance_count/40000 ; length(unique(parameterGP$theta[,1]))/40000
(40000-parameter_ideal_random$GP_reject_count)/40000 ;length(unique(parameterGP$theta[,1]))/40000


## Number of Auxiliary Variable Simulation
parameterGP$GP_acceptance_count
40000-parameter_ideal_random$GP_reject_count


## ESS
batchmeans::ess(avm_after_burnin)
batchmeans::ess(gp_after_burnin)
batchmeans::ess(sub4_after_burnin)


## EFF
(1-parameterGP$GP_acceptance_count/40000)/(1-(length(unique(parameterGP$theta[,1]))/40000))
(1-(40000-parameter_ideal_random$GP_reject_count)/40000)/(1-length(unique(parameterGP$theta[,1]))/40000)


## Density Comparison
colors <- c("black", "blue", "red")

for (i in 1:4) {
  dens_avm <- density(avm_after_burnin[, i])
  dens_gp <- density(gp_after_burnin[, i])
  dens_sub4 <- density(sub4_after_burnin[, i])
  
  y_max <- max(dens_avm$y, dens_gp$y, dens_sub4$y)
  
  plot(dens_avm, col = colors[1], lwd = 2,
       main = paste("Parameter", i),
       xlab = paste("Theta", i),
       ylab = "Density",
       ylim = c(0, y_max))
  
  lines(dens_gp, col = colors[2], lwd = 2)
  lines(dens_sub4, col = colors[3], lwd = 2)
  
  legend("topright", legend = c("AVM", "GP", "SUB4"),
         col = colors, lwd = 2, bty = "n", cex = 0.8)
}
