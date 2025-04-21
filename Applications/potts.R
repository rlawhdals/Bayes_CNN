#### load packages ####
source('packages.R')
sourceCpp('potts.cpp')

#### Potts example ####
# generate potts samples
grid = 32 # lattice size
nsamp = 1e4*5 # sample size
ncolor = 4  # number of status
beta = 0.8 # true parameters

set.seed(1)
init_mat <- matrix(sample(ncolor, grid*grid, replace=TRUE), nrow=grid, ncol=grid)
foo = packPotts(init_mat,ncolor)
potts_sim = unpackPotts(potts(foo, c(rep(0,ncolor), beta), nbatch=50)$final)
cycle = 10
stat = potts::calc_t_full(potts_sim,ncolor)[ncolor+1]

# calculate MPLE
t_stat<-calc_t(potts_sim,ncolor=ncolor)
t_stat
t_cache_mple <- generate_t_cache(potts_sim, ncolor, t_stat, grid*grid/4, 4,
                                 fourpixel.nonoverlap)
beta.initial<-c(rep(0,ncolor),0.3)
optim.mple <- optim(beta.initial[-1], composite.ll, gr=gr.composite.ll,
                    t_stat, t_cache_mple, method="BFGS",
                    control=list(fnscale=-1))

hessian_matrix <- optimHess(optim.mple$par, composite.ll, gr.composite.ll, t_stat, t_cache_mple)

m.hat = optim.mple$par[ncolor];m.hat
v.hat = -1/hessian_matrix[ncolor,ncolor];v.hat


##--------------------##
##### Baseline AVM #####
##--------------------##
COV = 1
sigma2 = 1
outer = 5e4
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1
beta_init = 0.1 # initial beta
ptm = proc.time()
potts_avm = potts_AVM(foo, ncolor, stat, COV=v.hat, beta_init, outer, cycle, updateCOV=F, 
                     sigma2, adaptInterval, adaptFactorExponent, adapIter)
potts.avm.time = proc.time() - ptm ; potts.avm.time
save(potts_avm, file='potts_avm.Rdata')
load('potts_avm.Rdata')
mean(potts_avm$theta[-(1:1e4)])

##--------------------------------------##
##### DA-AVM with Function emulation #####
##--------------------------------------##
### Make ABC particles ###
p=1
const = 3
D1 = matrix(c(m.hat - const*sqrt(v.hat), m.hat + const*sqrt(v.hat)),2,p)
num.point <- 3000 
th.potts <- matrix(0,num.point,p)
dim(th.potts)

Design <- lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
th.potts <- Design$design
hist(th.potts)
th.potts[,1] <- (D1[2,1]-D1[1,1])*th.potts[,1]+D1[1,1]
hist(th.potts)

N <- 1    
cycle <- 100
n_chain = 10
ptm <- proc.time()
potts_summat <- pAuxSamp_potts(foo, ncolor, cycle, c(th.potts), N, n_chain) # Generate summary statistics of auxiliary variables for given particles
s2time1 = proc.time()-ptm ; s2time1
dist = sqrt((stat-potts_summat)^2)
hist(dist,breaks = 100)
summary(dist)
eps = quantile(dist,0.1)
length(th.potts[which(dist < eps),])
m1 =  mean(th.potts[which(dist < eps),])
S1 = var(th.potts[which(dist < eps),])
m1 ; S1
thABC.potts = rnorm(200, m1, sqrt(S1))

D2 <- rbind(min(thABC.potts), max(thABC.potts))
D2
num.point <- 40 # desired number of ABC particles
th11 <- matrix(0,num.point,p)
Design <- lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=11)
th11 <- Design$design
th11[,1] <- (D2[2,1]-D2[1,1])*th11[,1]+D2[1,1]
hist(th11) # Final ABC particles

### Construct GP emulator ###
N <- 1000 # number of importance sampling estimate
cycle <- 100 # number of inner sampler
hat <- mean(thABC.potts)

ptm <- proc.time()
Sample = pResponse_potts(foo,ncolor,cycle,hat,N) # Generate summary statistics of auxiliary variables for given particles to construct IS estimate
IStime2 = proc.time()-ptm ; IStime2
dim(Sample)

y <- c() # estimate z(theta)
for(i in 1:num.point){  
  cal = Sample*(th11[i]-hat) 
  mx = max(cal)
  y[i] = mx + log( mean( exp(cal-mx) ) ) # importance sampling estimate
}

# unnormalized likelihood 
lhX <- c()
for(i in 1:num.point){ lhX[i] = stat%*%th11[i,] }

# full likelihood 
y <- lhX - y ; y
cor(th11,y)

## GP emulator
Niter <- 25000
theta <- matrix(hat,1)
ptm <- proc.time()
m <- km(~ ., design = matrix(th11,ncol=1), response = matrix(y[1:num.point],num.point,1),covtype = "matern3_2")
GPtime = proc.time()-ptm ; GPtime
m

x.point <- data.frame( t(theta[1,]) )
colnames(x.point) = 'design'

# Kriging step
pred.m <- predict(m,x.point, "UK")
lhXZ <- ypred <- pred.m$mean

COV <- var(th11)

# GLS bets coefficients 
beta.hat.gls <- c(coef(m)$trend1,coef(m)$trend2)

potts.gp = potts_GPEm(10000,c(theta),sqrt(COV),c(lhXZ),betahat=beta.hat.gls,
                      phihat=c(coef(m)$range,coef(m)$sd2),Designmat=c(th11),y,stat)
mean(potts.gp)

### DA-AVM with GP emulator ###
ptm = proc.time()
potts_gp = potts_DAAVM_GP(foo, ncolor, stat, COV, hat, outer, 10, updateCOV=F, 
                            sigma2, adaptInterval, adaptFactorExponent, adapIter,
                            lhXZ,beta.hat.gls,c(coef(m)$range,coef(m)$sd2),c(th11),y)
potts.gp.times = proc.time()-ptm ; potts.gp.times
save(potts_gp, file='potts_gp.Rdata')

mean(potts_gp$theta[-(1:(outer*0.2))])
plot(density(potts_gp$theta[-(1:outer*0.2)]))
ts.plot(potts_gp$theta[-(1:outer*0.2)])

potts_gp$First_acceptance_rate
potts_gp$Final_acceptance_rate

##-----------------------------------------##
##### DA-AVM with Frequentist estimator #####
##-----------------------------------------##
# m.hat is MPLE and v.hat is the inverse of negative hessian
ptm = proc.time()
potts_F = potts_DAAVM_F(foo, ncolor, stat, v.hat, m.hat, outer, cycle, updateCOV, 
                                  sigma2, adaptInterval, adaptFactorExponent, adapIter,mu = m.hat,covariance = v.hat)
potts.f.time = proc.time() - ptm ; potts.f.time
save(potts_F, file='potts_f.Rdata')


##-----------------------##
#### Result comparison ####
##-----------------------##
## Posterior mean
mean(potts_avm$theta[-(1:outer*0.2)])
mean(potts_gp$theta[-(1:outer*0.2)])
mean(potts_F$theta[-(1:outer*0.2)])

## Trace plot
par(mfrow=c(1,3))
ts.plot(potts_avm$theta[-(1:outer*0.2)], main='AVM')
abline(h=mean(potts_avm$theta[-(1:outer*0.2)]),col='red',lty=2)
ts.plot(potts_gp$theta[-(1:outer*0.2)], main='DA-AVM_GP')
abline(h=mean(potts_gp$theta[-(1:outer*0.2)]),col='red',lty=2)
ts.plot(potts_F$theta[-(1:outer*0.2)], main='DA-AVM_F')
abline(h=mean(potts_F$theta[-(1:outer*0.2)]),col='red',lty=2)

## a 95% HPD interval
round(HDInterval::hdi(potts_avm$theta[-(1:outer*0.2)]),2)
round(HDInterval::hdi(potts_gp$theta[-(1:outer*0.2)]),2)
round(HDInterval::hdi(potts_F$theta[-(1:outer*0.2)]),2)

## Acceptance rate
mean(potts_avm$Acceptance_rate[-(1:outer*0.2)])
potts_gp$First_acceptance_rate ; potts_gp$Final_acceptance_rate
potts_F$First_acceptance_rate ; potts_F$Final_acceptance_rate

## number of auxiliary variable simulations
outer * potts_gp$First_acceptance_rate
outer * potts_F$First_acceptance_rate

## ESS
batchmeans::ess(potts_avm$theta[-(1:outer*0.2)])
batchmeans::ess(potts_gp$theta[-(1:outer*0.2)])
batchmeans::ess(potts_F$theta[-(1:outer*0.2)])

## Eff
(1-potts_gp$First_acceptance_rate)/(1-potts_gp$Final_acceptance_rate)
(1-potts_F$First_acceptance_rate)/(1-potts_F$Final_acceptance_rate)

## Density comparison
par(mfrow=c(1,1))
plot(density(potts_avm$theta[-(1:outer*0.2)]),ylim=c(0,13),main='Density for beta')
lines(density(potts_gp$theta[-(1:outer*0.2)]),col='red')
lines(density(potts_F$theta[-(1:outer*0.2)]),col='blue')
points(mean(potts_avm$theta[-(1:outer*0.2)]),0,cex=1)
points(mean(potts_gp$theta[-(1:outer*0.2)]),0,cex=1,col='red')
points(mean(potts_F$theta[-(1:outer*0.2)]),0,cex=1.5,col='blue')
legend('topleft',legend=c('AVM', 'DA-AVM_GP', 'DA-AVM_F'), col=c('black','red','blue'),lty=c(1),bty='n')

save.image(file='potts.RData')
