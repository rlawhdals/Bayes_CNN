#### load packages ####
source('packages.R')
sourceCpp('ergm.cpp')

#### ERGM example ####
# load dataset
data(faux.mesa.high)
summary(faux.mesa.high)
plot(faux.mesa.high)
X = as.matrix(faux.mesa.high)
dim(X)
table(X)

# Define model formula
m_formula = faux.mesa.high ~ edges + nodematch("Grade",diff=T) + gwdegree(0.25,fixed=TRUE) + gwesp(0.25,fixed=TRUE)
grade <- faux.mesa.high %v% "Grade" 
sex <- faux.mesa.high %v% "Sex"  ; sex[sex=="M"] = 1 ; sex[sex=="F"] = 0 ;sex <- as.numeric(sex)
Xsummary = Summary(X, grade, sex)
Xsummary
stat9 = summary(m_formula) ; stat9

# calculate MPLE
mple.ergm <-ergm(m_formula,estimate="MPLE")
hat <- coef(mple.ergm) ; hat
COV <- solve(-mple.ergm$hessian)


##--------------------##
##### Baseline AVM #####
##--------------------##
thin = 1
outer = 5e4 * thin
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1

ptm <- proc.time()
ergm_avm = ergm_AVM(X, grade, sex, COV, matrix(hat,1), outer, 10,  adaptInterval, adaptFactorExponent, adapIter, thin)
avm_time = proc.time()-ptm ; avm_time
save(ergm_avm, file='ergm_avm.Rdata')


##--------------------------------------##
##### DA-AVM with Function emulation #####
##--------------------------------------##
### Make ABC particles ###
p <- 9 # dimension of parameter
std.m = as.vector(summary(mple.ergm)$coefficient[,2])

const = 3 # hyper parameter
D1 <- matrix(c(mple.ergm$coef[1] - const*std.m[1], mple.ergm$coef[1] + const*std.m[1],
                    mple.ergm$coef[2] - const*std.m[2], mple.ergm$coef[2] + const*std.m[2],
                    mple.ergm$coef[3] - const*std.m[3], mple.ergm$coef[3] + const*std.m[3],
                    mple.ergm$coef[4] - const*std.m[4], mple.ergm$coef[4] + const*std.m[4],
                    mple.ergm$coef[5] - const*std.m[5], mple.ergm$coef[5] + const*std.m[5],
                    mple.ergm$coef[6] - const*std.m[6], mple.ergm$coef[6] + const*std.m[6],
                    mple.ergm$coef[7] - const*std.m[7], mple.ergm$coef[7] + const*std.m[7],
                    mple.ergm$coef[8] - const*std.m[8], mple.ergm$coef[8] + const*std.m[8],
                    mple.ergm$coef[9] - const*std.m[9], mple.ergm$coef[9] + const*std.m[9]),2,p)

num.point <- 3000    
th9 <- matrix(0,num.point,p)
dim(th9)

Design <- lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
th9 <- Design$design
hist(th9[,1])
hist(th9[,2])

th9[,1] <- (D1[2,1]-D1[1,1])*th9[,1]+D1[1,1]
th9[,2] <- (D1[2,2]-D1[1,2])*th9[,2]+D1[1,2] 
th9[,3] <- (D1[2,3]-D1[1,3])*th9[,3]+D1[1,3]
th9[,4] <- (D1[2,4]-D1[1,4])*th9[,4]+D1[1,4]
th9[,5] <- (D1[2,5]-D1[1,5])*th9[,5]+D1[1,5]
th9[,6] <- (D1[2,6]-D1[1,6])*th9[,6]+D1[1,6] 
th9[,7] <- (D1[2,7]-D1[1,7])*th9[,7]+D1[1,7]
th9[,8] <- (D1[2,8]-D1[1,8])*th9[,8]+D1[1,8]
th9[,9] <- (D1[2,9]-D1[1,9])*th9[,9]+D1[1,9]
hist(th9[,1])

N <- 1    
cycle <- 1
n_chain = 10

ptm <- proc.time()
summat <- pAuxSamp_ergm(X, grade, sex, cycle,num = 10, numcore = 10, Designmat = th9) # Generate summary statistics of auxiliary variables for given particles
s2time = proc.time()-ptm ;s2time 
summat9 = summat$Mu
dim(summat9)

dist = sqrt(  apply(  ( matrix( rep(stat9,num.point), num.point, byrow = T) - summat9[1:num.point,] )^2, 1, sum ) )

hist(dist,breaks = 100)
summary(dist)
eps = quantile(dist,0.1)
nrow(th9[which(dist < eps),])

m9 =  apply(  th9[which(dist < eps),], 2, mean)
S9 = cov( th9[which(dist < eps),] )
thABC11 = mvrnorm( 400, m9, S9)
m9 ; S9
dim(thABC11)

D2 <- rbind(apply(thABC11,2,min),apply(thABC11,2,max))

num.point <- 400 # desired number of ABC particles
th9 <- matrix(0,num.point,p)

# Design for model 
ptm <- proc.time()
Design <- lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=11)
proc.time()- ptm
th9 <- Design$design
th9[,1] <- (D2[2,1]-D2[1,1])*th9[,1]+D2[1,1]
th9[,2] <- (D2[2,2]-D2[1,2])*th9[,2]+D2[1,2]
th9[,3] <- (D2[2,3]-D2[1,3])*th9[,3]+D2[1,3]
th9[,4] <- (D2[2,4]-D2[1,4])*th9[,4]+D2[1,4]
th9[,5] <- (D2[2,5]-D2[1,5])*th9[,5]+D2[1,5]
th9[,6] <- (D2[2,6]-D2[1,6])*th9[,6]+D2[1,6] 
th9[,7] <- (D2[2,7]-D2[1,7])*th9[,7]+D2[1,7]
th9[,8] <- (D2[2,8]-D2[1,8])*th9[,8]+D2[1,8]
th9[,9] <- (D2[2,9]-D2[1,9])*th9[,9]+D2[1,9]

### Construct GP emulator ###
N <- 1000      # number of importance sampling estimate
cycle <- 10     # number of inner sampler
hat9 <- apply(thABC11,2,mean)

ptm <- proc.time()
Sample9 = pResponse_ergm(X, cycle, hat9, N, 
                         grade, sex, 20) # Generate set of auxiliary variables for each particle point
IStime = proc.time()-ptm ; IStime

y9 <- c() # estimate z(theta)
for(i in 1:num.point){  
  cal9 = Sample9%*%(th9[i,]-hat9) 
  mx9 = max(cal9)
  y9[i] = mx9 + log( mean( exp(cal9-mx9) ) ) # importance sampling estimate
}

# unnormalized likelihood 
lhX9 <- c()
for(i in 1:num.point){ lhX9[i] = stat9%*%th9[i,] }

# full likelihood 
y9 <- lhX9 - y9

theta9 <- matrix(hat9,1)
ptm <- proc.time()
m9 <- km(~ ., design = th9[1:num.point,], response = matrix(y9[1:num.point],num.point,1),covtype = "matern3_2")
GPtime9 = proc.time()-ptm ; GPtime9

x.point9 <- data.frame( t(theta9[1,]) )

# Kriging, step 4 & 5
pred.m9 <- predict(m9,x.point9, "UK")
lhXZ9 <- ypred9 <- pred.m9$mean

# starting cov matrix for proposal
COV_th9 <- cov(thABC11)

# GLS bets coefficients 
beta.hat.gls9 <- c(coef(m9)$trend1,coef(m9)$trend2,coef(m9)$trend3,coef(m9)$trend4,
                   coef(m9)$trend5,coef(m9)$trend6,coef(m9)$trend7,coef(m9)$trend8,
                   coef(m9)$trend9,coef(m9)$trend10)

ptm = proc.time()
ergm_emul = ergm_GPEm(25000,theta9,COV_th9, lhXZ9,beta.hat.gls9,phi=coef(m9)$range,coef(m9)$sd2,th9,y9)
gp_multi_time9 = proc.time() - ptm ; gp_multi_time9
round(colMeans(ergm_emul$theta[-(1:1e4),]),2)

### DA-AVM with GP emulator ###
ptm <- proc.time()
ergm_gp <- ergm_DAAVM_GP(X=X, COV=COV_th9, theta=theta9, outer=outer, cycle=10,
                                        grade, sex, lhXZ=lhXZ9, betahat=beta.hat.gls9,
                                        phihat=coef(m9)$range,sigmasq = coef(m9)$sd, Designmat=th9, y=y9,
                                        adaptInterval, adaptFactorExponent, adapIter, thin)
ergm.gp.times = proc.time()-ptm ; ergm.gp.times
save(ergm_gp, file='ergm_gp.Rdata')


##-----------------------------------------##
##### DA-AVM with Frequentist estimator #####
##-----------------------------------------##
# Calculate MCML
ptm = proc.time()
m_fit_mcml <-ergm(m_formula,estimate="MLE")
mcmltime9 = proc.time()-ptm ; mcmltime9
hat9_ <- coef(m_fit_mcml)
COV9_ <- solve(-m_fit_mcml$hessian)

ptm <- proc.time()
ergm_F = ergm_DAAVM_F(X, COV9_, matrix(hat9_,1), outer = outer, cycle=10, mu = hat9_, covariance = COV9_,
                                         grade, sex, adaptInterval,adaptFactorExponent,adapIter,thin)
ergm.f.time = proc.time() - ptm ; ergm.f.time
save(ergm_F, file='ergm_f.Rdata')


##-----------------------##
#### Result comparison ####
##-----------------------##
## Posterior mean
round(colMeans(ergm_avm$theta[-(1:1e4),]),2)
round(colMeans(ergm_gp$theta[-(1:1e4),]),2)
round(colMeans(ergm_F$theta[-(1:1e4),]),2)

## a 95% HPD interval
apply(round(ergm_avm$theta[-(1:outer*0.2),],2),2,HDInterval::hdi)
apply(round(ergm_gp$theta[-(1:outer*0.2),],2),2,HDInterval::hdi)
apply(round(ergm_F$theta[-(1:outer*0.2),],2),2,HDInterval::hdi)

## Acceptance rate
ergm_avm$Accentance_rate
ergm_gp$First_acceptance_rate ; ergm_gp$Final_acceptance_rate
ergm_F$First_acceptance_rate ; ergm_F$Final_acceptance_rate

## number of auxiliary variable simulations
outer * ergm_gp$First_acceptance_rate
outer * ergm_F$First_acceptance_rate

## ESS
batchmeans::ess(ergm_avm$theta[-(1:outer*0.2)])
batchmeans::ess(ergm_gp$theta[-(1:outer*0.2)])
batchmeans::ess(ergm_F$theta[-(1:outer*0.2)])

## Eff
(1-ergm_gp$First_acceptance_rate)/(1-ergm_gp$Final_acceptance_rate)
(1-ergm_F$First_acceptance_rate)/(1-ergm_F$Final_acceptance_rate)

## Density comparison
par(mfrow=c(1,1))
plot(density(ergm_avm$theta[-(1:outer*0.2),1]),main='Density for theta1')
lines(density(ergm_gp$theta[-(1:outer*0.2),1]),col='red')
lines(density(ergm_F$theta[-(1:outer*0.2),1]),col='blue')
points(mean(ergm_avm$theta[-(1:outer*0.2),1]),0,cex=1)
points(mean(ergm_gp$theta[-(1:outer*0.2),1]),0,cex=1,col='red')
points(mean(ergm_F$theta[-(1:outer*0.2),1]),0,cex=1.5,col='blue')
legend('topleft',legend=c('AVM', 'DA-AVM_GP', 'DA-AVM_F'), col=c('black','red','blue'),lty=c(1),bty='n')

save.image(file='ergm.RData')
