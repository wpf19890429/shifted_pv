## Simulation 1: mu vector
library(glmgen)
setwd("~/Rstudio files")
source("utils.R")
source("shifted.pv.R")

#Fixed parameters
m<-5000
mu0.vec<-seq(from=1, to=2, by=0.2)
pis<-rep(0.01, m)
pis[1001:1200]<-0.8
pis[2001:2200]<-0.8
pis[3001:3200]<-0.8
pis[4001:4200]<-0.8
np<-length(mu0.vec)
q<-0.1
nrep<-100
locs<-1:m
bdw0<-50
lambda=2
rel_tol=1e-4

# Methods
bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)

law.or.fdr<-rep(0, np)
law.or.etp<-rep(0, np)
law.or.atp<-rep(0, np)

law.or.fdp<-matrix(rep(0, nrep*np), np, nrep)
law.or.ntp<-matrix(rep(0, nrep*np), np, nrep)
law.or.tp<-matrix(rep(0, nrep*np), np, nrep)

law.dd.fdr<-rep(0, np)
law.dd.etp<-rep(0, np)
law.dd.atp<-rep(0, np)

law.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
law.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)
law.dd.tp<-matrix(rep(0, nrep*np), np, nrep)

pvshifted.or.fdr<-rep(0, np)
pvshifted.or.etp<-rep(0, np)
pvshifted.or.atp<-rep(0, np)

pvshifted.or.fdp<-matrix(rep(0, nrep*np), np, nrep)
pvshifted.or.ntp<-matrix(rep(0, nrep*np), np, nrep)
pvshifted.or.tp<-matrix(rep(0, nrep*np), np, nrep)

pvshifted.dd.fdr<-rep(0, np)
pvshifted.dd.etp<-rep(0, np)
pvshifted.dd.atp<-rep(0, np)

pvshifted.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
pvshifted.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)
pvshifted.dd.tp<-matrix(rep(0, nrep*np), np, nrep)

for (i in 1:np)
{
  cat("\n", "iteration i= ", i, "\n", "iteration j=")
  mu0<-mu0.vec[i]
  
  for (j in 1:nrep)
  {
    cat(j)
    #Generate data
    theta<-rbinom(m, size=1, prob=pis)
    pii<-sum(theta)/m
    x0<-rnorm(m, mean=0, sd=1)
    x1<-rnorm(m, mean=mu0, sd=1)
    x<-(1-theta)*x0+theta*x1
    pv<-2*pnorm(-abs(x), 0, 1)
    
    bh.th<-bh.func(pv, 0.9)$th
    pitaos<-pis_1D.func(x=x, tau=bh.th, h=bdw0)
    
    #Iterative optimize parameter k and prior probability pis.hat
    kor<-likelihood(pv, pis)
    pi.hat<-epsest.func(x,0,1)
    ks<-likelihood(pv, pi.hat)
    kold<-0
    passcounter=1
    while(abs(kold-ks)>rel_tol && passcounter<=100) {
      f0 <- 1
      f1 <- ks*exp(-ks*pv)
      # Initial settings for algorithm
      drift <- 1
      beta_hat <- rep(0,m)
      prior_prob <- ilogit(beta_hat)
      objective_old <- sum(log(prior_prob*f1 + (1-prior_prob)*f0),na.rm=TRUE)
      while(drift > rel_tol) {
        # E step
        prior_prob <- ilogit(beta_hat)
        m1 <- prior_prob*f1
        m0 <- (1-prior_prob)*f0
        post_prob <- m1/(m1+m0)
        
        # M step
        weights <- prior_prob*(1-prior_prob)
        y <- beta_hat - (prior_prob - post_prob)/weights
        
        # Solve the 1D FusedLasso problem using the subroutine in glmgen
        fl0 <- trendfilter(drop(y), weights=drop(weights), k = 0, family='gaussian', lambda=lambda)
        
        beta_hat <- fl0$beta
        prior_prob <- ilogit(beta_hat)
        
        objective_new <- sum(log(prior_prob*f1 + (1-prior_prob)*f0),na.rm=TRUE)
        drift <- abs(objective_old - objective_new)/(abs(objective_old) + rel_tol)
        objective_old <- objective_new
      }
      kold<-ks
      ks<-likelihood(pv, prior_prob)
      passcounter<-passcounter+1
    }
   
    #calculate the fdp,ntp 
    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)
    
    law.or.res<-law.func(pv, pis, q)
    law.or.de<-law.or.res$de
    law.or.fdp[i, j]<-sum((1-theta)*law.or.de)/max(sum(law.or.de), 1)
    law.or.ntp[i, j]<-sum(theta*law.or.de)/sum(theta)
    law.or.tp[i, j]<-sum(theta*law.or.de)
    
    law.dd.res<-law.func(pv, pitaos, q)
    law.dd.de<-law.dd.res$de
    law.dd.fdp[i, j]<-sum((1-theta)*law.dd.de)/max(sum(law.dd.de), 1)
    law.dd.ntp[i, j]<-sum(theta*law.dd.de)/sum(theta)
    law.dd.tp[i, j]<-sum(theta*law.dd.de)
    
    pvshifted.or.de<-shifted.pv(pv, pis, q, kor)$de
    pvshifted.or.fdp[i, j]<-sum((1-theta)*pvshifted.or.de)/max(sum(pvshifted.or.de), 1)
    pvshifted.or.ntp[i, j]<-sum(theta*pvshifted.or.de)/sum(theta)
    pvshifted.or.tp[i, j]<-sum(theta*pvshifted.or.de)
    
    pvshifted.dd.de<-shifted.pv(pv, prior_prob, q ,ks)$de
    pvshifted.dd.fdp[i, j]<-sum((1-theta)*pvshifted.dd.de)/max(sum(pvshifted.dd.de), 1)
    pvshifted.dd.ntp[i, j]<-sum(theta*pvshifted.dd.de)/sum(theta)
    pvshifted.dd.tp[i, j]<-sum(theta*pvshifted.dd.de)
  }
  
  bh.fdr[i]<-mean(bh.fdp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])
  
  law.or.fdr[i]<-mean(law.or.fdp[i,])
  law.or.etp[i]<-mean(law.or.ntp[i,])	
  law.or.atp[i]<-mean(law.or.tp[i,])	
  
  law.dd.fdr[i]<-mean(law.dd.fdp[i,])
  law.dd.etp[i]<-mean(law.dd.ntp[i,])	
  law.dd.atp[i]<-mean(law.dd.tp[i,])
  
  pvshifted.or.fdr[i]<-mean(pvshifted.or.fdp[i,])
  pvshifted.or.etp[i]<-mean(pvshifted.or.ntp[i,])
  pvshifted.or.atp[i]<-mean(pvshifted.or.tp[i,])	
  
  pvshifted.dd.fdr[i]<-mean(pvshifted.dd.fdp[i,])
  pvshifted.dd.etp[i]<-mean(pvshifted.dd.ntp[i,])
  pvshifted.dd.atp[i]<-mean(pvshifted.dd.tp[i,])

}
output<-rbind(bh.fdr,law.or.fdr,law.dd.fdr,pvshifted.or.fdr,pvshifted.dd.fdr,
              bh.etp,law.or.etp,law.dd.etp,pvshifted.or.etp,pvshifted.dd.etp)
write.csv(output, file='1Dtest_muvec.csv')
fdr1.mthd<-cbind(pvshifted.or.fdr, pvshifted.dd.fdr, law.or.fdr, law.dd.fdr, bh.fdr)
etp1.mthd<-cbind(pvshifted.or.etp, pvshifted.dd.etp, law.or.etp, law.dd.etp, bh.etp)

#figure
pdf("1Dtest_muvec.pdf")
par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
matplot(mu0.vec, fdr1.mthd, type="o", pch=1:5, lwd=rep(2,5), main="FDR Comparison", xlab=expression(mu), ylab="FDR", ylim=c(0.01, 0.30))
legend("topleft", c("LASS.or","LASS.dd","LAWS.or", "LAWS.dd","BH"), pch=1:5, col=1:5, lwd=2)
matplot(mu0.vec, etp1.mthd, type="o", pch=1:5, lwd=rep(2,5), main="Power Comparison", xlab=expression(mu), ylab="Power",ylim=c(0.01, 1.00))
legend("topleft", c("LASS.or","LASS.dd","LAWS.or", "LAWS.dd",'BH'), pch=1:5, col=1:5, lwd=2)
dev.off()

