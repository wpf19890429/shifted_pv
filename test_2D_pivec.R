# Simulating spatial data
# S: 2D spatial domain with n1*n2 points

setwd("~/Rstudio files")
source("utils.R")
source("shifted.pv.R")
library('genlasso')
library('foreach')
library('doParallel')

n<-128
m<-n^2
D<-genlasso::getD2dSparse(n,n)
chol_factor <- Matrix::Cholesky(Matrix::crossprod(D) + Matrix::Diagonal(m), perm=TRUE, super=FALSE)
thetaS<-matrix(rep(0, m), n, n) # true states of nature
pis<-matrix(0, n, n)

# A and B: Two different signal area shapes
isA.func<-function(i, j)
{
  # decide whether (i, j) is in the rectangle A
  if(i>10 && i<=30 && j>10 && j<=30)
  {y <- 1}
  else if(i>50 && i<=65 && j>60 && j<=75)
  {y <- 1}
  else
  {y<-0}
  return(y)
}
isB.func<-function(i, j)
{
  if(sqrt((i-20)^2+(j-20)^2<=100))
  {y <- 1}
  else if(i>50 && i<=65 && j>60 && j<=75)
  {y <- 1}
  else
  {y<-0}
  return(y)
}

q<-0.1
pi0.vec<-seq(from=0.5, to=0.9, by=0.1)
pis<-matrix(0, n, n)
mu0<-1.5
np<-length(pi0.vec)
nrep<-100
locs<-1:m
loop_matrix = expand.grid(1:np, 1:nrep)


cl=makeCluster(18)
clusterEvalQ(cl,library('genlasso'))
valuelist=list('n','m','np','pi0.vec','q','locs')
clusterExport(cl,valuelist)
registerDoParallel(cl)
# Methods

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

law.or.fdr<-rep(0, np)
law.or.etp<-rep(0, np)

law.dd.fdr<-rep(0, np)
law.dd.etp<-rep(0, np)

shifted.or.fdr<-rep(0, np)
shifted.or.etp<-rep(0, np)

shifted.dd.fdr<-rep(0, np)
shifted.dd.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)

law.or.fdp<-matrix(rep(0, nrep*np), np, nrep)
law.or.ntp<-matrix(rep(0, nrep*np), np, nrep)

law.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
law.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)

shifted.or.fdp<-matrix(rep(0, nrep*np), np, nrep)
shifted.or.ntp<-matrix(rep(0, nrep*np), np, nrep)

shifted.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
shifted.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)

for_out1 <- foreach(iter = 1:(np*nrep), .combine='rbind') %dopar%
  {
    i<-as.numeric(unlist(loop_matrix[iter,][1]))
    j<-as.numeric(unlist(loop_matrix[iter,][2]))
    pi0<-pi0.vec[i]
    
    set.seed(iter+2000)
    for(i1 in 1:n)
    {
      for (j1 in 1:n)
      { if (isA.func(i1,j1)==1)
        {pis[i1,j1] = pi0
         thetaS[i1,j1]<-rbinom(1, 1, pi0)}
        else
        {pis[i1,j1] = 0.01
         thetaS[i1,j1]<-rbinom(1, 1, 0.01)}
      }
    }
    
    thetaS.vec<-c(thetaS)
    pis.vec = c(pis)
    
    pii<-sum(thetaS.vec)/m
    x0<-rnorm(m, mean=0, sd=1)
    x1<-rnorm(m, mean=mu0, sd=1)
    x.vec<-(1-thetaS.vec)*x0+thetaS.vec*x1
    x<-matrix(x.vec, n, n)
    pv.vec<-2*pnorm(-abs(x.vec), 0, 1)
    
    bh.th<-bh.func(pv.vec, 0.9)$th
    pis.hat<-pis_2D.func(x, tau=bh.th, h=4.5)
    kor<-likelihood(pv.vec,pis.vec)
    pi.est<-epsest.func(x.vec,0,1)
    ks<-likelihood(pv.vec, pi.est)
    lambda = .38
    rel_tol = 1e-4
    kold=0
    passcounter=1
    while(abs(ks-kold)>1e-4 && passcounter<=100) {
      f0 = 1
      f1 = ks*exp(-ks*pv.vec)
      beta_hat = rep(0,m)
      # Initialization
      fl0 = list(x = rep(mean(x.vec), m), # likelihood term
                 z = rep(0, m), # slack variable for likelihood
                 r = rep(0, nrow(D)), # penalty term
                 s = rep(0, nrow(D)), # slack variable for penalty
                 u_dual = rep(0,m), # scaled dual variable for constraint x = z
                 t_dual = rep(0,nrow(D))) # scaled dual variable for constraint r = s
      travel = 1
      prior_prob = ilogit(beta_hat)
      old_objective = sum(-log(prior_prob*f1 + (1-prior_prob)*f0)) + lambda * sum(abs(fl0$r))
      converged = FALSE
      times=1
      while(!converged && times<=100) {
        # E step
        m1 = prior_prob*f1
        m0 = (1-prior_prob)*f0
        post_prob = m1/(m1+m0)
        
        # M step: one ADMM iteration, analogous to a single Newton iteration
        weights = prior_prob*(1-prior_prob)
        y = beta_hat - (prior_prob - post_prob)/weights
        fl0 = fit_graphfusedlasso_cholcache(y, lambda=lambda, D=D, chol_factor=chol_factor, weights=weights,
                                            initial_values = fl0, rel_tol = rel_tol, alpha=1.8, adaptive=FALSE)
        beta_hat = fl0$x
        prior_prob = ilogit(beta_hat)
        
        # Check relative convergence
        new_objective = sum(-log(prior_prob*f1 + (1-prior_prob)*f0)) + lambda * sum(abs(fl0$r))
        travel = abs(new_objective - old_objective)
        old_objective = new_objective
        converged = {travel/(old_objective + rel_tol) < rel_tol}
        times=times+1
      }
      kold<-ks
      ks<-likelihood(pv.vec, prior_prob)
      passcounter<-passcounter+1
    }
    
    bh.res<-bh.func(pv.vec, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-thetaS.vec)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(thetaS.vec*bh.de)/sum(thetaS.vec)
    
    law.or.res<-law.func(pvs=pv.vec, pis.vec, q)
    law.or.de<-law.or.res$de
    law.or.fdp[i, j]<-sum((1-thetaS.vec)*law.or.de)/max(sum(law.or.de), 1)
    law.or.ntp[i, j]<-sum(thetaS.vec*law.or.de)/sum(thetaS.vec)
    
    law.dd.res<-law.func(pvs=pv.vec, pis.hat, q)
    law.dd.de<-law.dd.res$de
    law.dd.fdp[i, j]<-sum((1-thetaS.vec)*law.dd.de)/max(sum(law.dd.de), 1)
    law.dd.ntp[i, j]<-sum(thetaS.vec*law.dd.de)/sum(thetaS.vec)
    
    
    shifted.or.res<-shifted.pv(pv.vec, pis.vec, q, kor)
    shifted.or.de<-shifted.or.res$de
    shifted.or.fdp[i, j]<-sum((1-thetaS.vec)*shifted.or.de)/max(sum(shifted.or.de), 1)
    shifted.or.ntp[i, j]<-sum(thetaS.vec*shifted.or.de)/sum(thetaS.vec)
    
    shifted.dd.res<-shifted.pv(pv.vec, prior_prob, q, ks)
    shifted.dd.de<-shifted.dd.res$de
    shifted.dd.fdp[i, j]<-sum((1-thetaS.vec)*shifted.dd.de)/max(sum(shifted.dd.de), 1)
    shifted.dd.ntp[i, j]<-sum(thetaS.vec*shifted.dd.de)/sum(thetaS.vec)
    
    cbind(i,j,shifted.or.fdp[i, j],shifted.dd.fdp[i, j],law.or.fdp[i, j],law.dd.fdp[i, j],bh.fdp[i, j],
          shifted.or.ntp[i, j],shifted.dd.ntp[i, j],law.or.ntp[i, j],law.dd.ntp[i, j],bh.ntp[i, j])
  }

mean_out1<-aggregate(for_out1,by=list(for_out1[,1]),FUN='mean')
fdr1<-mean_out1[,4:8]
etp1<-mean_out1[,9:13]

write.csv(for_out1, file='2Dtest_pivec.csv')
pdf("2Dtest_pivec.pdf")

par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(pi0.vec, fdr1, type="o", pch=1:5, lwd=2, main="FDR Comparison", xlab=expression(pi[paste(2,d)]), ylab="FDR", ylim=c(0.01, 0.30))
legend("topleft", c("LASS.or","LASS.dd","LAWS.or", "LAWS.dd","BH"), pch=1:5, col=1:5, lwd=2)

matplot(pi0.vec, etp1, type="o", pch=1:5, lwd=2, main="Power Comparison", xlab=expression(pi[paste(2,d)]), ylab="Power",ylim=c(0.01, 1.00))
legend("topleft", c("LASS.or","LASS.dd","LAWS.or", "LAWS.dd","BH"), pch=1:5, col=1:5, lwd=2)

dev.off()

