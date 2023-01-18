setwd("~/Rstudio files")
source("utils.R")
source("shifted.pv.R")
library('genlasso')
zs<-read.csv("C:/Users/78469/Desktop/fmri_slice_zscores.csv")
zmatrix<-as.matrix(zs)
zvec<-as.vector(zmatrix)
m<-length(zvec)
D = genlasso::getD2dSparse(nrow(zmatrix),ncol(zmatrix))
chol_factor = Matrix::Cholesky(Matrix::crossprod(D) + Matrix::Diagonal(m), perm=TRUE, super=FALSE)
pv.vec<-2*pnorm(-abs(zvec), 0, 1)
#LAWS procedure
bh.th<-bh.func(pv.vec, 0.9)$th
pis.hat<-pis_2D.func(zmatrix, tau=bh.th, h=4.5)
#Iterative optimize with a guess of prior_prob
prior_prob<-epsest.func(zvec,0,1)
rel_tol = 1e-4
lambda=0.38
ks=1
kold=0
passcounter=1
while(abs(kold-ks)>rel_tol && passcounter<10) {
  kold<-ks
  ks<-likelihood(pv.vec, prior_prob)
  f0 = 1
  f1 = ks*exp(-ks*pv.vec)
  beta_hat = rep(0,m)
  fl0 = list(x = rep(mean(zvec), m), # likelihood term
             z = rep(0, m), # slack variable for likelihood
             r = rep(0, nrow(D)), # penalty term
             s = rep(0, nrow(D)), # slack variable for penalty
             u_dual = rep(0, m), # scaled dual variable for constraint x = z
             t_dual = rep(0,nrow(D))) # scaled dual variable for constraint r = s
  travel = 1
  prior_prob = ilogit(beta_hat)
  old_objective = sum(-log(prior_prob*f1 + (1-prior_prob)*f0))
  converged = FALSE
  times=1
  while(!converged && times<100) {

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
  passcounter<-passcounter+1
}
q.vec=seq(from=0,to=0.02,by=0.002)
law.de=rep(0,11)
shifted.de=rep(0,11)
bh.de=rep(0,11)
for (i in 1:11) {
  q=q.vec[i]
  law.res<-law.func(pv.vec, pis.hat, q)
  law.de[i]<-sum(law.res$de)
  shifted.res<-shifted.pv(pv.vec, prior_prob, q, ks)
  shifted.de[i]<-sum(shifted.res$de)
  bh.res<-bh.func(pv.vec,q)
  bh.de[i]<-sum(bh.res$de)
}

#figure
pdf("real_data.pdf")
plot(q.vec,shifted.de,type='o',pch=2,col=2,lwd=2,xlim = c(0.00,0.02),ylim=c(0,500),
     xlab = 'Target FDR level',ylab='Discoveries')
lines(q.vec,law.de,col=4,pch=4,lwd=2)
points(q.vec,law.de,pch=4,col=4,lwd=2)
lines(q.vec,bh.de,col=5,pch=5,lwd=2)
points(q.vec,bh.de,pch=5,col=5,lwd=2)
legend('topleft',c('LASS.dd','LAWS.dd','BH'),pch=c(2,4,5),col=c(2,4,5),lwd=2)
dev.off()