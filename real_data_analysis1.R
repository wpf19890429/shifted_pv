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
q=0.005
pv.vec<-2*pnorm(-abs(zvec), 0, 1)
#LAWS procedure
bh.th<-bh.func(pv.vec, 0.9)$th
pis.hat<-pis_2D.func(zmatrix, tau=bh.th, h=4.5)
#Iterative optimize with a guess of prior_prob
prior_prob<-epsest.func(zvec,0,1)
rel_tol = 1e-4
lambda=.38
ks=1
kold=0
passcounter=1
while(abs(kold-ks)>rel_tol && passcounter<100) {
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

law.res<-law.func(pv.vec, pis.hat, q)
law.de<-sum(law.res$de)
shifted.res<-shifted.pv(pv.vec, prior_prob, q, ks)
shifted.de<-sum(shifted.res$de)
bh.res<-bh.func(pv.vec,q)
bh.de<-sum(bh.res$de)
dematrix1<-matrix(shifted.res$de,128,128)
locs1<-which(dematrix1 == 1, arr.ind=T)
dematrix2<-matrix(law.res$de,128,128)
locs2<-which(dematrix2 == 1, arr.ind=T)
dematrix3<-matrix(bh.res$de,128,128)
locs3<-which(dematrix3 == 1, arr.ind=T)

#figure
pdf("real_data2.pdf")
mybreaks = c(seq(0,2, length=10),6,10)
mycols = c(grey(seq(1,0,length=11)^0.5))
par(mar=c(1,1,3,1),mfrow=c(2,2))
image(1:128, 1:128, abs(zmatrix), breaks=mybreaks, col=mycols,xlab=' ', ylab=' ',
      main="(a) z scores",cex.main=0.8, axes=FALSE,las=1)
image(1:128, 1:128, abs(zmatrix), breaks=mybreaks, col=mycols,xlab=' ', ylab=' ',
      main="(b) BH",cex.main=0.8, axes=FALSE,las=1)
points(locs3[,1],locs3[,2],col='brown3',pch=18)
image(1:128, 1:128, abs(zmatrix), breaks=mybreaks, col=mycols,xlab=' ', ylab=' ',
      main="(c) LAWS.dd",cex.main=0.8, axes=FALSE,las=1)
points(locs2[,1],locs2[,2],col='brown3',pch=18)
image(1:128, 1:128, abs(zmatrix), breaks=mybreaks, col=mycols,xlab=' ', ylab=' ',
      main="(d) LASS.dd",cex.main=0.8, axes=FALSE,las=1)
points(locs1[,1],locs1[,2],col='brown3',pch=18)
dev.off()