##### Certain functions in references LAWS and FDRsmoothing,Written by their authors
bh.func<-function(pv, q)
{ 
  # the input 
    # pv: the p-values
    # q: the FDR level
  # the output 
    # nr: the number of hypothesis to be rejected
    # th: the p-value threshold
    # de: the decision rule

  m=length(pv)
  st.pv<-sort(pv)   
  pvi<-st.pv/1:m
  de<-rep(0, m)
  if (sum(pvi<=q/m)==0)
  {
    k<-0
    pk<-1
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    de[which(pv<=pk)]<-1
  }
  y<-list(nr=k, th=pk, de=de)
  return (y)
}

law.func<-function(pvs, pis, q)
{
    ## implementing "spatial multiple testing by locally adaptive weighting"
	## Arguments
	 # pvs: p-values
	 # pis: conditional probabilities
	 # q: FDR level
	## Values
	 # de: the decision
	 # th: the threshold for weighted p-values
	 
	m<-length(pvs)
	nu<-10e-5
	pis[which(pis<nu)]<-nu # stabilization
	pis[which(pis>1-nu)]<-1-nu # stabilization
	ws<-pis/(1-pis)
    pws<-pvs/ws
    st.pws<-sort(pws)
    fdps<-sum(pis)*st.pws/(1:m)
    de<-rep(0, m)
    if(sum(fdps<=q)==0)
    {
    	k<-0
    	pwk<-1
    }
    else
    {
    	k<-max(which(fdps<=q))
    	pwk<-st.pws[k]
    	de[which(pws<=pwk)]<-1
    }
    y<-list(nr=k, th=pwk, de=de)
	return (y)
}	 

epsest.func <- function(x,u,sigma)
{
  # x is a vector of normal variables
  # u is the mean 
  # sigma is the standard deviation
  # the output is the estimated non-null proportion
  
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)

  epsest=NULL

  for (j in 1:length(tt)) { 

    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi

    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    } 
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

pis_1D.func<- function(x, tau=0.1, h=50)
{
  ## pis_est.func calculates the conditional proportions pis
  ## Arguments
   # x: z-values
   # tau: the screening threshold, which can be prespecified or chosen adaptively
   # bdw: bandwidth
  ## Values
   # pis: the conditional proportions

  m <- length(x)
  s <- 1:m # auxiliary variable
  pval <- 2*pnorm(-abs(x))
  p.est <-rep(0, m)
  for (i in 1:m) { 
  	kht<-dnorm(s-i, 0, h)
  	p.est[i]<-sum(kht[which(pval>=tau)])/((1-tau)*sum(kht))
  }
  p.est[which(p.est>1)] <-1
  return(1-p.est)
}

disvec.func<-function(dims, s)
{
   # disvec computes the distances of all points on a m1 times m2 spatial domain to a point s
  ## Arguments:
   # dims=c(d1, d2): the dimensions
   # s=c(s1, s2): a spatial point
  ## Values:
   # a vector of distances
   
   m<-dims[1]*dims[2]
   dis.vec<-rep(0, m)
   for(i in 1:dims[1])
   {
   	dis.vec[((i-1)*dims[2]+1):(i*dims[2])]<-sqrt((i-s[1])^2+(1:dims[2]-s[2])^2)
   }
   return(dis.vec) 
}

pis_2D.func<- function(x, tau=0.1, h=10)
{
  ## pis_2D.func calculates the conditional proportions pis
  ## Arguments
   # x: a matrix of z-values
   # tau: the screening threshold, which can be prespecified or chosen adaptively
   # bdw: bandwidth
  ## Values
   # pis: conditional proportions

  dims<-dim(x)
  m<-dims[1]*dims[2]
  x.vec<-c(t(x))
  pv.vec<-2*pnorm(-abs(x.vec))
  scr.idx<-which(pv.vec>=tau)
  p.est<-matrix(rep(0, m), dims[1], dims[2])  
  
  for (i in 1:dims[1]) 
  {
  	for (j in 1:dims[2]) 
  	{
  		s<-c(i, j)
  		dis.vec<-disvec.func(dims, s)
  		kht<-dnorm(dis.vec, 0, h)
    	p.est[i,j]<-min(1-1e-5, sum(kht[scr.idx])/((1-tau)*sum(kht)))
  	}
  }
  pis.est<-1-c(p.est)
  return(pis.est)
}
fit_graphfusedlasso_cholcache = function(y, lambda, D, chol_factor = NULL, weights=NULL, initial_values = NULL, iter_max = 10000, rel_tol = 1e-4, alpha=1.0, inflate=2, adaptive=FALSE) {
  require(Matrix)
  
  n = length(y)
  m = nrow(D)
  a = 2*lambda # step-size parameter
  
  if(missing(weights)) {
    weights = rep(1, n)
  }
  
  # Check if we need a Cholesky decomp of system involving graph Laplacian
  if(missing(chol_factor)) {
    L = Matrix::crossprod(D)
    chol_factor = Matrix::Cholesky(L + Matrix::Diagonal(n))
  }
  
  # Initialize primal and dual variables from warm start
  if(missing(initial_values)) {
    x = rep(0, n) # likelihood term
    z = rep(0, n) # slack variable for likelihood
    r = rep(0, m) # penalty term
    s = rep(0, m) # slack variable for penalty
    u_dual = rep(0,n) # scaled dual variable for constraint x = z
    t_dual = rep(0,m) # scaled dual variable for constraint r = s
  } else {
    x = initial_values$x
    z = initial_values$z
    r = initial_values$r
    s = initial_values$s
    t_dual = initial_values$t_dual
    u_dual = initial_values$u_dual
  }
  
  primal_trace = NULL
  dual_trace = NULL
  converged = FALSE
  counter = 0
  while(!converged & counter < iter_max) {
    
    # Update x
    x = {weights * y + a*(z - u_dual)}/{weights + a}
    x_accel = alpha*x + (1-alpha)*z
    
    # Update constraint term r
    arg = s - t_dual
    if(adaptive) {
      local_lambda = 1/{1+(lambda)*abs(arg)}  # Minimax-concave penalty instead?
    } else {
      local_lambda = lambda
    }
    r = softthresh(arg, local_lambda/a)
    r_accel = alpha*r + (1-alpha)*s
    
    # Projection to constraint set
    arg = x_accel + u_dual + Matrix::crossprod(D, r_accel + t_dual)
    z_new = drop(Matrix::solve(chol_factor, arg))
    s_new = as.numeric(D %*% z_new)
    dual_residual_u = a*(z_new - z)
    dual_residual_t = a*(s_new - s)
    z = z_new
    s = s_new
    
    # Dual update
    primal_residual_x = x_accel - z
    primal_residual_r = r_accel - s
    u_dual = u_dual + primal_residual_x
    t_dual = t_dual + primal_residual_r
    
    # Check convergence
    primal_resnorm = sqrt(mean(c(primal_residual_x, primal_residual_r)^2))
    dual_resnorm = sqrt(mean(c(dual_residual_u, dual_residual_t)^2))
    if(dual_resnorm < rel_tol && primal_resnorm < rel_tol) {
      converged=TRUE
    }
    primal_trace = c(primal_trace, primal_resnorm)
    dual_trace = c(dual_trace, dual_resnorm)
    counter = counter+1
    
    # Update step-size parameter based on norm of primal and dual residuals
    if(primal_resnorm > 5*dual_resnorm) {
      a = inflate*a
      u_dual = u_dual/inflate
      t_dual = t_dual/inflate
    } else if(dual_resnorm > 5*primal_resnorm) {
      a = a/inflate
      u_dual = inflate*u_dual
      t_dual = inflate*t_dual
    }
  }
  list(x=x, r=r, z=z, s=s, u_dual=u_dual, t_dual=t_dual,
       primal_trace = primal_trace, dual_trace=dual_trace, counter=counter)
}
softthresh = function(x, lambda) {
  return(sign(x)*pmax(0, abs(x) - lambda))
}