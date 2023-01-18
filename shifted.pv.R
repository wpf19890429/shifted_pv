shifted.pv<-function(pv, pis, q , ks)
{

  # the input 
    # pv: the p-values
    # pis: conditional probabilities
    # q: the FDR level
  # the output 
    # nr: the number of hypothesis to be rejected
    # th: the shifted p-value threshold
    # de: the decision rule


#### shifted p_value

	m<-length(pv)
	p_shifted<-pv+log((1-pis)/pis)/ks

	ranked_p_shifted<-p_shifted[order(p_shifted)]
	ranking<-order(p_shifted)

	for(k in m:1)
	{
		c_star<-ranked_p_shifted[k]-log((1-pis)/pis)[ranking]/ks
		c_star<-adjust(c_star) #adjust c_star in [0,1]
		FDP_value<-sum(c_star*(1-pis[ranking]),na.rm=TRUE)/k
		if(FDP_value<q)
			break
	} 

	th<-ranked_p_shifted[k]
      pv.shifted.de<-rep(0, m)
	pv.shifted.de[which(p_shifted<ranked_p_shifted[k])]<-1
	nr<-sum(pv.shifted.de)

	return(list(th=th, nr=nr, de=pv.shifted.de))

}

adjust<-function(c, min=1, max=0)
{
  length(c)
  c_min<-(c+min)/2-abs(c-min)/2
  c_min_max<-(c_min+max)/2+abs(c_min-max)/2
  return(c_min_max)
}
#the log-likelihood function
llk<-function(k,pv,pis) {
  sum(log(pis*k/exp(k*pv)+1-pis))
}
#optimization for ks
likelihood<-function(pv,pis) {
  lm1<-optimize(f=llk, interval=c(0,1000), pv=pv, pis=pis, maximum=T)
  lm1$maximum
}

flogit = function(x) log(x/(1-x))
ilogit = function(x) 1/{1+exp(-x)}