model{
	for(i in 1:Nobs){
		psi[i] ~ dnorm(mu.psi[i],tau.psi)
		mu.psi[i] <- psi.r[shrubD[i]] + ((psi.s[shrubD[i]]-psi.r[shrubD[i]])/((1+((a[shrubD[i]]*h.mpa[i])^n[shrubD[i]]))^(1-(1/n[shrubD[i]])))
	}
	for(j in 1:NshrubD){
		alpha.mpa[j]
		#convert to cm
		alpah.cm[j] <- alpha.mpa[j] * 10197.1621
		n[j]
	}
	tau.psi <- pow(sig.psi,-2)
	sig.psi ~ dunif(0,10)
}