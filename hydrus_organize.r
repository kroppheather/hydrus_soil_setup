######################################
#####organizes data for model run#####
######################################
library(R2OpenBUGS)

####van genutchen data####
#model directory
modDI <- c("c:\\Users\\hkropp\\Google Drive\\hydrus\\van_genut\\run7")
chain1 <- read.bugs(paste0(modDI,"\\CODAchain1.txt"))
chain2 <- read.bugs(paste0(modDI,"\\CODAchain2.txt"))
chain3 <- read.bugs(paste0(modDI,"\\CODAchain3.txt"))

chains <- rbind(chain1[[1]],chain2[[1]],chain3[[1]])	

parmS <- data.frame(Mean= apply(chains,2,mean),	
			pc2.5=apply (chains, 2, quantile,probs=0.025)	 ,
			pc97.5=apply (chains, 2, quantile,probs=0.975))
			
			