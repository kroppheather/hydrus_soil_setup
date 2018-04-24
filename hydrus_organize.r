######################################
#####organizes data for model run#####
######################################
library(R2OpenBUGS)
library(plyr)
####van genutchen data####
#model directory
modDI <- c("c:\\Users\\hkropp\\Google Drive\\hydrus\\van_genut\\run7")
#model output
chain1 <- read.bugs(paste0(modDI,"\\CODAchain1.txt"))
chain2 <- read.bugs(paste0(modDI,"\\CODAchain2.txt"))
chain3 <- read.bugs(paste0(modDI,"\\CODAchain3.txt"))
chains <- rbind(chain1[[1]],chain2[[1]],chain3[[1]])	

#get parameter info
parmS <- data.frame(Mean= apply(chains,2,mean),	
			pc2.5=apply (chains, 2, quantile,probs=0.025)	 ,
			pc97.5=apply (chains, 2, quantile,probs=0.975))
			
parmG <- parmS[rownames(parmS)=="alpha.cm"|rownames(parmS)=="n",]		
#note in the model residual S was 0.01 since that was the lowest
#value in the dataset
#and qs is 0.41 since it is the same across all soil types

####soil texture parameters####

#read in soil texture data
datT <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\texture.csv")

#read in soil moisture curves
datM <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\soil_mrc.csv")
colnames(datM)[2] <- "Depth"
#read in hydrus soil catalog
#contains parameters from
#Carsel and Parrish [1988] for each soil texture classification
datH <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus_soil_cat.csv")

#change soil type name to match
colnames(datH)[1] <- "texture"

#aggregate by soil layer then join to hydrus layer



#join texture to characteristics for each shrub
#soilT <- join(datT,datH, by="texture", type="left")