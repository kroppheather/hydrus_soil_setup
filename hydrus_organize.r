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
datT$depthI <- ifelse(datT$Depth=="0-10",1,
				ifelse(datT$Depth=="10-20",2,
				ifelse(datT$Depth=="20-30",3,4)))
				
				
#read in hydrus soil catalog
#contains parameters from
#Carsel and Parrish [1988] for each soil texture classification
datH <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus_soil_cat.csv")

#change soil type name to match
colnames(datH)[1] <- "texture"

#aggregate by soil layer then join to hydrus layer
sandL <- aggregate(datT$sand.,by=list(datT$depthI),FUN="mean")
clayL <- aggregate(datT$clay.,by=list(datT$depthI),FUN="mean")
colnames(sandL) <- c("depthI","sand")
colnames(clayL) <- c("depthI","clay")

textI <- data.frame(depthI=seq(1,4),
			texture=c("loamy sand","loamy sand","sandy loam","sandy loam"))

#join texture to characteristics for each shrub
textI <- join(textI,datH, by="texture", type="left")

textI$Ks.cm.min <- (textI$Ks.cm.day/24)/60

####met data for ET ####
datMET <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\MCFD_met.csv")

