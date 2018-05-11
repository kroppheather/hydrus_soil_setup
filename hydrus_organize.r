######################################
#####organizes data for model run#####
######################################
library(R2OpenBUGS)
library(plyr)
library(lubridate)
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
datMET <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\MCFD_cleaned_all.csv")

#extract date information
dateMET <- as.Date(datMET$Date, "%m/%d/%Y")
datMET$doy <- yday(dateMET)
datMET$year <- year(dateMET)

#convert time to hour
datMET$hour <- round_any(datMET$Time/100,.5)

#gap fill missing data
tempR <- lm(datMET$temp_4507_C~datMET$temp_4617_C)

datMET$temp4507_gap <- ifelse(is.na(datMET$temp_4507_C),
						tempR$coefficients[1]+(tempR$coefficients[2]*datMET$temp_4617_C),
						datMET$temp_4507_C)



#calculate priestly taylor
#esat
#ta degrees C
e.sat <- function(Ta){
	0.61121*exp((17.502*Ta)/(Ta+240.97))

}


#Ta degrees C and e.sat kpa
DELTA <- function(Ta,e.sat){
	(17.502*240.97*e.sat)/((Ta+240.97)^2)

}

#PET calc
#Rn W/m2

PET.WM <- function(Rn,DELTA){
	1.26*(DELTA/(DELTA+0.066))*Rn

}

#latent heat of vaporization
lambda <- function(Ta){
	2495000-((0.00236*Ta)*1000000)

}


#convert PET from W/m2 to M/s
#by dividing by the density of water and the latent heat of vaporization
PET.MS <- function(PET.WM, lambda){
	PET.WM/(1000*lambda)

}


#temp_4507_C
#SR_4703_wsqm

#calculations
datMET$e.sat <- e.sat(datMET$temp4507_gap)
datMET$DELTA <- DELTA(datMET$temp4507_gap,datMET$e.sat)
datMET$PER.WM <- PET.WM(datMET$SR_4703_wsqm,datMET$DELTA)
datMET$lambda <- lambda(datMET$temp4507_gap)
datMET$PET.MS <- PET.MS(datMET$PER.WM,datMET$lambda)
#mm/hr
plot(datMET$year+(datMET$doy/366)+(datMET$hour/24),datMET$PET.MS*100*60)

plot(datMET$year+(datMET$doy/366)+(datMET$hour/24),datMET$PER.WM)