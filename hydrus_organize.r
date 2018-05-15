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

#read in Pima station radiation for gap fill
datFeb15 <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\Feb 2015.csv")						
datFeb15$doy <-	rep(seq(69,54, by=-1),each=49)					
datFeb15$year <- rep(2015, dim(datFeb15)[1])
datFeb15$hour <- rep(c(NA, seq(24,.5,by=-.5)),times=16)						
#read in Pima station april						
datApril16 <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\April 2016.csv")
datApril16$doy <-	rep(seq(214,102, by=-1),each=49)					
datApril16$year <- rep(2016, dim(datApril16)[1])
datApril16$hour <- rep(c(NA, seq(24,.5,by=-.5)),times=113)		
#read in Pima station august					
datAug16 <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\August 2016.csv")
datAug16$doy <-	rep(seq(254,239, by=-1),each=49)					
datAug16$year <- rep(2016, dim(datAug16)[1])
datAug16$hour <- rep(c(NA, seq(24,.5,by=-.5)),times=16)	
#read in south mountain station
datSept16 <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\sept 2016 south mountain.csv")
datSept16 $doy <-	rep(seq(251,250, by=-1),each=49)					
datSept16 $year <- rep(2016, dim(datSept16)[1])
datSept16 $hour <- rep(c(NA, seq(24,.5,by=-.5)),times=2)	
	

PimaAll <- rbind(datFeb15,datApril16,datAug16)
PimaAll <- PimaAll[,-1]
#now join
datMET<- join(datMET,PimaAll,by=c("doy","year","hour"),type="left")
#gapfill datMET

datMET$SR_gap <- ifelse(is.na(datMET$SR_4703_wsqm),
					datMET$pima_wmsq,datMET$SR_4703_wsqm)
					
#now join South mountain park
datMET<- join(datMET,datSept16,by=c("doy","year","hour"),type="left")					

datMET$SR_gap2 <- ifelse(is.na(datMET$SR_gap),
					datMET$south_mountain_wsqm,datMET$SR_gap)					

datMET[is.na(datMET$SR_gap2),]						
#convert precip to cm
datMET$ppt_4500_cm <-datMET$ppt_4500_in*2.54
#convert precip to mm
datMET$ppt_4500_mm <- datMET$ppt_4500_cm*10
#gapfill precipitation data with site
precipR <- lm(datMET$ppt_4500_mm~datMET$ppt_site_mm)						
plot(datMET$ppt_site_mm,datMET$ppt_4500_mm)						
#gapfill precip
datMET$ppt_4500_gap_mm <- ifelse(is.na(datMET$ppt_4500_mm),
					0+(precipR$coefficients[2]*datMET$ppt_site_mm),
					datMET$ppt_4500_mm)
								
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
datMET$PET.MS <- PET.MS(datMET$SR_gap2,datMET$lambda)
#cm/hr
plot(datMET$year+(datMET$doy/366)+(datMET$hour/24),datMET$PET.MS*100*60)

plot(datMET$year+(datMET$doy/366)+(datMET$hour/24),datMET$PER.WM)


#covert precip back to cm
#aggregate to hourly
hhR <- which(datMET$hour-floor(datMET$hour)==0)
precipHH <- datMET$ppt_4500_gap_mm[hhR]+datMET$ppt_4500_gap_mm[hhR-1]


datMET$PET.cmhr <- datMET$PET.MS*100*60*60



#put together data frame for hydrus
#Time=hr
#precip cm /hour
#potET cm /hour
#LAI assumed to be 0.64
#used hydrus default extinction 0.39
#need to calculate PET root vs soil for reading in manually
metOut <- data.frame(Time=seq(1,length(precipHH)),
			Precip=precipHH,	
			rsoil = round(datMET$PET.cmhr[hhR]*exp(-0.39*0.64),6),
			hCritA = rep(99000,length(precipHH)),
			rroot = round(datMET$PET.cmhr[hhR]*(1-exp(-0.39*0.64)),6))
			
			
			
write.table(metOut, "c:\\Users\\hkropp\\Google Drive\\hydrus\\met_output.csv",
			sep=",",row.names=FALSE)	


write.table(datMET, "c:\\Users\\hkropp\\Google Drive\\hydrus\\met_All_file.csv",
			sep=",",row.names=FALSE)	



#########################################
#read in root data

rootDat <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\root extract.csv")

rcounts <- rootDat$depth.end-rootDat$depth.start	
rootInc <- data.frame(depthI = seq(0,80),
			rootP=c(0,rep(rootDat$r.perc,times=rcounts)))
			
sum(rootInc$rootP)			
write.table(rootInc,"c:\\Users\\hkropp\\Google Drive\\hydrus\\rootAll_out.csv",
			sep=",", row.names=FALSE)		