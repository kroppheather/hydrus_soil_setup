library(plyr)
library(R2OpenBUGS)
library(coda)

#set model directory
modD <- "c:\\Users\\hkropp\\Google Drive\\hydrus\\van_genut\\run1"

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

#join texture to characteristics for each shrub
soilT <- join(datT,datH, by="texture", type="left")

#conversion of megapascals to cm
# mPA * (cm/MPa)
#1MPa= 10197.1621 cm

#create a shrub id
#depth id
ids <- unique(data.frame(ShrubID=datT$ShrubID, Depth=datT$Depth))

ids$shrubD <- seq(1,dim(ids)[1])
ids$shrubI <- ifelse(ids$ShrubID=="1a",1,
				ifelse(ids$ShrubID=="2a",2,
				ifelse(ids$ShrubID=="3a",3,
				ifelse(ids$ShrubID=="4a",4,
				ifelse(ids$ShrubID=="5a",5,
				ifelse(ids$ShrubID=="6a",6,
				ifelse(ids$ShrubID=="1b",7,
				ifelse(ids$ShrubID=="2b",8,
				ifelse(ids$ShrubID=="3b",9,
				ifelse(ids$ShrubID=="4b",10,
				ifelse(ids$ShrubID=="5b",11,12)))))))))))
ids$depthI <- ifelse(ids$Depth=="0-10",1,
				ifelse(ids$Depth=="10-20",2,
				ifelse(ids$Depth=="20-30",3,4)))	

#join ids to data frames

soilTx <- join(soilT,ids,by=c("Depth","ShrubID"),type="left")
soilMx <- join(datM, ids,by=c("Depth","ShrubID"),type="left")				

#commbine texture data into soil Mx to get bulk density calculations
soilAllx <- join(soilMx,soilTx, by=c("Depth","ShrubID","shrubD","shrubI","depthI"),
				type="left")
#bd from saxton and rawls				
soilAllx$bd <- (1-soilAllx$qs)*2.65

#calculate vwc
soilAllx$vwc <- soilAllx$bd*soilAllx$rwc

#start by reading in psi as data based on texture and see how model runs
datalist <- list(Nobs=dim(soilAllx)[1],
				psi=soilAllx$vwc,
				psi.r=soilTx$qr,
				psi.s=soilTx$qs,
				h.mpa=abs(soilAllx$wp),
				shrubD=soilAllx$shrubD,
				NshrubD=dim(ids)[1])
				
				
#starting values
startV <- list(list(n=rnorm(dim(ids)[1],1.8,.1),alpha.cm=runif(dim(ids)[1],.15,.2)),
				list(n=rnorm(dim(ids)[1],2.2,.1),alpha.cm=runif(dim(ids)[1],.1,.15)),
				list(n=rnorm(dim(ids)[1],2.6,.1),alpha.cm=runif(dim(ids)[1],.01,.1)))

params <- c("alpha.mpa","alpha.cm","n","sig.psi")				
				
bugs(data=datalist, inits=startV,parameters.to.save=params,
             n.iter=5000, n.chains=3, n.burnin=2000, n.thin=10,
             model.file="c:\\Users\\hkropp\\Documents\\GitHub\\hydrus_soil_setup\\van_genut.txt",
			 codaPkg=TRUE,
             OpenBUGS.pgm="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe",
			 debug=TRUE,
             working.directory=paste0(modD))				
