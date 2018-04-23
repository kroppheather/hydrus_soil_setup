library(plyr)
library(R2OpenBUGS)
library(coda)

#set model directory
modD <- "c:\\Users\\hkropp\\Google Drive\\hydrus\\van_genut\\run8"

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


#look at soil across depth

soilT$depthI <- ifelse(soilT$Depth=="0-10",1,
				ifelse(soilT$Depth=="10-20",2,
				ifelse(soilT$Depth=="20-30",3,4)))
				
depthR <- aggregate(soilT$qr,by=list(soilT$depthI),FUN="mean")
colnames(depthR) <- c("depthI","qr")
depthS <- aggregate(soilT$qs,by=list(soilT$depthI),FUN="mean")
colnames(depthS) <- c("depthI","qs")

depthT <- aggregate(soilT$sand.,by=list(soilT$depthI),FUN="mean")
colnames(depthT) <- c("depthI","sand")
depthC <- aggregate(soilT$clay.,by=list(soilT$depthI),FUN="mean")
colnames(depthC) <- c("depthI","clay")			

depthK <- aggregate(soilT$Ks.cm.day,by=list(soilT$depthI),FUN="mean")
colnames(depthK) <- c("depthI","KS")	

theta <- function(psi.r,psi.s,alpha.mpa,h.mpa,n){


psi.r + ((psi.s-psi.r)/((1+((alpha.mpa*h.mpa)^n))^(1-(1/n))))
}	


# 1MPa =1019.7 cm of water




#convert MPa to cm of water

plot(soilAllx$vwc[soilAllx$depthI==1],abs(soilAllx$wp[soilAllx$depthI==1])*1019.7, pch=19,
		col="cornflowerblue", xlab="volumetric soil moisture", ylab ="water potential (- cm)")
points(soilAllx$vwc[soilAllx$depthI==2],abs(soilAllx$wp[soilAllx$depthI==2])*1019.7, pch=19,
		col="darkgreen")		
points(soilAllx$vwc[soilAllx$depthI==3],abs(soilAllx$wp[soilAllx$depthI==3])*1019.7, pch=19,
		col="tomato3")	
points(soilAllx$vwc[soilAllx$depthI==4],abs(soilAllx$wp[soilAllx$depthI==4])*1019.7, pch=19,
		col="darkorchid3")

points(theta(0.01,mean(soilT$qs),
		.2,seq(0,150000),1.35),	seq(0,150000),type="l")
#start by reading in psi as data based on texture and see how model runs
datalist <- list(Nobs=dim(soilAllx)[1],
				psi=soilAllx$vwc,
				psi.r=rep(0.01,4),
				psi.s=depthS$qs,
				h.cm=abs(soilAllx$wp)*1019.7,
				depth=soilAllx$depthI,
				Ndepth=4)
				
				
#starting values
startV <- list(list(n=rnorm(4,1.8,.1),alpha.cm=runif(4,.15,.20)),
				list(n=rnorm(4,2.2,.1),alpha.cm=runif(4,.10,.15)),
				list(n=rnorm(4,2.6,.1),alpha.cm=runif(4,.20,.25)))

params <- c("alpha.cm","n","sig.psi","rep.psi")				
				
bugs(data=datalist, inits=startV,parameters.to.save=params,
             n.iter=5000, n.chains=3, n.burnin=2000, n.thin=10,
             model.file="c:\\Users\\hkropp\\Documents\\GitHub\\hydrus_soil_setup\\van_genut.txt",
			 codaPkg=TRUE,
             OpenBUGS.pgm="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe",
			 debug=TRUE,
             working.directory=paste0(modD))				

			 
#start plot that reads in all parameters
#read in coda
chain1 <- read.bugs("c:\\Users\\hkropp\\Google Drive\\hydrus\\van_genut\\run3\\CODAchain1.txt")
chain2 <- read.bugs("c:\\Users\\hkropp\\Google Drive\\hydrus\\van_genut\\run3\\CODAchain2.txt")
chain3 <- read.bugs("c:\\Users\\hkropp\\Google Drive\\hydrus\\van_genut\\run3\\CODAchain3.txt")


chains <- rbind(chain1[[1]],chain2[[1]],chain3[[1]])	

parmMean <- apply(chains,2,mean)	
parm2.5 <- apply (chains, 2, quantile,probs=0.025)	 
parm97.5 <- apply (chains, 2, quantile,probs=0.975)

OUT <- data.frame(parms =colnames(chains), Mean=parmMean,pc2.5=parm2.5,
				pc97.5=parm97.5)
				
#evaluate model fit and parameters
for(i in 1:dim(ids)[1]){
	jpeg(paste0(modD,"\\plots\\soil\\",ids$ShrubID[i],"_",ids$depthI[i],".jpg"),
			quality=100,width=1000,height=1000)
		plot(abs(soilAllx$wp[soilAllx$shrubD==ids$shrubD[i]])*1019.7,
				soilAllx$vwc[soilAllx$shrubD==ids$shrubD[i]], pch=19,
				col="cornflowerblue",cex=2)
		points(seq(0,150000),theta(0.01,mean(soilT$qs),
		.14,seq(0,150000),1.3),	type="l", lwd=3, col="tomato3")	
		
	dev.off()		
}
		