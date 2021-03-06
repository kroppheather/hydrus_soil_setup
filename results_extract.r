#read in hydrus output

hydOut <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus\\runs\\run_final\\all_extend_hcrit_out\\Obs_Node_out.csv")
HHmod <- read.table("c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus\\runs\\run_final\\all_extend_hcrit_out\\Obs_Node_out.txt", sep = "")
datMET <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\met_All_file.csv")
hhR <- which(datMET$hour-floor(datMET$hour)==0)
thetaO <- data.frame(depth0=hydOut[,3])
for(i in 2:61){
	thetaO <- cbind(thetaO, hydOut[,i*3])
	colnames(thetaO)[i] <- paste0("depth",i-1)
}
hhO <- data.frame(depth0=HHmod[,2])
hhID <- seq(5,181,by=3)
for(i in 1:length(hhID)){
	hhO <- cbind(hhO, HHmod[,hhID[i]])
	colnames(hhO)[i+1] <- paste0("depth",i+1)
}

#add time info

swcOut <- data.frame(doy=datMET$doy[hhR],year=datMET$year[hhR],hour=datMET$hour[hhR],
			thetaO)
			
hhOut <-data.frame(doy=datMET$doy[hhR],year=datMET$year[hhR],hour=datMET$hour[hhR],		
			hhO)
plot(swcOut$depth30, type="l", col="cornflowerblue")
plot(rowMeans(swcOut[,4:34]), type="l", col="cornflowerblue")
plot(rowMeans(swcOut[swcOut$year>=2015,4:34]), type="l", col="cornflowerblue")
plot(swcOut$depth60, type="l", col="cornflowerblue")
plot(swcOut$depth60[swcOut$year>=2015], type="l", col="cornflowerblue")

plot(rowMeans(hhOut[,4:34])/10197, type="l", col="cornflowerblue")
plot(hhOut[,31]/10197, type="l", col="cornflowerblue")
write.table(swcOut,
			"c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus\\runs\\run_final\\all_extend_hcrit_out\\VWC_soil_out.csv",
			sep=",", row.names=FALSE)
			
write.table(hhOut,
			"c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus\\runs\\run_final\\all_extend_hcrit_out\\HH_soil_out.csv",
			sep=",", row.names=FALSE)			