#read in hydrus output

hydOut <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus\\runs\\run_final\\all_out\\Obs_Node_out.csv")

datMET <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\met_All_file.csv")
hhR <- which(datMET$hour-floor(datMET$hour)==0)
thetaO <- data.frame(depth0=hydOut[,3])
for(i in 2:61){
	thetaO <- cbind(thetaO, hydOut[,i*3])
	colnames(thetaO)[i] <- paste0("depth",i-1)
}

#add time info

swcOut <- data.frame(doy=datMET$doy[hhR],year=datMET$year[hhR],hour=datMET$hour[hhR],
			thetaO)
			
plot(swcOut$depth30, type="l", col="cornflowerblue")
plot(swcOut$depth60[swcOut$year==2015], type="l", col="cornflowerblue")

write.table(swcOut,
			"c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus\\runs\\run_final\\all_out\\VWC_soil_out.csv",
			sep=",", row.names=FALSE)