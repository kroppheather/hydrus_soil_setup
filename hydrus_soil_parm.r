
#read in soil texture data
datT <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\texture.csv")

#read in soil moisture curves
datM <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\soil_mrc.csv")

#read in hydrus soil catalog
#contains parameters from
#Carsel and Parrish [1988] for each soil texture classification
datH <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\hydrus_soil_cat.csv")


#first get texture for values out of soil hydraulic table
#see if there is any significant difference in soil between shrubs and depths
#check how many groups
unique(data.frame(depth=datT$Depth, shrub=datT$ShrubID))

modT <- lm(datT$sand.~datT$ShrubID+datT$Depth)
anova(modT)

SandMean <- aggregate(datT$sand., by=list(datT$ShrubID), FUN="mean")
SandSD <- aggregate(datT$sand., by=list(datT$ShrubID), FUN="sd")
plot(SandMean$Group.1, SandMean$x, pch=19, ylim=c(70,90))
arrows(seq(1,12), SandMean$x-SandSD$x,seq(1,12), SandMean$x+SandSD$x, code=0)
abline(h=mean(datT$sand.))

modTC <- lm(datT$clay.~datT$ShrubID+datT$Depth)
anova(modTC)

ClayMean <- aggregate(datT$clay., by=list(datT$ShrubID), FUN="mean")
ClaySD <- aggregate(datT$clay., by=list(datT$ShrubID), FUN="sd")
plot(ClayMean$Group.1, ClayMean$x, pch=19, ylim=c(0,20))
arrows(seq(1,12), ClayMean$x-ClaySD$x,seq(1,12), ClayMean$x+ClaySD$x, code=0)
abline(h=mean(datT$clay.))

ClayD <- aggregate(datT$clay., by=list(datT$Depth), FUN="mean")
ClaySDD <- aggregate(datT$clay., by=list(datT$Depth), FUN="sd")


plot(ClayD$Group.1, ClayD$x, pch=19, ylim=c(0,20))
arrows(seq(1,4), ClayD$x-ClaySDD$x,seq(1,4), ClayD$x+ClaySDD$x, code=0)
abline(h=mean(datT$clay.))







