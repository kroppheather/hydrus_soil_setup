library(plyr)
#read in soil texture data
datT <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\texture.csv")

#read in soil moisture curves
datM <- read.csv("c:\\Users\\hkropp\\Google Drive\\hydrus\\soil_mrc.csv")

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


