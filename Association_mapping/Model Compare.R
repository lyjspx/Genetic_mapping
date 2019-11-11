getwd()
setwd("/Users/desert/OneDrive - North Dakota University System/Research/ProteinAndTanSpot")

allFile=list.files()

naiveModel=read.table(allFile[9],header = T,na.strings = "NaN")
naiveModel[1:5,1:5]
pModel=read.table(allFile[16],header = T,na.strings = "NaN")
kModel=read.table(allFile[5],header = T,na.strings = "NaN")
pkModel=read.table(allFile[15],header = T,na.strings = "NaN")

expPValue=ppoints(39483/3)

library(hydroGOF)
MsdSedVblupNaive=mse(sort(naiveModel[naiveModel$Trait=="SedVblup",]$p),expPValue)
MsdSedVblupP=mse(sort(pModel[pModel$Trait=="SedVblup",]$p),expPValue)
MsdSedVblupK=mse(sort(kModel[kModel$Trait=="SedVblup",]$p),expPValue)
MsdSedVblupPK=mse(sort(pkModel[pkModel$Trait=="SedVblup",]$p),expPValue)


