#######
setwd("/Users/desert/OneDrive - North Dakota University System/Research/ProteinAndTanSpot")

pvalue = read.table("PC5_K_statistics.txt",header = T, na.strings = "NaN")
pvalue[1:5,1:5]

ToxCqValue= p.adjust(pvalue$p[pvalue$Trait=="ToxC"], method = "fdr")
ND12qValue= p.adjust(pvalue$p[pvalue$Trait=="ND12"], method = "fdr")
SedVblupqValue= p.adjust(pvalue$p[pvalue$Trait=="SedVblup"], method = "fdr")

qvalue = c(ToxCqValue,SedVblupqValue,ND12qValue)

pvalue$q = qvalue

write.csv(pvalue,"pvalue.csv")
