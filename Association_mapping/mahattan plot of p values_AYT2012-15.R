setwd("/Users/xuehuili/Desktop/Xuehui/teaching/PLSC749 Applied plant molecular breeding/labs/association mapping")

####################################
#load AYT12-15 qualityT data
p.value.quality=read.delim("p values_AYT12-15Ind699_GBS50_qualityT_Nmodel.txt", header=T, na.string="na")
dim(p.value.quality)
p.value.quality[1:5,]

p.WT=p.value.quality[,c(1:3,4)]
colnames(p.WT)=c("SNP", "CHR", "BP", "P")
p.KWT=p.value.quality[,c(1:3,5)]
colnames(p.KWT)=c("SNP", "CHR", "BP", "P")
p.SED=p.value.quality[,c(1:3,6)]
colnames(p.SED)=c("SNP", "CHR", "BP", "P")
p.EXT=p.value.quality[,c(1:3,7)]
colnames(p.EXT)=c("SNP", "CHR", "BP", "P")
p.PRO=p.value.quality[,c(1:3,8)]
colnames(p.PRO)=c("SNP", "CHR", "BP", "P")
p.COL=p.value.quality[,c(1:3,9)]
colnames(p.COL)=c("SNP", "CHR", "BP", "P")


p.value.P.quality=read.delim("p values_AYT12-15Ind699_GBS50_qualityT_Pmodel.txt", header=T, na.string="na")
dim(p.value.P.quality)
p.value.P.quality[1:5,]

p.WT.P=p.value.P.quality[,c(1:3,4)]
colnames(p.WT.P)=c("SNP", "CHR", "BP", "P")
p.KWT.P=p.value.P.quality[,c(1:3,5)]
colnames(p.KWT.P)=c("SNP", "CHR", "BP", "P")
p.SED.P=p.value.P.quality[,c(1:3,6)]
colnames(p.SED.P)=c("SNP", "CHR", "BP", "P")
p.EXT.P=p.value.P.quality[,c(1:3,7)]
colnames(p.EXT.P)=c("SNP", "CHR", "BP", "P")
p.PRO.P=p.value.P.quality[,c(1:3,8)]
colnames(p.PRO.P)=c("SNP", "CHR", "BP", "P")
p.COL.P=p.value.P.quality[,c(1:3,9)]
colnames(p.COL.P)=c("SNP", "CHR", "BP", "P")


p.value.K.quality=read.delim("p values_AYT12-15Ind699_GBS50_qualityT_Kmodel.txt", header=T, na.string="na")
dim(p.value.K.quality)
p.value.K.quality[1:5,]

p.WT.K=p.value.K.quality[,c(1:3,4)]
colnames(p.WT.K)=c("SNP", "CHR", "BP", "P")
p.KWT.K=p.value.K.quality[,c(1:3,5)]
colnames(p.KWT.K)=c("SNP", "CHR", "BP", "P")
p.SED.K=p.value.K.quality[,c(1:3,6)]
colnames(p.SED.K)=c("SNP", "CHR", "BP", "P")
p.EXT.K=p.value.K.quality[,c(1:3,7)]
colnames(p.EXT.K)=c("SNP", "CHR", "BP", "P")
p.PRO.K=p.value.K.quality[,c(1:3,8)]
colnames(p.PRO.K)=c("SNP", "CHR", "BP", "P")
p.COL.K=p.value.K.quality[,c(1:3,9)]
colnames(p.COL.K)=c("SNP", "CHR", "BP", "P")


p.value.PK.quality=read.delim("p values_AYT12-15Ind699_GBS50_qualityT_PKmodel.txt", header=T, na.string="na")
dim(p.value.PK.quality)
p.value.PK.quality[1:5,]

p.WT.PK=p.value.PK.quality[,c(1:3,4)]
colnames(p.WT.PK)=c("SNP", "CHR", "BP", "P")
p.KWT.PK=p.value.PK.quality[,c(1:3,5)]
colnames(p.KWT.PK)=c("SNP", "CHR", "BP", "P")
p.SED.PK=p.value.PK.quality[,c(1:3,6)]
colnames(p.SED.PK)=c("SNP", "CHR", "BP", "P")
p.EXT.PK=p.value.PK.quality[,c(1:3,7)]
colnames(p.EXT.PK)=c("SNP", "CHR", "BP", "P")
p.PRO.PK=p.value.PK.quality[,c(1:3,8)]
colnames(p.PRO.PK)=c("SNP", "CHR", "BP", "P")
p.COL.PK=p.value.PK.quality[,c(1:3,9)]
colnames(p.COL.PK)=c("SNP", "CHR", "BP", "P")

## load package "qqman"
tiff(filename="Figure 1b. manhattan plots for AYT12-15 qualityT.tiff",width=10, height=8,units='in',res=300)
tiff(filename="Figure 1c. manhattan plots for AYT12-15 qualityT_Langdon.tiff",width=10, height=8,units='in',res=300)
tiff(filename="Figure 1d. manhattan plots for AYT12-15 qualityT_Williston.tiff",width=10, height=8,units='in',res=300)
par(mfrow=c(4,1))
par(mfrow=c(1,1), las=1)
par(mfrow=c(3,2), las=1)
manhattan(p.WT.PK, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Test weight")
manhattan(p.SED.PK, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Sedimentation volume")
manhattan(p.EXT.PK, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Semolina extraction rate")
manhattan(p.PRO.PK, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Semolina protein content")
manhattan(p.COL.PK, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Semolina color")
manhattan(p.KWT.PK, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Thousand kernel weight")
dev.off()

## estimate expected p values from observed p values
exp_pvalue.WT=-log(ppoints(length(p.WT$P)))
exp_pvalue.WT.K=-log(ppoints(length(p.WT.K$P)))
exp_pvalue.WT.P=-log(ppoints(length(p.WT.P$P)))
exp_pvalue.WT.PK=-log(ppoints(length(p.WT.PK$P)))

length(exp_pvalue.WT.K)
MSD.WT=sum((p.WT$P-exp_pvalue.WT)^2)/length(exp_pvalue.WT)
MSD.WT
MSD.WT.K=sum((p.WT.K$P-exp_pvalue.WT.K)^2)/length(exp_pvalue.WT.K)
MSD.WT.K
MSD.WT.Q=sum((p.WT.Q$P-exp_pvalue.WT.Q)^2)/length(exp_pvalue.WT.Q)
MSD.WT.Q
MSD.WT.QK=sum((p.WT.QK$P-exp_pvalue.WT.QK)^2)/length(exp_pvalue.WT.QK)
MSD.WT.QK

tiff(filename="Figure 4a. manhattan plots (p-value) for WT_Infiltration.mean_TASSEL.tiff",width=9, height=8,units='in',res=300)
par(mfrow=c(4,1))
manhattan(p.WT, main="WT (TASSEL:Naive) (MSD=1.68)", suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"))
manhattan(p.WT.K, main="WT (TASSEL:K) (MSD=1.33)", suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"))
manhattan(p.WT.P, main="WT (TASSEL:P) (MSD=1.60)", suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"))
manhattan(p.WT.PK, main="WT (TASSEL:P+K) (MSD=1.34)", suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"))
dev.off()


# estimate FDR
fdr.WT=p.adjust(p.value.quality[,4], method=c("fdr"))
fdr.WT=cbind(p.value.quality[1:3], fdr.WT)
colnames(fdr.WT)=c("SNP", "CHR", "BP", "P")
manhattan(fdr.WT, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"))

#load AYT12-15 agronomicT data
p.value.Agron=read.delim("p values_AYT12-15Ind699_GBS50_agronT_PK.txt", header=T, na.string="na")
dim(p.value.Agron)
p.value.Agron[1:5,]

p.GY=p.value.Agron[,c(1:3,4)]
colnames(p.GY)=c("SNP", "CHR", "BP", "P")
p.DTH=p.value.Agron[,c(1:3,5)]
colnames(p.DTH)=c("SNP", "CHR", "BP", "P")
p.HT=p.value.Agron[,c(1:3,6)]
colnames(p.HT)=c("SNP", "CHR", "BP", "P")
p.GY_L=p.value.Agron[,c(1:3,7)]
colnames(p.GY_L)=c("SNP", "CHR", "BP", "P")
p.GY_W=p.value.Agron[,c(1:3,9)]
colnames(p.GY_W)=c("SNP", "CHR", "BP", "P")
p.HT_L=p.value.Agron[,c(1:3,8)]
colnames(p.HT_L)=c("SNP", "CHR", "BP", "P")
p.HT_W=p.value.Agron[,c(1:3,10)]
colnames(p.HT_W)=c("SNP", "CHR", "BP", "P")

tiff(filename="Figure. manhattan plots for AYT12-15 agronomicT.tiff",width=9, height=8,units='in',res=300)
par(mfrow=c(4,1))
par(mfrow=c(1,1), las=1)
par(mfrow=c(2,2), las=1)
manhattan(p.GY, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Grain yield")
manhattan(p.DTH, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Date to heading")
manhattan(p.HT, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Plant height")
manhattan(p.GY_L, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Grain yield at Langdon")
manhattan(p.GY_W, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Grain yield at Williston")
manhattan(p.HT_L, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Grain yield at Langdon")
manhattan(p.HT_W, suggestiveline=F, genomewideline = F, col = c("blue4", "orange3"), chrlabs=c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"), main="Grain yield at Williston")
dev.off()
