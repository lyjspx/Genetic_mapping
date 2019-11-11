#Set working directory under Mac
setwd("/Users/desert/Google Drive/CodeCenter/rdata")

#load packages "qtl"
DM3=read.cross("csvr", "/Users/desert/Google Drive/Alfalfa/Alfalfa_tetra_map", "AH_map_pheno_r.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "bc")

#trait1 = matDupy12OCT
#trait2 = matPepsi13July
#trait3 = matPepsi14Sep
#trait4 = matPepsi14Nov

str(DM3)
nind(DM3)
nchr(DM3)
totmar(DM3)
nmar(DM3)
plot(DM3)
plotMap(DM3)
plotMissing(DM3)

##########################################################
#Spearman's correlations among the phenotypes
DM3.pheno = DM3$pheno[complete.cases(DM3$pheno[,2:5]),]
str(DM3.pheno)
mean(DM3.pheno[,2])
min(DM3.pheno[,2])
max(DM3.pheno[,2])

mean(DM3.pheno[,3])
min(DM3.pheno[,3])
max(DM3.pheno[,3])

mean(DM3.pheno[,4])
min(DM3.pheno[,4])
max(DM3.pheno[,4])

mean(DM3.pheno[,5])
min(DM3.pheno[,5])
max(DM3.pheno[,5])

cor(DM3.pheno[,2], DM3.pheno[,3])
cor.test(DM3.pheno[,2], DM3.pheno[,3], method="spearman")
cor.test(DM3.pheno[,2], DM3.pheno[,3], method="pearson")

cor(DM3.pheno[,2], DM3.pheno[,4])
cor.test(DM3.pheno[,2], DM3.pheno[,4], method="spearman")
cor.test(DM3.pheno[,2], DM3.pheno[,4], method="pearson")

cor(DM3.pheno[,2], DM3.pheno[,5])
cor.test(DM3.pheno[,2], DM3.pheno[,5], method="spearman")
cor.test(DM3.pheno[,2], DM3.pheno[,5], method="pearson")

cor(DM3.pheno[,3], DM3.pheno[,4])
cor.test(DM3.pheno[,3], DM3.pheno[,4], method="spearman")
cor.test(DM3.pheno[,3], DM3.pheno[,4], method="pearson")

cor(DM3.pheno[,3], DM3.pheno[,5])
cor.test(DM3.pheno[,3], DM3.pheno[,5], method="spearman")
cor.test(DM3.pheno[,3], DM3.pheno[,5], method="pearson")

cor(DM3.pheno[,4], DM3.pheno[,5])
cor.test(DM3.pheno[,4], DM3.pheno[,5], method="spearman")
cor.test(DM3.pheno[,4], DM3.pheno[,5], method="pearson")

boxplot(DM3.pheno[,2:5], ylab="Score")



##########################################################
#remove samples with a lot of missing marker values

nind(DM3)
par(mfrow=c(1,1), las=1)
plot(ntyped(DM3), xlab="Sample", ylab="No. typed markers", main="No. typed markers by samples") 
DM3 <- subset(DM3, ind=(ntyped(DM3)>(totmar(DM3)/2)))
nind(DM3)

#identify duplicate individuals
# cg <- comparegeno(DM3)
# str(cg)
# hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes") 
# rug(cg[lower.tri(cg)])
# wh <- which(cg > 0.9, arr=TRUE)
# wh <- wh[wh[,1] < wh[,2],]
# wh
# genotype=pull.geno(DM3)
# str(genotype)
# dim(genotype)
# table(genotype[70,], genotype[81,])
# 
# gc=genClones(DM3, tol=0.8)
# gc$cgd
# DM3=fixClones(DM3, gc$cgd, consensus = TRUE)
# nind(DM3)
#################
#plot genetic map of DM3
tiff(filename="Fig. Genetic map of DM3.tiff", width=8, height=6,units='in',res=300)
plotMap(DM3, main="Genetic map of DM3" )
dev.off()

write.table(summaryMap(DM3),file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")


###################################################################
#single QTL genome scan
#calculate conditional probability if using "imp" method
DM3.prob <- calc.genoprob(DM3, step=1, error.prob=0.01)
str(DM3.prob)

#simulate missing genotype if using "imp" method
DM3=jittermap(DM3, amount=1e-6)
DM3.sim <- sim.geno(DM3, step=2, n.draws=128, err=0.001)
str(DM3.sim)

##1.1 single-QTL genome scan for trait1
##scanone wiht method "imp"
out.imp.trait1 <- scanone(DM3.sim, method="imp", pheno.col=2)
summary(out.imp.trait1,threshold = 3.23)
plot(out.imp.trait1)

operm.imp.trait1 <- scanone(DM3.sim, pheno.col=2, method="imp", n.perm=1000)
summary(operm.imp.trait1, alpha=0.05) #3.23
summary(operm.imp.trait1, alpha=0.01) #4.03

plot(out.imp.trait1,main = "QTL mapping of matDupy12OCT in DM3")
abline(h=3.23, col=3, lty=2)
text(750,3.5,"1000 permutation at p=0.05, threshold = 3.23",col = 3)

lodint(out.imp.trait1,chr="3D",drop=1.5)
lodint(out.imp.trait1,chr="7B",drop=1.5)
lodint(out.imp.trait1,chr="7D",drop=1.5)
drop1.5 = rbind(lodint(out.imp.trait1,chr="3D",drop=1.5),
                lodint(out.imp.trait1,chr="7B",drop=1.5),
                lodint(out.imp.trait1,chr="7D",drop=1.5))

write.table(drop1.5,file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")


chr.trait1 <- c("3D","7B","7D")
pos.trait1 <- c(66.8,47.1,48.4)
qtl.trait1 <- makeqtl(DM3.sim, chr.trait1, pos.trait1)
my.formula=y~Q1
out.fq.trait1 <- fitqtl(DM3.sim, pheno.col=2, qtl=qtl.trait1, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.trait1



# 
# #two-QTL genome scan
# ##1.2 two-QTL genome scan with "scantwo" for trait1
# out2.imp.DM3trait1 <- scantwo(DM3.prob, pheno.col=2, method="imp")
# operm2.imp.DM3trait1 <- scantwo(DM3.prob, pheno.col=2, method="imp", n.perm=100)
# str(operm2.imp.DM3trait1)
# summary(operm2.imp.DM3trait1, alpha=0.01) #7.36 5.78 5.21 6.23 3.86 3.83
# summary(operm2.imp.DM3trait1, alpha=0.05) #7.08 5.47 4.8 5.52 3.56 3.54
# summary(out2.imp.DM3trait1, thresholds=c(7.08, 5.47, 4.8, 5.52, 3.56), what="int")
# 
# 
# #stepwise selection of multiple-QTL model
# ##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
# calc.penalties(operm2.imp.DM3trait1, alpha=0.05)
# 
# ##1.4 stepwise selection for trait1
# chr.trait1 <- c("S3D")
# pos.trait1 <- c(66.76)
# qtl.trait1 <- makeqtl(DM3.sim, chr.trait1, pos.trait1)
# qtl.trait1
# my.formula=y~Q1
# stepout.DM3trait1=stepwiseqtl(DM3.sim, pheno.col=2, additive.only=FALSE, penalties=c(3.55, 4.79, 1.92), qtl = qtl.trait1, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
# stepout.DM3trait1
# 
# chr.trait1 <- c("S3D")
# pos.trait1 <- c(66.76)
# qtl.trait1 <- makeqtl(DM3.sim, chr.trait1, pos.trait1)
# my.formula=y~Q1
# out.fq.trait1 <- fitqtl(DM3.sim, pheno.col=5, qtl=qtl.trait1, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
# out.fq.trait1
# 
# 
# 
# 
# stepout.DM3trait1=stepwiseqtl(DM3.prob, pheno.col=2, additive.only=FALSE, penalties=c(3.73, 5.1, 1.93), qtl = qtl.trait1, verbose=FALSE, method="imp", keeptrace = TRUE)
# stepout.DM3trait1
# 
# names(attributes(stepout.DM3trait1))
# trait1.trace=attr(stepout.DM3trait1, "trace")
# par(mfrow=c(3,4))
# for(i in seq(along=trait1.trace))
#   plotModel(trait1.trace[[i]], chronly=TRUE,
#             main=paste(i, ": pLOD =",
#                        round(attr(trait1.trace[[i]], "pLOD"), 2)))
# 
# chr.trait1 <- c("1A", "1B", "5A", "7A")
# pos.trait1 <- c(144, 52.876, 176, 76)
# qtl.trait1 <- makeqtl(DM3.sim, chr.trait1, pos.trait1)
# my.formula=y~Q1 + Q2 + Q3 + Q4 +Q1:Q4
# out.fq.trait1 <- fitqtl(DM3.sim, pheno.col=2, qtl=qtl.trait1, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
# out.fq.trait1
# 
# rqtl.trait1 =refineqtl(DM3.sim, pheno.col=2, qtl=qtl.trait1, method="imp", formula=my.formula)
# rqtl.trait1
# 
# tiff(filename="Fig. QTL for trait1.tiff", width=8, height=6,units='in',res=300)
# plot(out.imp.trait1, chr=c("1A", "1B", "5A", "7A"), col=4, ylim=c(0,85), ylab="LOD", lwd=1)
# plotLodProfile(rqtl.trait1, col=2, add=T, qtl.labels = FALSE, lwd=1)
# abline(h=4.58, col=3, lty=2)
# dev.off()
# 
# mar.1A.trait1=find.marker(DM3.sim, chr="1A",pos=144)
# mar.1A.trait1
# mar.1B.trait1=find.marker(DM3.sim, chr="1B",pos=52.876)
# mar.1B.trait1
# mar.5A.trait1=find.marker(DM3.sim, chr="5A",pos=176)
# mar.5A.trait1
# mar.7A.trait1=find.marker(DM3.sim, chr="7A",pos=76)
# mar.7A.trait1
# 
# bayesint(rqtl.trait1, qtl.index=1, prob=0.99)
# bayesint(rqtl.trait1, qtl.index=2, prob=0.99)
# bayesint(rqtl.trait1, qtl.index=3, prob=0.99)
# bayesint(rqtl.trait1, qtl.index=4, prob=0.99)
# 
# lodint(rqtl.trait1,qtl.index=1, drop=1.5)
# lodint(rqtl.trait1,qtl.index=2, drop=1.5)
# lodint(rqtl.trait1,qtl.index=3, drop=1.5)
# lodint(rqtl.trait1,qtl.index=4, drop=1.5)
# 
# out.fq.trait1 <- fitqtl(DM3.sim, pheno.col=2, qtl=rqtl.trait1, method="imp", formula=my.formula, get.ests=TRUE)
# out.fq.trait1
# 
# tiff(filename="Fig. QTL interaction for trait1.tiff", width=6, height=6,units='in',res=300)
# effectplot(DM3.sim, pheno.col=2, mname1=mar.1A.trait1, geno1=c("Rusty", "PI466979"), ylim=c(2,7), mname2=mar.7A.trait1, geno2=c("Rusty", "PI466979"), xlab="7A@76.0", ylab="Resistance Score", legend.lab="1A@144.0")
# dev.off()


###########################################
####Trait2#################
##scanone wiht method "imp"
out.imp.trait2 <- scanone(DM3.sim, method="imp", pheno.col=3)
summary(out.imp.trait2,threshold = 3.27)
plot(out.imp.trait2)

operm.imp.trait2 <- scanone(DM3.sim, pheno.col=2, method="imp", n.perm=1000)
summary(operm.imp.trait2, alpha=0.05) #3.27
summary(operm.imp.trait2, alpha=0.01) 

plot(out.imp.trait2,main = "QTL mapping of matPepsi13July in DM3")
abline(h=3.27, col=3, lty=2)
text(750,3.5,"1000 permutation at p=0.05, threshold = 3.27",col = 3)


lodint(out.imp.trait2,chr="7B",drop=1.5)
lodint(out.imp.trait2,chr="8A",drop=1.5)
drop1.5 = rbind(lodint(out.imp.trait2,chr="8A",drop=1.5),
                lodint(out.imp.trait2,chr="7B",drop=1.5))

write.table(drop1.5,file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")


chr.trait2 <- c("7B","8A")
pos.trait2 <- c(32.5,50.3)
qtl.trait2 <- makeqtl(DM3.sim, chr.trait2, pos.trait2)
my.formula=y~Q1
out.fq.trait2 <- fitqtl(DM3.sim, pheno.col=3, qtl=qtl.trait2, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.trait2



# #two-QTL genome scan
# ##1.2 two-QTL genome scan with "scantwo" for trait2
# out2.imp.DM3trait2 <- scantwo(DM3.prob, pheno.col=3, method="imp")
# operm2.imp.DM3trait2 <- scantwo(DM3.prob, pheno.col=3, method="imp", n.perm=100)
# str(operm2.imp.DM3trait2)
# summary(operm2.imp.DM3trait2, alpha=0.01) #1% 12.2 11.2 9.15 8.34 6.98 4.11
# summary(operm2.imp.DM3trait2, alpha=0.05) #5% 11.2 10.4 8.3 7.43 6.32 3.88
# summary(out2.imp.DM3trait2, thresholds=c(11.2, 10.4, 8.3, 7.43, 6.32), what="int")
# 
# #stepwise selection of multiple-QTL model
# ##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
# calc.penalties(operm2.imp.DM3trait2, alpha=0.05)
# 
# ##1.4 stepwise selection for trait2
#  
# 
# chr.trait2 <- c("S7B")
# pos.trait2 <- c(32.547)
# qtl.trait2 <- makeqtl(DM3.sim, chr.trait2, pos.trait2)
# my.formula=y~Q1
# out.fq.trait2 <- fitqtl(DM3.sim, pheno.col=3, qtl=qtl.trait2, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
# out.fq.trait2
# 
# 
# names(attributes(stepout.DM3trait2))
# trait2.trace=attr(stepout.DM3trait2, "trace")
# par(mfrow=c(3,4))
# for(i in seq(along=trait2.trace))
#   plotModel(trait2.trace[[i]], chronly=TRUE,
#             main=paste(i, ": pLOD =",
#                        round(attr(trait2.trace[[i]], "pLOD"), 2)))
# 
# rqtl.trait2 =refineqtl(DM3.sim, pheno.col=3, qtl=qtl.trait2, method="imp", formula=my.formula)
# rqtl.trait2
# 
# tiff(filename="Fig. QTL for trait2.tiff", width=8, height=6,units='in',res=300)
# plot(out.imp.trait2, chr=c("1A", "1B", "5A", "7A"), col=4, ylim=c(0,85), ylab="LOD", lwd=1)
# plotLodProfile(rqtl.trait2, col=2, add=T, qtl.labels = FALSE, lwd=1)
# abline(h=4.58, col=3, lty=2)
# dev.off()
# 
# mar.S7B.trait2=find.marker(DM3.sim, chr="S7B",pos=32.547)
# mar.S7B.trait2
# 
# bayesint(rqtl.trait2, qtl.index=1, prob=0.99)
# 
# lodint(rqtl.trait2,qtl.index=1, drop=1.5)
# 
# out.fq.trait2 <- fitqtl(DM3.sim, pheno.col=3, qtl=rqtl.trait2, method="imp", formula=my.formula, get.ests=TRUE)
# out.fq.trait2
# 
# # tiff(filename="Fig. QTL interaction for trait2.tiff", width=6, height=6,units='in',res=300)
# # effectplot(DM3.sim, pheno.col=3, mname1=mar.1A.trait2, geno1=c("Rusty", "PI466979"), ylim=c(2,7), mname2=mar.7A.trait2, geno2=c("Rusty", "PI466979"), xlab="7A@76.0", ylab="Resistance Score", legend.lab="1A@144.0")
# # dev.off()

###########################################
####Trait3#################
##scanone wiht method "imp"
out.imp.trait3 <- scanone(DM3.sim, method="imp", pheno.col=4)
summary(out.imp.trait3,threshold = 3.33)
plot(out.imp.trait3)

operm.imp.trait3 <- scanone(DM3.sim, pheno.col=4, method="imp", n.perm=1000)
summary(operm.imp.trait3, alpha=0.05) #3.33
summary(operm.imp.trait3, alpha=0.01) #4.03
#No QTL

#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for trait3
out2.imp.DM3trait3 <- scantwo(DM3.prob, pheno.col=3, method="imp")
operm2.imp.DM3trait3 <- scantwo(DM3.prob, pheno.col=3, method="imp", n.perm=100)
str(operm2.imp.DM3trait3)
summary(operm2.imp.DM3trait3, alpha=0.01) #1% 12.2 11.2 9.15 8.34 6.98 4.11
summary(operm2.imp.DM3trait3, alpha=0.05) #5% 11.2 10.4 8.3 7.43 6.32 3.88
summary(out2.imp.DM3trait3, thresholds=c(11.2, 10.4, 8.3, 7.43, 6.32), what="int")

#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.DM3trait3, alpha=0.05)

##1.4 stepwise selection for trait3


chr.trait3 <- c("S7B")
pos.trait3 <- c(32.547)
qtl.trait3 <- makeqtl(DM3.sim, chr.trait3, pos.trait3)
my.formula=y~Q1
out.fq.trait3 <- fitqtl(DM3.sim, pheno.col=3, qtl=qtl.trait3, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.trait3


names(attributes(stepout.DM3trait3))
trait3.trace=attr(stepout.DM3trait3, "trace")
par(mfrow=c(3,4))
for(i in seq(along=trait3.trace))
  plotModel(trait3.trace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =",
                       round(attr(trait3.trace[[i]], "pLOD"), 2)))

rqtl.trait3 =refineqtl(DM3.sim, pheno.col=3, qtl=qtl.trait3, method="imp", formula=my.formula)
rqtl.trait3

tiff(filename="Fig. QTL for trait3.tiff", width=8, height=6,units='in',res=300)
plot(out.imp.trait3, chr=c("1A", "1B", "5A", "7A"), col=4, ylim=c(0,85), ylab="LOD", lwd=1)
plotLodProfile(rqtl.trait3, col=2, add=T, qtl.labels = FALSE, lwd=1)
abline(h=4.58, col=3, lty=2)
dev.off()

mar.S7B.trait3=find.marker(DM3.sim, chr="S7B",pos=32.547)
mar.S7B.trait3

bayesint(rqtl.trait3, qtl.index=1, prob=0.99)

lodint(rqtl.trait3,qtl.index=1, drop=1.5)

# tiff(filename="Fig. QTL interaction for trait3.tiff", width=6, height=6,units='in',res=300)
# effectplot(DM3.sim, pheno.col=3, mname1=mar.1A.trait3, geno1=c("Rusty", "PI466979"), ylim=c(2,7), mname2=mar.7A.trait3, geno2=c("Rusty", "PI466979"), xlab="7A@76.0", ylab="Resistance Score", legend.lab="1A@144.0")
# dev.off()

###########################################
####trait4#################
##scanone wiht method "imp"
out.imp.trait4 <- scanone(DM3.sim, method="imp", pheno.col=5)
summary(out.imp.trait4,threshold = 3.27)
plot(out.imp.trait4)

operm.imp.trait4 <- scanone(DM3.sim, pheno.col=5, method="imp", n.perm=1000)
summary(operm.imp.trait4, alpha=0.05) #3.27
summary(operm.imp.trait4, alpha=0.01) 

summary(out.imp.trait4,threshold=c(3.27))

plot(out.imp.trait4,main = "QTL mapping of matPepsi14Nov in DM3")
abline(h=3.27, col=3, lty=2)
text(750,3.5,"1000 permutation at p=0.05, threshold = 3.27",col = 3)

lodint(out.imp.trait4,chr="3D",drop=1.5)
lodint(out.imp.trait4,chr="7B",drop=1.5)
lodint(out.imp.trait4,chr="7D",drop=1.5)
drop1.5 = rbind(lodint(out.imp.trait4,chr="3D",drop=1.5),
                lodint(out.imp.trait4,chr="7B",drop=1.5),
                lodint(out.imp.trait4,chr="7D",drop=1.5))

write.table(drop1.5,file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")


chr.trait4 <- c("3D","7B","7D")
pos.trait4 <- c(69.4,53.8,50.4)
qtl.trait4 <- makeqtl(DM3.sim, chr.trait4, pos.trait4)
my.formula=y~Q3
out.fq.trait4 <- fitqtl(DM3.sim, pheno.col=5, qtl=qtl.trait4, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.trait4

# #two-QTL genome scan
# ##1.2 two-QTL genome scan with "scantwo" for trait4
# out2.imp.DM3trait4 <- scantwo(DM3.prob, pheno.col=5, method="imp")
# operm2.imp.DM3trait4 <- scantwo(DM3.prob, pheno.col=5, method="imp", n.perm=100)
# str(operm2.imp.DM3trait4)
# summary(operm2.imp.DM3trait4, alpha=0.01) #1% 7.62 6.28 5.43 6.6 3.35 4.27
# summary(operm2.imp.DM3trait4, alpha=0.05) #5% 6.97 5.57 4.92 5.67 3.22 3.47
# summary(out2.imp.DM3trait4, thresholds=c(6.97,5.57,4.92,5.67,3.22), what="int")
# 
# #stepwise selection of multiple-QTL model
# ##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
# calc.penalties(operm2.imp.DM3trait4, alpha=0.05) #3.468869 4.919449 2.102937
# 
# ##1.4 stepwise selection for trait4
# chr.trait4 <- c("S3D")
# pos.trait4 <- c(70.20)
# qtl.trait4 <- makeqtl(DM3.sim, chr.trait4, pos.trait4)
# qtl.trait4
# my.formula=y~Q1
# stepout.DM3trait4=stepwiseqtl(DM3.sim, pheno.col=5, additive.only=FALSE, penalties=c(3.47, 4.92, 2.10), qtl = qtl.trait4, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
# stepout.DM3trait4
# 
# 
# 
# chr.trait4 <- c("S3D","S7A","S7B")
# pos.trait4 <- c(70.204,54.869,53.841)
# qtl.trait4 <- makeqtl(DM3.sim, chr.trait4, pos.trait4)
# my.formula=y~Q1+Q2+Q3
# out.fq.trait4 <- fitqtl(DM3.sim, pheno.col=5, qtl=qtl.trait4, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
# out.fq.trait4
# 
# 
# names(attributes(stepout.DM3trait4))
# trait4.trace=attr(stepout.DM3trait4, "trace")
# par(mfrow=c(3,4))
# for(i in seq(along=trait4.trace))
#   plotModel(trait4.trace[[i]], chronly=TRUE,
#             main=paste(i, ": pLOD =",
#                        round(attr(trait4.trace[[i]], "pLOD"), 2)))
# 
# rqtl.trait4 =refineqtl(DM3.sim, pheno.col=5, qtl=qtl.trait4, method="imp", formula=my.formula)
# rqtl.trait4
# 
# tiff(filename="Fig. QTL for trait4.tiff", width=8, height=6,units='in',res=300)
# plot(out.imp.trait4, chr=c("1A", "1B", "5A", "7A"), col=4, ylim=c(0,85), ylab="LOD", lwd=1)
# plotLodProfile(rqtl.trait4, col=2, add=T, qtl.labels = FALSE, lwd=1)
# abline(h=4.58, col=3, lty=2)
# dev.off()
# 
# mar.S7B.trait4=find.marker(DM3.sim, chr="S7B",pos=32.547)
# mar.S7B.trait4
# 
# bayesint(rqtl.trait4, qtl.index=1, prob=0.99)
# 
# lodint(rqtl.trait4,qtl.index=1, drop=1.5)
# 
# # tiff(filename="Fig. QTL interaction for trait4.tiff", width=6, height=6,units='in',res=300)
# # effectplot(DM3.sim, pheno.col=3, mname1=mar.1A.trait4, geno1=c("Rusty", "PI466979"), ylim=c(2,7), mname2=mar.7A.trait4, geno2=c("Rusty", "PI466979"), xlab="7A@76.0", ylab="Resistance Score", legend.lab="1A@144.0")
# # dev.off()
# 
