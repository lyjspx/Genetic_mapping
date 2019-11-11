#Set working directory under Mac
setwd("C:/Users/liu.yuan/Google Drive/Alfalfa/Alfalfa_tetra_map")

#load packages "qtl"
DM5=read.cross("csvr", "E:/Google Drive/Alfalfa/Alfalfa_tetra_map", "HB_map_pheno_r.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "bc")

#trait1 = matDupy12OCT
#trait2 = matPepsi13July
#trait3 = matPepsi14Sep
#trait4 = matPepsi14Nov

str(DM5)
nind(DM5)
nchr(DM5)
totmar(DM5)
nmar(DM5)
plot(DM5)
plotMap(DM5)
plotMissing(DM5)
summaryMap(DM5)

write.table(summary.map(DM5),file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")

##########################################################
#Spearman's correlations among the phenotypes
DM5.pheno = DM5$pheno
str(DM5.pheno)
mean(DM5.pheno[,2])
min(DM5.pheno[,2])
max(DM5.pheno[,2])

mean(DM5.pheno[,3])
min(DM5.pheno[,3])
max(DM5.pheno[,3])

mean(DM5.pheno[,4])
min(DM5.pheno[,4])
max(DM5.pheno[,4])

mean(DM5.pheno[,5])
min(DM5.pheno[,5])
max(DM5.pheno[,5])

cor(DM5.pheno[,2], DM5.pheno[,3])
cor.test(DM5.pheno[,2], DM5.pheno[,3], method="spearman")
cor.test(DM5.pheno[,2], DM5.pheno[,3], method="pearson")

cor(DM5.pheno[,2], DM5.pheno[,4])
cor.test(DM5.pheno[,2], DM5.pheno[,4], method="spearman")
cor.test(DM5.pheno[,2], DM5.pheno[,4], method="pearson")

cor(DM5.pheno[,3], DM5.pheno[,4])
cor.test(DM5.pheno[,3], DM5.pheno[,4], method="spearman")
cor.test(DM5.pheno[,3], DM5.pheno[,4], method="pearson")
boxplot(DM5.pheno[,2:4], ylab="Score")

##########################################################
#remove samples with a lot of missing marker values

nind(DM5)
par(mfrow=c(1,1), las=1)
plot(ntyped(DM5), xlab="Sample", ylab="No. typed markers", main="No. typed markers by samples") 
DM5 <- subset(DM5, ind=(ntyped(DM5)>(totmar(DM5)/2)))
nind(DM5)

#identify duplicate individuals
# cg <- comparegeno(DM5)
# str(cg)
# hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes") 
# rug(cg[lower.tri(cg)])
# wh <- which(cg > 0.9, arr=TRUE)
# wh <- wh[wh[,1] < wh[,2],]
# wh
# genotype=pull.geno(DM5)
# str(genotype)
# dim(genotype)
# table(genotype[70,], genotype[81,])
# 
# gc=genClones(DM5, tol=0.8)
# gc$cgd
# DM5=fixClones(DM5, gc$cgd, consensus = TRUE)
# nind(DM5)
#################
#plot genetic map of DM5
tiff(filename="Fig. Genetic map of DM5.tiff", width=8, height=6,units='in',res=300)
plotMap(DM5, main="Genetic map of DM5" )
dev.off()

###################################################################
#single QTL genome scan
#calculate conditional probability if using "imp" method
DM5.prob <- calc.genoprob(DM5, step=1, error.prob=0.01)
str(DM5.prob)

#simulate missing genotype if using "imp" method
DM5=jittermap(DM5, amount=1e-6)
DM5.sim <- sim.geno(DM5, step=2, n.draws=128, err=0.001)
str(DM5.sim)

##1.1 single-QTL genome scan for trait1
##scanone wiht method "imp"
out.imp.trait1 <- scanone(DM5.sim, method="imp", pheno.col=2)
summary(out.imp.trait1)
plot(out.imp.trait1)

operm.imp.trait1 <- scanone(DM5.sim, pheno.col=2, method="imp", n.perm=1000)
summary(operm.imp.trait1, alpha=0.05) #3.3
summary(operm.imp.trait1, alpha=0.01) #4

###No QTL


###########################################
####Trait2#################
##scanone wiht method "imp"
out.imp.trait2 <- scanone(DM5.sim, method="imp", pheno.col=3)
summary(out.imp.trait2,threshold=3.18)
plot(out.imp.trait2)

operm.imp.trait2 <- scanone(DM5.sim, pheno.col=3, method="imp", n.perm=1000)
summary(operm.imp.trait2, alpha=0.05) #3.18
summary(operm.imp.trait2, alpha=0.01) #


plot(out.imp.trait2,main = "QTL mapping of matPepsi13July in DM5")
abline(h=3.18, col=3, lty=2)
text(1750,3.5,"1000 permutation at p=0.05, threshold = 3.18",col = 3)


lodint(out.imp.trait2,chr="1C",drop=1.5)
lodint(out.imp.trait2,chr="1D",drop=1.5)
lodint(out.imp.trait2,chr="2C",drop=1.5)
lodint(out.imp.trait2,chr="2D",drop=1.5)

drop1.5 = rbind(lodint(out.imp.trait2,chr="1C",drop=1.5),
                lodint(out.imp.trait2,chr="1D",drop=1.5),
                lodint(out.imp.trait2,chr="2C",drop=1.5),
                lodint(out.imp.trait2,chr="2D",drop=1.5)
                )

write.table(drop1.5,file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")


chr.trait2 <- c("2D")
pos.trait2 <- c(36.055)
qtl.trait2 <- makeqtl(DM5.sim, chr.trait2, pos.trait2)
my.formula=y~Q1
out.fq.trait2 <- fitqtl(DM5.sim, pheno.col=3, qtl=qtl.trait2, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.trait2



#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for trait2
out2.imp.DM5trait2 <- scantwo(DM5.prob, pheno.col=3, method="imp")
operm2.imp.DM5trait2 <- scantwo(DM5.prob, pheno.col=3, method="imp", n.perm=100)
str(operm2.imp.DM5trait2)
summary(operm2.imp.DM5trait2, alpha=0.01) #1% 11.4 9.88 8.56 8.16 7.47 4.49
summary(operm2.imp.DM5trait2, alpha=0.05) #5%  9.75 8.07 5.88 7.15 6.34 3.69
summary(out2.imp.DM5trait2, thresholds=c(9.75, 8.07, 5.88, 7.15, 6.34), what="int") #No pair

#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.DM5trait2, alpha=0.05)

##1.4 stepwise selection for trait2
 
chr.trait2 <- c("1")
pos.trait2<- c(5.104)
qtl.trait2 <- makeqtl(DM5.sim, chr.trait2, pos.trait2)
qtl.trait2
my.formula=y~Q1
stepout.DM5trait2=stepwiseqtl(DM5.sim, pheno.col=3, additive.only=FALSE, penalties=c(3.69, 5.88, 4.37), qtl = qtl.trait2, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.DM5trait2


#only QTL in chromosome 1 left

chr.trait2 <- c("1")
pos.trait2 <- c(6)
qtl.trait2 <- makeqtl(DM5.sim, chr.trait2, pos.trait2)
my.formula=y~Q1
out.fq.trait2 <- fitqtl(DM5.sim, pheno.col=3, qtl=qtl.trait2, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.trait2


names(attributes(stepout.DM5trait2))
trait2.trace=attr(stepout.DM5trait2, "trace")
par(mfrow=c(3,4))
for(i in seq(along=trait2.trace))
  plotModel(trait2.trace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =",
                       round(attr(trait2.trace[[i]], "pLOD"), 2)))

rqtl.trait2 =refineqtl(DM5.sim, pheno.col=3, qtl=qtl.trait2, method="imp", formula=my.formula)
rqtl.trait2

tiff(filename="Fig. QTL for trait2.tiff", width=8, height=6,units='in',res=300)
plot(out.imp.trait2, chr=c("1A", "1B", "5A", "7A"), col=4, ylim=c(0,85), ylab="LOD", lwd=1)
plotLodProfile(rqtl.trait2, col=2, add=T, qtl.labels = FALSE, lwd=1)
abline(h=4.58, col=3, lty=2)
dev.off()

mar.S7B.trait2=find.marker(DM5.sim, chr="S7B",pos=32.547)
mar.S7B.trait2
bayesint(rqtl.trait2, qtl.index=1, prob=0.99)

lodint(rqtl.trait2,qtl.index=1, drop=1.5)


# tiff(filename="Fig. QTL interaction for trait2.tiff", width=6, height=6,units='in',res=300)
# effectplot(DM5.sim, pheno.col=3, mname1=mar.1A.trait2, geno1=c("Rusty", "PI466979"), ylim=c(2,7), mname2=mar.7A.trait2, geno2=c("Rusty", "PI466979"), xlab="7A@76.0", ylab="Resistance Score", legend.lab="1A@144.0")
# dev.off()

###########################################
####Trait3#################
##scanone wiht method "imp"
out.imp.trait3 <- scanone(DM5.sim, method="imp", pheno.col=4)
summary(out.imp.trait3)
plot(out.imp.trait3)

operm.imp.trait3 <- scanone(DM5.sim, pheno.col=4, method="imp", n.perm=1000)
summary(operm.imp.trait3, alpha=0.05) #3.22
summary(operm.imp.trait3, alpha=0.01) #3.88
#No QTL

#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for trait3
out2.imp.DM5trait3 <- scantwo(DM5.prob, pheno.col=3, method="imp")
operm2.imp.DM5trait3 <- scantwo(DM5.prob, pheno.col=3, method="imp", n.perm=100)
str(operm2.imp.DM5trait3)
summary(operm2.imp.DM5trait3, alpha=0.01) #1% 12.2 11.2 9.15 8.34 6.98 4.11
summary(operm2.imp.DM5trait3, alpha=0.05) #5% 11.2 10.4 8.3 7.43 6.32 3.88
summary(out2.imp.DM5trait3, thresholds=c(11.2, 10.4, 8.3, 7.43, 6.32), what="int")

#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.DM5trait3, alpha=0.05)

##1.4 stepwise selection for trait3


chr.trait3 <- c("S7B")
pos.trait3 <- c(32.547)
qtl.trait3 <- makeqtl(DM5.sim, chr.trait3, pos.trait3)
my.formula=y~Q1
out.fq.trait3 <- fitqtl(DM5.sim, pheno.col=3, qtl=qtl.trait3, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.trait3


names(attributes(stepout.DM5trait3))
trait3.trace=attr(stepout.DM5trait3, "trace")
par(mfrow=c(3,4))
for(i in seq(along=trait3.trace))
  plotModel(trait3.trace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =",
                       round(attr(trait3.trace[[i]], "pLOD"), 2)))

rqtl.trait3 =refineqtl(DM5.sim, pheno.col=3, qtl=qtl.trait3, method="imp", formula=my.formula)
rqtl.trait3

tiff(filename="Fig. QTL for trait3.tiff", width=8, height=6,units='in',res=300)
plot(out.imp.trait3, chr=c("1A", "1B", "5A", "7A"), col=4, ylim=c(0,85), ylab="LOD", lwd=1)
plotLodProfile(rqtl.trait3, col=2, add=T, qtl.labels = FALSE, lwd=1)
abline(h=4.58, col=3, lty=2)
dev.off()

mar.S7B.trait3=find.marker(DM5.sim, chr="S7B",pos=32.547)
mar.S7B.trait3

bayesint(rqtl.trait3, qtl.index=1, prob=0.99)

lodint(rqtl.trait3,qtl.index=1, drop=1.5)

# tiff(filename="Fig. QTL interaction for trait3.tiff", width=6, height=6,units='in',res=300)
# effectplot(DM5.sim, pheno.col=3, mname1=mar.1A.trait3, geno1=c("Rusty", "PI466979"), ylim=c(2,7), mname2=mar.7A.trait3, geno2=c("Rusty", "PI466979"), xlab="7A@76.0", ylab="Resistance Score", legend.lab="1A@144.0")
# dev.off()

###########################################
####trait4#################
##scanone wiht method "imp"
out.imp.trait4 <- scanone(DM5.sim, method="imp", pheno.col=5)
summary(out.imp.trait4)
plot(out.imp.trait4)

operm.imp.trait4 <- scanone(DM5.sim, pheno.col=5, method="imp", n.perm=1000)
summary(operm.imp.trait4, alpha=0.05) #3.28
summary(operm.imp.trait4, alpha=0.01) #3.84
#No QTL

summary(out.imp.trait4,threshold=c(3.33))
#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for trait4
out2.imp.DM5trait4 <- scantwo(DM5.prob, pheno.col=5, method="imp")
operm2.imp.DM5trait4 <- scantwo(DM5.prob, pheno.col=5, method="imp", n.perm=100)
str(operm2.imp.DM5trait4)
summary(operm2.imp.DM5trait4, alpha=0.01) #1% 7.62 6.28 5.43 6.6 3.35 4.27
summary(operm2.imp.DM5trait4, alpha=0.05) #5% 6.97 5.57 4.92 5.67 3.22 3.47
summary(out2.imp.DM5trait4, thresholds=c(6.97,5.57,4.92,5.67,3.22), what="int")

#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.DM5trait4, alpha=0.05) #3.468869 4.919449 2.102937

##1.4 stepwise selection for trait4
chr.trait4 <- c("S3D")
pos.trait4 <- c(70.20)
qtl.trait4 <- makeqtl(DM5.sim, chr.trait4, pos.trait4)
qtl.trait4
my.formula=y~Q1
stepout.DM5trait4=stepwiseqtl(DM5.sim, pheno.col=5, additive.only=FALSE, penalties=c(3.47, 4.92, 2.10), qtl = qtl.trait4, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.DM5trait4



chr.trait4 <- c("S3D","S7A","S7B")
pos.trait4 <- c(70.204,54.869,53.841)
qtl.trait4 <- makeqtl(DM5.sim, chr.trait4, pos.trait4)
my.formula=y~Q1+Q2+Q3
out.fq.trait4 <- fitqtl(DM5.sim, pheno.col=5, qtl=qtl.trait4, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.trait4


names(attributes(stepout.DM5trait4))
trait4.trace=attr(stepout.DM5trait4, "trace")
par(mfrow=c(3,4))
for(i in seq(along=trait4.trace))
  plotModel(trait4.trace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =",
                       round(attr(trait4.trace[[i]], "pLOD"), 2)))

rqtl.trait4 =refineqtl(DM5.sim, pheno.col=5, qtl=qtl.trait4, method="imp", formula=my.formula)
rqtl.trait4

tiff(filename="Fig. QTL for trait4.tiff", width=8, height=6,units='in',res=300)
plot(out.imp.trait4, chr=c("1A", "1B", "5A", "7A"), col=4, ylim=c(0,85), ylab="LOD", lwd=1)
plotLodProfile(rqtl.trait4, col=2, add=T, qtl.labels = FALSE, lwd=1)
abline(h=4.58, col=3, lty=2)
dev.off()

mar.S7B.trait4=find.marker(DM5.sim, chr="S7B",pos=32.547)
mar.S7B.trait4

bayesint(rqtl.trait4, qtl.index=1, prob=0.99)

lodint(rqtl.trait4,qtl.index=1, drop=1.5)

# tiff(filename="Fig. QTL interaction for trait4.tiff", width=6, height=6,units='in',res=300)
# effectplot(DM5.sim, pheno.col=3, mname1=mar.1A.trait4, geno1=c("Rusty", "PI466979"), ylim=c(2,7), mname2=mar.7A.trait4, geno2=c("Rusty", "PI466979"), xlab="7A@76.0", ylab="Resistance Score", legend.lab="1A@144.0")
# dev.off()

