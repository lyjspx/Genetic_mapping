#Set working directory under Mac
setwd("/Users/desert/Google Drive/CodeCenter/rdata")

#load packages "qtl"
RP979=read.cross("csvr", "/Users/desert/Google Drive/CodeCenter/rdata", "RP979_rotated_Rqtl_6624mapped_JoinMap_tan_spot.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "riself")

str(RP979)
nind(RP979)
nchr(RP979)
totmar(RP979)
nmar(RP979)
plot(RP979)
plotMap(RP979)
plotMissing(RP979)
summaryMap(RP979)

##########################################################
#Spearman's correlations among the phenotypes
str(RP979$pheno)
RP979.pheno=RP979$pheno[complete.cases(RP979$pheno[,2:4]),]
str(RP979.pheno)
mean(RP979.pheno[,2])
min(RP979.pheno[,2])
max(RP979.pheno[,2])

mean(RP979.pheno[,3])
min(RP979.pheno[,3])
max(RP979.pheno[,3])

mean(RP979.pheno[,4])
min(RP979.pheno[,4])
max(RP979.pheno[,4])

cor(RP979.pheno[,2], RP979.pheno[,3])
cor.test(RP979.pheno[,2], RP979.pheno[,3], method="spearman")
cor.test(RP979.pheno[,2], RP979.pheno[,3], method="pearson")

cor(RP979.pheno[,2], RP979.pheno[,4])
cor.test(RP979.pheno[,2], RP979.pheno[,4], method="spearman")
cor.test(RP979.pheno[,2], RP979.pheno[,4], method="pearson")

cor(RP979.pheno[,3], RP979.pheno[,4])
cor.test(RP979.pheno[,3], RP979.pheno[,4], method="spearman")
cor.test(RP979.pheno[,3], RP979.pheno[,4], method="pearson")
boxplot(RP979.pheno[,2:4], ylab="Score")

##########################################################
#remove samples with a lot of missing marker values
par(mfrow=c(1,1), las=1)
plot(ntyped(RP979), xlab="Sample", ylab="No. typed markers", main="No. typed markers by samples") 
RP979 <- subset(RP979, ind=(ntyped(RP979)>2800))
nind(RP979)

#identify duplicate individuals
cg <- comparegeno(RP979)
str(cg)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes") 
rug(cg[lower.tri(cg)])
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh
genotype=pull.geno(RP979)
str(genotype)
dim(genotype)
table(genotype[70,], genotype[81,])

gc=genClones(RP979, tol=0.8)
gc$cgd
RP979=fixClones(RP979, gc$cgd, consensus = TRUE)

#################
#plot genetic map of RP979
tiff(filename="Fig. Genetic map of RP979.tiff", width=8, height=6,units='in',res=300)
plotMap(RP979, main="Genetic map of RP979" )
dev.off()

###################################################################
#single QTL genome scan
#calculate conditional probability if using "imp" method
RP979.prob <- calc.genoprob(RP979, step=1, error.prob=0.01)
str(RP979.prob)

#simulate missing genotype if using "imp" method
RP979=jittermap(RP979, amount=1e-6)
RP979.sim <- sim.geno(RP979, step=2, n.draws=128, err=0.001)
str(RP979.sim)

##1.1 single-QTL genome scan for TTKSK
##scanone wiht method "imp"
out.imp.TTKSK <- scanone(RP979.sim, method="imp", pheno.col=2)
summary(out.imp.TTKSK)
plot(out.imp.TTKSK)
#lodint(out.imp.TTKSK,chr="1B",drop=1.5)

operm.imp.TTKSK <- scanone(RP979.sim, pheno.col=4, method="imp", n.perm=1000)
summary(operm.imp.TTKSK, alpha=0.05)
summary(operm.imp.TTKSK, alpha=0.01)

#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for TTKSK
out2.imp.RP979TTKSK <- scantwo(RP979.prob, pheno.col=2, method="imp")
operm2.imp.RP979TTKSK <- scantwo(RP979.prob, pheno.col=2, method="imp", n.perm=100)
str(operm2.imp.RP979TTKSK)
summary(operm2.imp.RP979TTKSK, alpha=0.01)
summary(operm2.imp.RP979TTKSK, alpha=0.05)
summary(out2.imp.RP979TTKSK, thresholds=c(7.29, 5.7, 4.97, 6.08, 3.04), what="int")

##multiple QTL mapping and forward/backward selection using "imp" method
###1.3 forward/backward selection for TTKSK, this can be easily done using "stepwiseqtl"
chr.TTKSK <- c("1B")
pos.TTKSK <- c(52.876)
qtl.TTKSK <- makeqtl(RP979.sim, chr.TTKSK, pos.TTKSK)
qtl.TTKSK
my.formula=y~Q1
out.fq.TTKSK <- fitqtl(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, formula=my.formula, get.ests=TRUE)
out.fq.TTKSK

out.aq.TTKSK=addqtl(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, formula=my.formula)
max(out.aq.TTKSK)
plot(out.aq.TTKSK)
summary(out.aq.TTKSK, threshold=c(3))

qtl.TTKSK=addtoqtl(RP979.sim, qtl.TTKSK, c("7A","5A"), c(74.5, 129.1))
qtl.TTKSK

my.formula=y~Q1 + Q2 + Q3
out.aq.TTKSK=addqtl(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, formula=my.formula)
plot(out.aq.TTKSK)
summary(out.aq.TTKSK, threshold=c(3))

qtl.TTKSK=addtoqtl(RP979.sim, qtl.TTKSK, c("3B"), c(284))
qtl.TTKSK

my.formula=y~Q1 + Q2 + Q3 + Q4
out.fq.TTKSK <- fitqtl(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, formula=my.formula, get.ests=TRUE)
out.fq.TTKSK

addint(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, formula=my.formula)

out.fq.TTKSK <- fitqtl(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, formula=my.formula, get.ests=TRUE)
out.fq.TTKSK

rqtl.TTKSK <- refineqtl(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, formula=my.formula)
rqtl.TTKSK
plotLodProfile(rqtl.TTKSK, ylim=c(0,90))
abline(h=3.9, col=3, lty=2)

out.fq.TTKSK <- fitqtl(RP979.sim, pheno.col=2, qtl=rqtl.TTKSK, formula=my.formula, get.ests=TRUE)
out.fq.TTKSK

#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.RP979TTKSK, alpha=0.05)

##1.4 stepwise selection for TTKSK
chr.TTKSK <- c("1B")
pos.TTKSK <- c(52.876)
qtl.TTKSK <- makeqtl(RP979.sim, chr.TTKSK, pos.TTKSK)
qtl.TTKSK
my.formula=y~Q1
stepout.RP979TTKSK=stepwiseqtl(RP979.sim, pheno.col=2, additive.only=FALSE, penalties=c(3.55, 4.97, 2.15), qtl = qtl.TTKSK, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP979TTKSK
stepout.RP979TTKSK=stepwiseqtl(RP979.prob, pheno.col=2, additive.only=FALSE, penalties=c(3.55, 4.97, 2.15), qtl = qtl.TTKSK, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP979TTKSK

names(attributes(stepout.RP979TTKSK))
TTKSK.trace=attr(stepout.RP979TTKSK, "trace")
par(mfrow=c(4,4))
for(i in seq(along=TTKSK.trace))
  plotModel(TTKSK.trace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =",
                       round(attr(TTKSK.trace[[i]], "pLOD"), 2)))

chr.TTKSK <- c("1A", "1B", "5A", "7A")
pos.TTKSK <- c(144, 52.876, 176, 76)
qtl.TTKSK <- makeqtl(RP979.sim, chr.TTKSK, pos.TTKSK)
my.formula=y~Q1 + Q2 + Q3 + Q4 +Q1:Q4
out.fq.TTKSK <- fitqtl(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.TTKSK

rqtl.TTKSK =refineqtl(RP979.sim, pheno.col=2, qtl=qtl.TTKSK, method="imp", formula=my.formula)
rqtl.TTKSK

tiff(filename="Fig. QTL for TTKSK.tiff", width=8, height=6,units='in',res=300)
plot(out.imp.TTKSK, chr=c("1A", "1B", "5A", "7A"), col=4, ylim=c(0,85), ylab="LOD", lwd=1)
plotLodProfile(rqtl.TTKSK, col=2, add=T, qtl.labels = FALSE, lwd=1)
abline(h=4.58, col=3, lty=2)
dev.off()

mar.1A.TTKSK=find.marker(RP979.sim, chr="1A",pos=144)
mar.1A.TTKSK
mar.1B.TTKSK=find.marker(RP979.sim, chr="1B",pos=52.876)
mar.1B.TTKSK
mar.5A.TTKSK=find.marker(RP979.sim, chr="5A",pos=176)
mar.5A.TTKSK
mar.7A.TTKSK=find.marker(RP979.sim, chr="7A",pos=76)
mar.7A.TTKSK

bayesint(rqtl.TTKSK, qtl.index=1, prob=0.99)
bayesint(rqtl.TTKSK, qtl.index=2, prob=0.99)
bayesint(rqtl.TTKSK, qtl.index=3, prob=0.99)
bayesint(rqtl.TTKSK, qtl.index=4, prob=0.99)

lodint(rqtl.TTKSK,qtl.index=1, drop=1.5)
lodint(rqtl.TTKSK,qtl.index=2, drop=1.5)
lodint(rqtl.TTKSK,qtl.index=3, drop=1.5)
lodint(rqtl.TTKSK,qtl.index=4, drop=1.5)

out.fq.TTKSK <- fitqtl(RP979.sim, pheno.col=2, qtl=rqtl.TTKSK, method="imp", formula=my.formula, get.ests=TRUE)
out.fq.TTKSK

tiff(filename="Fig. QTL interaction for TTKSK.tiff", width=6, height=6,units='in',res=300)
effectplot(RP979.sim, pheno.col=2, mname1=mar.1A.TTKSK, geno1=c("Rusty", "PI466979"), ylim=c(2,7), mname2=mar.7A.TTKSK, geno2=c("Rusty", "PI466979"), xlab="7A@76.0", ylab="Resistance Score", legend.lab="1A@144.0")
dev.off()

##################################################
###########################################
#2.1 single-QTL genome scan for TRTTF
out.imp.TRTTF <- scanone(RP979.sim, method="imp", pheno.col=3)
summary(out.imp.TRTTF)
plot(out.imp.TRTTF)
lodint(out.imp.TRTTF,chr="1B",drop=1.5)
lodint(out.imp.TRTTF,chr="2B",drop=1.5)

operm.imp.TRTTF <- scanone(RP979.sim, pheno.col=3, method="imp", n.perm=1000)
summary(operm.imp.TRTTF, alpha=0.05)
summary(operm.imp.TRTTF, alpha=0.01)

#2.2 two-QTL genome scan with "scantwo" for TRTTF
out2.imp.RP979TRTTF <- scantwo(RP979.prob, pheno.col=3, method="imp")
operm2.imp.RP979TRTTF <- scantwo(RP979.prob, pheno.col=3, method="imp", n.perm=100)
str(operm2.imp.RP979TRTTF)
summary(operm2.imp.RP979TRTTF, alpha=0.01)
summary(operm2.imp.RP979TRTTF, alpha=0.05)

summary(out2.imp.RP979TRTTF, thresholds=c(7.18, 5.48, 4.73, 5.68, 2.86), what="int")
plot(out2.imp.RP979TRTTF)

###2.3 forward/backward selection for TRTTF
chr.TRTTF <- c("1B", "2B")
pos.TRTTF <- c(52.876, 116.7)
qtl.TRTTF <- makeqtl(RP979.sim, chr.TRTTF, pos.TRTTF)
qtl.TRTTF
my.formula=y~Q1 + Q2
out.fq.TRTTF <- fitqtl(RP979.sim, pheno.col=3, qtl=qtl.TRTTF, formula=my.formula, get.ests=TRUE)
out.fq.TRTTF

out.aq.TRTTF=addqtl(RP979.sim, pheno.col=3, qtl=qtl.TRTTF, formula=my.formula)
max(out.aq.TRTTF)
plot(out.aq.TRTTF)
summary(out.aq.TRTTF, threshold=c(3))

addint(RP979.sim, qtl=qtl.TRTTF, formula=my.formula)

out.ap.TRTTF=addpair(RP979.sim, pheno.col=3, chr="1A", qtl=qtl.TRTTF, formula=my.formula)
max(out.ap.TRTTF)
plot(out.ap.TRTTF)

rqtl.TRTTF <- refineqtl(RP979.sim, pheno.col=3, qtl=qtl.TRTTF, formula=my.formula)
rqtl.TRTTF
plot(out.imp.TRTTF, chr=c("1B", "2B"), ylim=c(0,20))
plotLodProfile(rqtl.TRTTF, ylim=c(0,90), col=2,add=T)
abline(h=3.9, col=3, lty=2)

out.fq.TRTTF <- fitqtl(RP979.sim, pheno.col=3, qtl=rqtl.TRTTF, formula=my.formula, get.ests=TRUE)
out.fq.TRTTF

##2.4 stepwise selection for TRTTF
calc.penalties(operm2.imp.RP979TRTTF, alpha=0.05)

chr.TRTTF <- c("1B", "2B")
pos.TRTTF <- c(52.876, 116.7)
qtl.TRTTF <- makeqtl(RP979.sim, chr.TRTTF, pos.TRTTF)
qtl.TRTTF
my.formula=y~Q1 + Q2
stepout.RP979TRTTF=stepwiseqtl(RP979.sim, pheno.col=3, additive.only=FALSE, penalties=c(3.31, 4.73, 2.17), qtl = qtl.TRTTF, formula=y~Q1 + Q2,  max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP979TRTTF
stepout.RP979TRTTF=stepwiseqtl(RP979.prob, pheno.col=3, additive.only=FALSE, penalties=c(3.31, 4.73, 2.17), qtl = qtl.TRTTF, formula=y~Q1 + Q2, max.qtl=6, refine.locations=TRUE, scan.pairs=TRUE, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP979TRTTF
names(attributes(stepout.RP979TRTTF))
TRTTF.trace=attr(stepout.RP979TRTTF, "trace")

par(mfrow=c(1,1))
plotLodProfile(stepout.RP979TRTTF)
par(mfrow=c(4,4))
for(i in seq(along=TRTTF.trace))
  plotModel(TRTTF.trace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =",
                       round(attr(TRTTF.trace[[i]], "pLOD"), 2)))

chr.TRTTF <- c("1B", "2B")
pos.TRTTF <- c(52.876, 116.701)
qtl.TRTTF <- makeqtl(RP979.sim, chr.TRTTF, pos.TRTTF)
qtl.TRTTF
my.formula=y~Q1 + Q2
out.fq.TRTTF <- fitqtl(RP979.sim, pheno.col=3, qtl=rqtl.TRTTF, formula=y~Q1 + Q2, get.ests=TRUE)
rqtl.TRTTF =refineqtl(RP979.sim, pheno.col=3, qtl=qtl.TRTTF, method="imp", formula=y~Q1 + Q2)
rqtl.TRTTF

mar.1B.TRTTF=find.marker(RP979.sim, chr="1B",pos=52.876)
mar.1B.TRTTF
mar.2B.TRTTF=find.marker(RP979.sim, chr="2B",pos=116.701)
mar.2B.TRTTF

tiff(filename="Fig. QTL for TRTTF.tiff", width=8, height=6,units='in',res=300)
plot(out.imp.TRTTF, chr=c("1B", "2B"), col=4, ylim=c(0,20), ylab="LOD", lwd=1)
plotLodProfile(rqtl.TRTTF, add=T, col=2, qtl.labels = FALSE, lwd=1)
abline(h=4.08, col=3, lty=2)
dev.off()

lodint(rqtl.TRTTF,qtl.index=1, drop=1.5)
lodint(rqtl.TRTTF,qtl.index=2, drop=1.5)

out.fq.TRTTF <- fitqtl(RP979.sim, pheno.col=3, qtl=rqtl.TRTTF, formula=y~Q1 + Q2, get.ests=TRUE)
out.fq.TRTTF <- fitqtl(RP979.sim, pheno.col=3, qtl=rqtl.TRTTF, drop=FALSE, formula=y~Q1 + Q2, get.ests=TRUE)
out.fq.TRTTF

#############################################################
##3.1 single-QTL genome scan for TMLKC
out.imp.TMLKC <- scanone(RP979.sim, method="imp", pheno.col=4)
summary(out.imp.TMLKC)
plot(out.imp.TMLKC)

operm.imp.TMLKC <- scanone(RP979.sim, pheno.col=4, method="imp", n.perm=1000)
summary(operm.imp.TMLKC, alpha=0.05)
summary(operm.imp.TMLKC, alpha=0.01)

#3.2 two-QTL genome scan with "scantwo" for TMLKC
out2.imp.RP979TMLKC <- scantwo(RP979.prob, pheno.col=4, method="imp")
operm2.imp.RP979TMLKC <- scantwo(RP979.prob, pheno.col=4, method="imp", n.perm=100)
summary(operm2.imp.RP979TMLKC, alpha=0.01)
summary(operm2.imp.RP979TMLKC, alpha=0.05)

summary(out2.imp.RP979TMLKC, thresholds=c(7.05, 5.57, 4.78, 6.16, 3.17), what="int")
plot(out2.imp.RP979TMLKC)

###3.3 forward/backward selection for TMLKC
chr.TMLKC <- c("1B")
pos.TMLKC <- c(52.876)
qtl.TMLKC <- makeqtl(RP979.sim, chr.TMLKC, pos.TMLKC)
qtl.TMLKC
my.formula=y~Q1
out.fq.TMLKC <- fitqtl(RP979.sim, pheno.col=4, qtl=qtl.TMLKC, formula=my.formula, get.ests=TRUE)
out.fq.TMLKC

out.aq.TMLKC=addqtl(RP979.sim, pheno.col=4, qtl=qtl.TMLKC, formula=my.formula)
max(out.aq.TMLKC)
plot(out.aq.TMLKC)
summary(out.aq.TMLKC, threshold=c(3))

qtl.TMLKC=addtoqtl(RP979.sim, qtl.TMLKC, c("7B","5A"), c(180.2, 140.4))
qtl.TMLKC

my.formula=y~Q1 + Q2 + Q3
out.aq.TMLKC=addqtl(RP979.sim, pheno.col=4, qtl=qtl.TMLKC, formula=my.formula)
plot(out.aq.TMLKC)
summary(out.aq.TMLKC, threshold=c(3))

qtl.TMLKC=addtoqtl(RP979.sim, qtl.TMLKC, c("3B"), c(284))
qtl.TMLKC

my.formula=y~Q1 + Q2 + Q3
out.fq.TMLKC <- fitqtl(RP979.sim, pheno.col=4, qtl=qtl.TMLKC, formula=my.formula, get.ests=TRUE)
out.fq.TMLKC

addint(RP979.sim, qtl=qtl.TMLKC, formula=my.formula)

out.ap.TMLKC=addpair(RP979.sim, pheno.col=4, chr="1A", qtl=qtl.TMLKC, formula=my.formula)
max(out.ap.TMLKC)
plot(out.ap.TMLKC)

rqtl.TMLKC <- refineqtl(RP979.sim, pheno.col=4, qtl=qtl.TMLKC, formula=my.formula)
rqtl.TMLKC

out.fq.TMLKC <- fitqtl(RP979.sim, pheno.col=4, qtl=rqtl.TMLKC, formula=my.formula, get.ests=TRUE)
out.fq.TMLKC

##3.4 stepwise selection for TMLKC
calc.penalties(operm2.imp.RP979TMLKC, alpha=0.05)

chr.TMLKC <- c("1B")
pos.TMLKC <- c(52.876)
qtl.TMLKC <- makeqtl(RP979.sim, chr.TMLKC, pos.TMLKC)
qtl.TMLKC
my.formula=y~Q1

stepout.RP979TMLKC=stepwiseqtl(RP979.sim, pheno.col=4, additive.only=FALSE, penalties=c(3.56, 4.78, 2.02), qtl = qtl.TMLKC, formula=y~Q1,  max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP979TMLKC
stepout.RP979TMLKC=stepwiseqtl(RP979.prob, pheno.col=4, additive.only=FALSE, penalties=c(3.56, 4.78, 2.02), qtl = qtl.TMLKC, formula=y~Q1, max.qtl=6, refine.locations=TRUE, scan.pairs=TRUE, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP979TMLKC
names(attributes(stepout.RP979TMLKC))
TMLKC.trace=attr(stepout.RP979TMLKC, "trace")

par(mfrow=c(1,1))
plotLodProfile(stepout.RP979TMLKC)
abline(h=4.5, col=3, lty=2)
par(mfrow=c(4,4))
for(i in seq(along=TMLKC.trace))
  plotModel(TMLKC.trace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =",
                       round(attr(TMLKC.trace[[i]], "pLOD"), 2)))

chr.TMLKC <- c("1B", "2B", "5A", "7A", "7B")
pos.TMLKC <- c(52.876, 114, 83.645, 54.943, 186)
qtl.TMLKC <- makeqtl(RP979.sim, chr.TMLKC, pos.TMLKC)
qtl.TMLKC
my.formula=y~Q1 + Q2 + Q3 + Q4 + Q5
out.fq.TMLKC <- fitqtl(RP979.sim, pheno.col=4, qtl=qtl.TMLKC, formula=y~Q1 + Q2 + Q3 + Q4 + Q5, get.ests=TRUE)
out.fq.TMLKC
rqtl.TMLKC =refineqtl(RP979.sim, pheno.col=4, qtl=qtl.TMLKC, method="imp", formula=y~Q1 + Q2 + Q3 + Q4 + Q5)
rqtl.TMLKC

mar.1B.TMKLC=find.marker(RP979.sim, chr="1B",pos=52.876)
mar.1B.TMKLC
mar.2B.TMKLC=find.marker(RP979.sim, chr="2B",pos=114)
mar.2B.TMKLC
mar.5A.TMKLC=find.marker(RP979.sim, chr="5A",pos=83.645)
mar.5A.TMKLC
mar.7A.TMKLC=find.marker(RP979.sim, chr="7A",pos=54.943)
mar.7A.TMKLC
mar.7B.TMKLC=find.marker(RP979.sim, chr="7B",pos=186)
mar.7B.TMKLC

tiff(filename="Fig. QTL for TMLKC.tiff", width=8, height=6,units='in',res=300)
plot(out.imp.TMLKC, chr=c("1B", "2B", "5A", "7A", "7B"), col=4, ylim=c(0,55), ylab="LOD", lwd=1)
plotLodProfile(rqtl.TMLKC, col=2, qtl.labels = FALSE, add=T, lwd=1)
abline(h=4.19, col=3, lty=2)
dev.off()

lodint(rqtl.TMLKC,qtl.index=1, drop=1.5)
lodint(rqtl.TMLKC,qtl.index=2, drop=1.5)
lodint(rqtl.TMLKC,qtl.index=3, drop=1.5)
lodint(rqtl.TMLKC,qtl.index=4, drop=1.5)
lodint(rqtl.TMLKC,qtl.index=5, drop=1.5)

out.fq.TMLKC <- fitqtl(RP979.sim, pheno.col=4, qtl=rqtl.TMLKC, formula=y~y~Q1 + Q2 + Q3 + Q4 + Q5, get.ests=TRUE)
out.fq.TMLKC <- fitqtl(RP979.sim, pheno.col=4, qtl=rqtl.TMLKC, drop=FALSE, formula=y~y~Q1 + Q2 + Q3 + Q4 + Q5, get.ests=TRUE)
out.fq.TMLKC


###############################
####################################
##composite interval mapping
###stepwise regression is used to select number of markers as cofactors
out.imp.RP979TTKSK.cim <- cim(RP979.sim, pheno.col=2, n.marcovar=4, window=5, method="imp", map.function=c("kosambi"))
summary(out.imp.RP979TTKSK.cim, threshold=c(3))
operm.imp.RP979TTKSK.cim <- cim(RP979.sim, pheno.col=2, n.marcovar=4, method="imp", map.function=c("kosambi"), n.perm=100)
summary(operm.imp.RP979TTKSK.cim)
plot(out.imp.RP979TTKSK.cim, show.marker.names=F, main="RP979_TTKSK QTLs_cim", ylim=c(0,80))
add.cim.covar(out.imp.RP979TTKSK.cim)
abline(h=3,col=3)

out.imp.RP979TRTTF.cim <- cim(RP979.sim, pheno.col=3, n.marcovar=2, window=2, method="imp", map.function=c("kosambi"))
summary(out.imp.RP979TRTTF.cim)
operm.imp.RP979TRTTF.cim <- cim(RP979.sim, pheno.col=3, n.marcovar=2, method="imp", map.function=c("kosambi"), n.perm=100)
summary(operm.imp.RP979TRTTF.cim)
plot(out.imp.RP979TRTTF.cim, show.marker.names=F, main="RP979_TRTTF QTLs_cim")

out.imp.RP979TMLKC.cim <- cim(RP979.sim, pheno.col=4, n.marcovar=5, window=5, method=c("imp"), map.function=c("kosambi"))
summary(out.imp.RP979TMLKC.cim)
operm.imp.RP979TMLKC.cim <- cim(RP979.sim, pheno.col=4, n.marcovar=5, method=c("imp"), map.function=c("kosambi"), n.perm=100)
summary(operm.imp.RP979TMLKC.cim)
plot(out.imp.RP979TMLKC.cim, show.marker.names=F, main="RP979_TMLKC QTLs_cim")
add.cim.covar(out.imp.RP979TMLKC.cim)

tiff(filename="Fig. QTL of stem rust resistance in RP979.tiff", width=12, height=8,units='in',res=300)
plot(out.imp.RP979TTKSK.cim, show.marker.names=F, ylim=c(0,50), ylab="LOD", chr=c("1B", "2B", "7B"))
plot(out.imp.RP979TRTTF.cim, show.marker.names=F, col=2, add=T, chr=c("1B", "2B", "7B"))
plot(out.imp.RP979TMLKC.cim, show.marker.names=F, col=4, add=T, chr=c("1B", "2B", "7B"))
abline(h=3.99, col=3, lty=2)
legend("topright", c("TTKSK", "TRTTF","TMLKC"), col = c(1,2,4), lty=c(1,1,1), text.col = 1, cex=1.2, bg = 0)
dev.off()

##estimate 1-lod interval
lodint(out.imp.RP979TTKSK.cim,chr="1B",drop=1)
bayesint(out.imp.RP979TTKSK.cim, chr="1B",prob=0.99)

lodint(out.imp.RP979TRTTF.cim,chr="1B",drop=1)
bayesint(out.imp.RP979TRTTF.cim, chr="1B",prob=0.99)
lodint(out.imp.RP979TRTTF.cim,chr="2B",drop=1)
bayesint(out.imp.RP979TRTTF.cim, chr="2B",prob=0.99)

lodint(out.imp.RP979TMLKC.cim,chr="1B",drop=1)
bayesint(out.imp.RP979TMLKC.cim, chr="1B",prob=0.99)
lodint(out.imp.RP979TMLKC.cim,chr="7B",drop=1)
bayesint(out.imp.RP979TMLKC.cim, chr="7B",prob=0.99)

##estimated QTL effects
max(out.imp.RP979TTKSK.cim)
mar1=find.marker(RP979,chr="1B",pos=52.876)
mar1
plotPXG(RP979,pheno.col=2,marker=mar1)
effectplot(RP979,pheno.col=2,mname1=mar1, ylab="TTKSK (score)" )

max(out.imp.RP979TRTTF.cim)
mar2=find.marker(RP979,chr="1B",pos=52.876)
mar2
plotPXG(RP979,pheno.col=3,marker=mar2)
effectplot(RP979,pheno.col=3,mname1=mar2, ylab="TRTTF (score)" )

mar3=find.marker(RP979,chr="2B",pos=116.831)
mar3
plotPXG(RP979,pheno.col=3,marker=mar3)
effectplot(RP979,pheno.col=3,mname1=mar3, ylab="TRTTF (score)" )

max(out.imp.RP979TMLKC.cim)
mar4=find.marker(RP979,chr="1B",pos=52.876)
mar4
plotPXG(RP979,pheno.col=4,marker=mar4)
effectplot(RP979,pheno.col=4,mname1=mar4, ylab="TMLKC (score)" )

max(out.imp.RP979TMLKC.cim)
mar5=find.marker(RP979,chr="7B",pos=186)
mar5
plotPXG(RP979,pheno.col=4,marker=mar5)
effectplot(RP979,pheno.col=4,mname1=mar5, ylab="TMLKC (score)" )

#################################################
