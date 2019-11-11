#Set working directory under Mac
setwd("/Users/xuehui.li/Desktop/XL/research/crops/durum_wheat/p6_genome wide nested association mapping of tan spot in durum/RP979")

#load packages "qtl"
RP979=read.cross("csvr", "/Users/desert/OneDrive - North Dakota University System/Research/Tan_Spot_Project", "Final_map_RP979.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "riself")

str(RP979)
nind(RP979)
nchr(RP979)
totmar(RP979)
nmar(RP979)
plot(RP979)
plotMap(RP979,main = "Genetic map of RP979")
plotMissing(RP979)
summaryMap(RP979)

hist(RP979$pheno$X86124Mean,breaks = 3, labels = T,ylim = c(0,150),xlim = c(0,5),right = T,main = "RP979 with 86-124",xlab = "lesion type")###left open, right close


#write.table(summary.map(RP979),file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")

##########################################################
#Spearman's correlations among the phenotypes
str(RP979$pheno)
RP979.pheno=RP979$pheno[complete.cases(RP979$pheno[,2:5]),]
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

mean(RP979.pheno[,5])
min(RP979.pheno[,5])
max(RP979.pheno[,5])

cor(RP979.pheno[,3], RP979.pheno[,6], method="pearson")
cor.test(RP979.pheno[,5], RP979.pheno[,8], method="pearson")

cor(RP979.pheno[,2], RP979.pheno[,4], method="pearson")
cor.test(RP979.pheno[,2], RP979.pheno[,4], method="pearson")

cor(RP979.pheno[,3], RP979.pheno[,4], method="pearson")
cor.test(RP979.pheno[,3], RP979.pheno[,4], method="pearson")
boxplot(RP979.pheno[,2:4], ylab="Score")

##########################################################
#################
#plot genetic map of RP979
tiff(filename="Fig. Genetic map of RP979.tiff", width=8, height=6,units='in',res=300)
plotMap(RP979, main="Genetic map of RP979" )
dev.off()

###################################################################
#single QTL genome scan
#calculate conditional probability if using "hk" method
RP979.prob <- calc.genoprob(RP979, step=1, error.prob=0.01)
str(RP979.prob)

#simulate missing genotype if using "imp" method
RP979=jittermap(RP979, amount=1e-6)
RP979.sim <- sim.geno(RP979, step=2, n.draws=128, err=0.001)
str(RP979.sim)

##1.1 single-QTL genome scan for 86
##scanone wiht method "imp"
out.imp.86R1 <- scanone(RP979.sim, method="imp", pheno.col=2)
summary(out.imp.86R1)
plot(out.imp.86R1)
operm.imp.86R1 <- scanone(RP979.sim, pheno.col=2, method="imp", n.perm=1000)

out.imp.86R2 <- scanone(RP979.sim, method="imp", pheno.col=3)
summary(out.imp.86R2)
plot(out.imp.86R2)
operm.imp.86R2 <- scanone(RP979.sim, pheno.col=3, method="imp", n.perm=1000)

out.imp.86R3 <- scanone(RP979.sim, method="imp", pheno.col=4)
summary(out.imp.86R3)
plot(out.imp.86R3)
operm.imp.86R3 <- scanone(RP979.sim, pheno.col=4, method="imp", n.perm=1000)

out.imp.86mean <- scanone(RP979.sim, method="imp", pheno.col=5)
summary(out.imp.86mean)
plot(out.imp.86mean, main = "RP979")
operm.imp.86mean <- scanone(RP979.sim, pheno.col=5, method="imp", n.perm=1000)
summary(operm.imp.86mean)


## 1-lod interval
# lodint(out.imp.86R2,chr="5A",drop=1.5)
# lodint(out.imp.86mean,chr="5A",drop=1.5)

#QTL effect
# max(out.imp.86)
# mar=find.marker(RP979.sim, chr="5A", pos=204.5)
# plotPXG(RP979.sim, pheno.col=3, marker=mar)
# effectplot(RP979.sim, pheno.col=3, mname1=mar)

#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for 86
out2.imp.86mean <- scantwo(RP979.sim, pheno.col=5, method="imp")
plot(out2.imp.86mean)

operm2.imp.86mean <- scantwo(RP979.sim, pheno.col=5, method="imp", n.perm=100)
str(operm2.imp.86mean)
summary(operm2.imp.86mean, alpha=0.01)
summary(operm2.imp.86mean, alpha=0.05)
summary(out2.imp.86mean, thresholds=c(7.02,5.18, 4.32, 5.68, 3.21), what="int") #no pairs


##multiple QTL mapping and forward/backward selection using "imp" method


#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.86mean, alpha=0.05)

##1.4 stepwise selection for 86124mean
chr.86 <- c("5B")
pos.86 <- c(18)
qtl.86 <- makeqtl(RP979.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1
stepout.RP979.86mean=stepwiseqtl(RP979.sim, pheno.col=5, additive.only=FALSE, penalties=c(3.08, 4.29, 1.89), qtl = qtl.86, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP979.86mean

chr.86 <- c("1A","3B","5B","5B","7A","7A")
pos.86 <- c(141.026,44.691,19.181,154.404,82.000,106.820)
qtl.86 <- makeqtl(RP979.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1+Q2+Q3+Q4+Q5
out.fq.86mean <- fitqtl(RP979.sim, pheno.col=5, qtl=qtl.86, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.86mean

rqtl.86mean =refineqtl(RP979.sim, pheno.col=5, qtl=qtl.86, method="imp", formula=my.formula)
rqtl.86mean
par(mfrow=c(1,1))
plotLodProfile(rqtl.86mean, col=2, qtl.labels = FALSE, lwd=1)

plot(out.imp.86mean,chr = c("1A","3B","5B","7A"),col=4, ylim=c(0,15), ylab="LOD", lwd=1,main="QTL mapping of 86-124 isolate on RP979")
plotLodProfile(rqtl.86mean, col=2, add=T, qtl.labels = TRUE, lwd=1)
abline(h=3.6, col=3, lty=2)
text(180,4,"1000 permutation at p=0.05, threshold = 3.6",col = 3)
legend(340,15,c("Standard Interval Mapping","Multiple QTL Method"),col=c(4,2),lwd=1)

mar.5A.86=find.marker(RP979.sim, chr="5A",pos=238)
mar.5A.86

bayesint(rqtl.86, qtl.index=1, prob=0.99)
bayesint(rqtl.86, qtl.index=2, prob=0.99)
bayesint(rqtl.86, qtl.index=3, prob=0.99)
bayesint(rqtl.86, qtl.index=4, prob=0.99)

lodint(rqtl.86,qtl.index=1, drop=1.5)
lodint(rqtl.86,qtl.index=2, drop=1.5)
lodint(rqtl.86,qtl.index=3, drop=1.5)
lodint(rqtl.86,qtl.index=4, drop=1.5)

out.fq.86 <- fitqtl(RP979.sim, pheno.col=2, qtl=rqtl.86, method="imp", formula=my.formula, get.ests=TRUE)
out.fq.86

tiff(filename="Fig. QTL effect for 86.tiff", width=6, height=6,units='in',res=300)
effectplot(RP979.sim, pheno.col=5, mname=mar.5A.86)
dev.off()

##################################################
##estimated QTL effects
max(out.imp.RP97986.cim)
mar1=find.marker(RP979,chr="5A",pos=204.55)
mar1
plotPXG(RP979,pheno.col=5,marker=mar1)
effectplot(RP979,pheno.col=5,mname1=mar1, ylab="86 (score)" )

#################################################




##################################################
##estimated QTL effects
max(out.imp.RP97986.cim)
mar1=find.marker(RP979,chr="5A",pos=204.55)
mar1
plotPXG(RP979,pheno.col=5,marker=mar1)
effectplot(RP979,pheno.col=5,mname1=mar1, ylab="86 (score)" )