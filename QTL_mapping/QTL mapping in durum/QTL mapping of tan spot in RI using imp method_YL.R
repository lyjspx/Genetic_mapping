
#load packages "qtl"
RPIum=read.cross("csvr", "/Users/desert/OneDrive - North Dakota University System/Research/Tan_Spot_Project", "Final_map_RIumillo.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "riself")

str(RPIum)
nind(RPIum)
nchr(RPIum)
totmar(RPIum)
nmar(RPIum)
plot(RPIum)
plotMap(RPIum,main = "Genetic map of RPIum")
plotMissing(RPIum)
summaryMap(RPIum)


hist(RPIum$pheno$X86124Mean,breaks = 4, labels = T,ylim = c(0,90),xlim = c(0,5),right = T,main = "RIum with 86-124",xlab = "lesion type")###left open, right close

#write.table(summary.map(RPIum),file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")

##########################################################
#Spearman's correlations among the phenotypes
str(RPIum$pheno)
RPIum.pheno=RPIum$pheno[complete.cases(RPIum$pheno[,2:5]),]
str(RPIum.pheno)
mean(RPIum.pheno[,2])
min(RPIum.pheno[,2])
max(RPIum.pheno[,2])

mean(RPIum.pheno[,3])
min(RPIum.pheno[,3])
max(RPIum.pheno[,3])

mean(RPIum.pheno[,4])
min(RPIum.pheno[,4])
max(RPIum.pheno[,4])

mean(RPIum.pheno[,5])
min(RPIum.pheno[,5])
max(RPIum.pheno[,5])



cor(RPIum.pheno[,3], RPIum.pheno[,6], method="pearson")
cor.test(RPIum.pheno[,5], RPIum.pheno[,8], method="pearson")

cor(RPIum.pheno[,2], RPIum.pheno[,4], method="pearson")
cor.test(RPIum.pheno[,2], RPIum.pheno[,4], method="pearson")

cor(RPIum.pheno[,3], RPIum.pheno[,4], method="pearson")
cor.test(RPIum.pheno[,3], RPIum.pheno[,4], method="pearson")
boxplot(RPIum.pheno[,2:4], ylab="Score")

##########################################################
#################
#plot genetic map of RPIum
# tiff(filename="Fig. Genetic map of RPIum.tiff", width=8, height=6,units='in',res=300)
# plotMap(RPIum, main="Genetic map of RPIum" )
# dev.off()

###################################################################
#single QTL genome scan
#calculate conditional probability if using "hk" method
RPIum.prob <- calc.genoprob(RPIum, step=1, error.prob=0.01)
str(RPIum.prob)

#simulate missing genotype if using "imp" method
RPIum=jittermap(RPIum, amount=1e-6)
RPIum.sim <- sim.geno(RPIum, step=2, n.draws=128, err=0.001)
str(RPIum.sim)

##1.1 single-QTL genome scan for 86
##scanone wiht method "imp"
out.imp.86R1 <- scanone(RPIum.sim, method="imp", pheno.col=2)
summary(out.imp.86R1)
plot(out.imp.86R1)
operm.imp.86R1 <- scanone(RPIum.sim, pheno.col=2, method="imp", n.perm=1000)

out.imp.86R2 <- scanone(RPIum.sim, method="imp", pheno.col=3)
summary(out.imp.86R2)
plot(out.imp.86R2)
operm.imp.86R2 <- scanone(RPIum.sim, pheno.col=3, method="imp", n.perm=1000)

out.imp.86R3 <- scanone(RPIum.sim, method="imp", pheno.col=4)
summary(out.imp.86R3)
plot(out.imp.86R3)
operm.imp.86R3 <- scanone(RPIum.sim, pheno.col=4, method="imp", n.perm=1000)

out.imp.86mean <- scanone(RPIum.sim, method="imp", pheno.col=5)
summary(out.imp.86mean)
plot(out.imp.86mean, main = "RI")
  operm.imp.86mean <- scanone(RPIum.sim, pheno.col=5, method="imp", n.perm=1000) #3.32
summary(operm.imp.86mean)



## 1-lod interval
# lodint(out.imp.86R2,chr="5A",drop=1.5)
# lodint(out.imp.86mean,chr="5A",drop=1.5)

#QTL effect
# max(out.imp.86)
# mar=find.marker(RPIum.sim, chr="5A", pos=204.5)
# plotPXG(RPIum.sim, pheno.col=3, marker=mar)
# effectplot(RPIum.sim, pheno.col=3, mname1=mar)

#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for 86
out2.imp.86mean <- scantwo(RPIum.sim, pheno.col=5, method="imp")
plot(out2.imp.86mean)

operm2.imp.86mean <- scantwo(RPIum.sim, pheno.col=5, method="imp", n.perm=100)
str(operm2.imp.86mean)
summary(operm2.imp.86mean, alpha=0.01)
summary(operm2.imp.86mean, alpha=0.05)
summary(out2.imp.86mean, thresholds=c(6.44,5.27, 4.9,4.85, 2.47), what="int") #no pairs



##multiple QTL mapping and forward/backward selection using "imp" method


#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.86mean, alpha=0.05)
##1.4 stepwise selection for 86124mean
chr.86 <- c("3A")
pos.86 <- c(21.1)
qtl.86 <- makeqtl(RPIum.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1
stepout.RPIum.86mean=stepwiseqtl(RPIum.sim, pheno.col=5, additive.only=FALSE, penalties=c(2.8, 4.89, 2.43), qtl = qtl.86, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RPIum.86mean

chr.86 <- c("3A","3B")
pos.86 <- c(21.125,72.621)
qtl.86 <- makeqtl(RPIum.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1+Q2
out.fq.86mean <- fitqtl(RPIum.sim, pheno.col=5, qtl=qtl.86, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.86mean

rqtl.86mean =refineqtl(RPIum.sim, pheno.col=5, qtl=qtl.86, method="imp", formula=my.formula)
rqtl.86mean
par(mfrow=c(1,1))
plotLodProfile(rqtl.86mean, col=2, qtl.labels = T, lwd=1)

plot(out.imp.86mean,chr = c("3A"),col=4, ylim=c(0,15), ylab="LOD", lwd=1,main="QTL mapping of 86-124 isolate on RIum")
plotLodProfile(rqtl.86mean, col=2, add=T, qtl.labels = TRUE, lwd=1)
abline(h=3.16, col=3, lty=2)
text(180,4,"1000 permutation at p=0.05, threshold = 3.16",col = 3)
legend(150,15,c("Standard Interval Mapping","Multiple QTL Method"),col=c(4,2),lwd=1)

mar.5A.86=find.marker(RPIum.sim, chr="3A",pos=21.1)
mar.5A.86

bayesint(rqtl.86mean, qtl.index=1, prob=0.99)
bayesint(rqtl.86mean, qtl.index=2, prob=0.99)


lodint(rqtl.86mean,qtl.index=1, drop=1.5)
lodint(rqtl.86,qtl.index=2, drop=1.5)

out.fq.86 <- fitqtl(RPIum.sim, pheno.col=5, qtl=rqtl.86mean, method="imp", formula=my.formula, get.ests=TRUE)
out.fq.86

tiff(filename="Fig. QTL effect for 86.tiff", width=6, height=6,units='in',res=300)
effectplot(RPIum.sim, pheno.col=5, mname=mar.5A.86)
dev.off()

##################################################
##estimated QTL effects
mar1=find.marker(RPIum,chr="3A",pos=21.1)
mar1
plotPXG(RPIum,pheno.col=5,marker=mar1)
mar_effect = effectplot(RPIum,pheno.col=5,mname1=mar1, ylab="86 (score)" )
mar_effect$Means

#################################################


