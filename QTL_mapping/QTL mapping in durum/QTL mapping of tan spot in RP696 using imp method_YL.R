#Set working directory under Mac
setwd("/Users/xuehui.li/Desktop/XL/research/crops/durum_wheat/p6_genome wide nested association mapping of tan spot in durum/RP696")

#load packages "qtl"
RP696=read.cross("csvr", "/Users/desert/OneDrive - North Dakota University System/Research/Tan_Spot_Project", "Final_map_RP696_GBS.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "riself")

str(RP696)
nind(RP696)
nchr(RP696)
totmar(RP696)
nmar(RP696)
plot(RP696)
plotMap(RP696,main = "Genetic map of RP696")
plotMissing(RP696)
summaryMap(RP696)

hist(RP696$pheno$X86124mean,breaks = 4, labels = T,ylim = c(0,100),xlim = c(0,5),right = T,main = "RP696 with 86-124",xlab = "lesion type")###left open, right close
#write.table(summary.map(RP696),file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")

##########################################################
#Spearman's correlations among the phenotypes
str(RP696$pheno)
RP696.pheno=RP696$pheno[complete.cases(RP696$pheno[,2:5]),]
str(RP696.pheno)
mean(RP696.pheno[,2])
min(RP696.pheno[,2])
max(RP696.pheno[,2])

mean(RP696.pheno[,3])
min(RP696.pheno[,3])
max(RP696.pheno[,3])

mean(RP696.pheno[,4])
min(RP696.pheno[,4])
max(RP696.pheno[,4])

mean(RP696.pheno[,5])
min(RP696.pheno[,5])
max(RP696.pheno[,5])

cor(RP696.pheno[,3], RP696.pheno[,6], method="pearson")
cor.test(RP696.pheno[,5], RP696.pheno[,8], method="pearson")

cor(RP696.pheno[,2], RP696.pheno[,4], method="pearson")
cor.test(RP696.pheno[,2], RP696.pheno[,4], method="pearson")

cor(RP696.pheno[,3], RP696.pheno[,4], method="pearson")
cor.test(RP696.pheno[,3], RP696.pheno[,4], method="pearson")
boxplot(RP696.pheno[,2:4], ylab="Score")

##########################################################
#################
#plot genetic map of RP696
tiff(filename="Fig. Genetic map of RP696.tiff", width=8, height=6,units='in',res=300)
plotMap(RP696, main="Genetic map of RP696" )
dev.off()

###################################################################
#single QTL genome scan
#calculate conditional probability if using "hk" method
RP696.prob <- calc.genoprob(RP696, step=1, error.prob=0.01)
str(RP696.prob)

#simulate missing genotype if using "imp" method
RP696=jittermap(RP696, amount=1e-6)
RP696.sim <- sim.geno(RP696, step=2, n.draws=128, err=0.001)
str(RP696.sim)

##1.1 single-QTL genome scan for 86
##scanone wiht method "imp"
out.imp.86R1 <- scanone(RP696.sim, method="imp", pheno.col=2)
summary(out.imp.86R1)
plot(out.imp.86R1)
operm.imp.86R1 <- scanone(RP696.sim, pheno.col=2, method="imp", n.perm=1000)

out.imp.86R2 <- scanone(RP696.sim, method="imp", pheno.col=3)
summary(out.imp.86R2)
plot(out.imp.86R2, main = "RP696 Rep2")
operm.imp.86R2 <- scanone(RP696.sim, pheno.col=3, method="imp", n.perm=1000)

out.imp.86R3 <- scanone(RP696.sim, method="imp", pheno.col=4)
summary(out.imp.86R3)
plot(out.imp.86R3,main = "RP696 Rep3")
operm.imp.86R3 <- scanone(RP696.sim, pheno.col=4, method="imp", n.perm=1000)

out.imp.86mean <- scanone(RP696.sim, method="imp", pheno.col=5)
summary(out.imp.86mean)
plot(out.imp.86mean, main = "RP696 mean")
operm.imp.86mean <- scanone(RP696.sim, pheno.col=5, method="imp", n.perm=1000)
summary(operm.imp.86mean)

## 1-lod interval
lodint(out.imp.86R1,chr="5A",drop=1.5)
 lodint(out.imp.86R2,chr="5A",drop=1.5)
 lodint(out.imp.86mean,chr="5A",drop=1.5)

#QTL effect
# max(out.imp.86)
# mar=find.marker(RP696.sim, chr="5A", pos=204.5)
# plotPXG(RP696.sim, pheno.col=3, marker=mar)
# effectplot(RP696.sim, pheno.col=3, mname1=mar)

#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for 86
out2.imp.86R2 <- scantwo(RP696.sim, pheno.col=3, method="imp")
plot(out2.imp.86R2)
operm2.imp.86R2 <- scantwo(RP696.sim, pheno.col=3, method="imp", n.perm=100)
str(operm2.imp.86R2)
plot(operm2.imp.86R2)
summary(operm2.imp.86R2, alpha=0.01)
summary(operm2.imp.86R2, alpha=0.05)
summary(out2.imp.86R2, thresholds=c(6.7,4.96, 4.29, 5.24, 2.8), what="int")


out2.imp.86mean <- scantwo(RP696.sim, pheno.col=5, method="imp")
plot(out2.imp.86mean)
operm2.imp.86mean <- scantwo(RP696.sim, pheno.col=5, method="imp", n.perm=100)
str(operm2.imp.86mean)
summary(operm2.imp.86mean, alpha=0.01)
summary(operm2.imp.86mean, alpha=0.05)
summary(out2.imp.86mean, thresholds=c(6.29,4.93, 4.31, 4.92, 2.67), what="int") #No pairs of loci


##multiple QTL mapping and forward/backward selection using "imp" method


#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.86mean, alpha=0.05)
#2.97 4.31 1.96

##1.4 stepwise selection for 86124mean
chr.86 <- c("3A")
pos.86 <- c(8.54)
qtl.86 <- makeqtl(RP696.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1
stepout.RP696.86mean=stepwiseqtl(RP696.sim, pheno.col=5, additive.only=FALSE, penalties=c(2.97, 4.31, 1.95), qtl = qtl.86, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP696.86mean

chr.86 <- c("2B","3A","5A","5A")
pos.86 <- c(10, 7.945,53.099,193.559)
qtl.86 <- makeqtl(RP696.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1+Q2+Q3+Q4
out.fq.86mean <- fitqtl(RP696.sim, pheno.col=5, qtl=qtl.86, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.86mean

rqtl.86mean =refineqtl(RP696.sim, pheno.col=5, qtl=qtl.86, method="imp", formula=my.formula)
rqtl.86mean

plot(out.imp.86mean,chr = c("2B","3A","5A"),col=4, ylim=c(0,15), ylab="LOD", lwd=1,main="QTL mapping of 86-124 isolate on RP696")
plotLodProfile(rqtl.86mean, col=2, add=T, qtl.labels = TRUE, lwd=1,show.marker.names = T)
abline(h=3.31, col=3, lty=2)
text(350,4,"1000 permutation at p=0.05, threshold = 3.31",col = 3)
legend(340,15,c("Standard Interval Mapping","Multiple QTL Method"),col=c(4,2),lwd=1)

mar.5A.86=find.marker(RP696.sim, chr="5A",pos=238)
mar.5A.86

bayesint(rqtl.86, qtl.index=1, prob=0.99)
bayesint(rqtl.86, qtl.index=2, prob=0.99)
bayesint(rqtl.86, qtl.index=3, prob=0.99)
bayesint(rqtl.86, qtl.index=4, prob=0.99)

lodint(rqtl.86mean,qtl.index=1, drop=1.5)
lodint(rqtl.86mean,qtl.index=2, drop=1.5)
lodint(rqtl.86mean,qtl.index=3, drop=1.5)
lodint(rqtl.86mean,qtl.index=4, drop=1.5)

out.fq.86 <- fitqtl(RP696.sim, pheno.col=2, qtl=rqtl.86, method="imp", formula=my.formula, get.ests=TRUE)
out.fq.86

tiff(filename="Fig. QTL effect for 86.tiff", width=6, height=6,units='in',res=300)
effectplot(RP696.sim, pheno.col=5, mname=mar.5A.86)
dev.off()

#################################################

##estimated QTL effects
mar_86_1=find.marker(RP696,chr="2B",pos=10)
mar_86_1
plotPXG(RP696,pheno.col=5,marker=mar_86_1)
mar_86_1_effect =effectplot(RP696,pheno.col=5,mname1=mar_86_1, ylab="86 (score)" )
mar_86_1_effect$Means

mar_86_2=find.marker(RP696,chr="3A",pos=7.9)
mar_86_2
plotPXG(RP696,pheno.col=5,marker=mar_86_2)
mar_86_2_effect =effectplot(RP696,pheno.col=5,mname1=mar_86_2, ylab="86 (score)" )
mar_86_2_effect$Means

mar_86_3=find.marker(RP696,chr="5A",pos=53.1)
mar_86_3
plotPXG(RP696,pheno.col=5,marker=mar_86_3)
mar_86_3_effect =effectplot(RP696,pheno.col=5,mname1=mar_86_3, ylab="86 (score)" )
mar_86_3_effect$Means

mar_86_4=find.marker(RP696,chr="5A",pos=193.6)
mar_86_4
x = plotPXG(RP696,pheno.col=5,marker=mar_86_4)
mar_86_4_effect =effectplot(RP696,pheno.col=5,mname1=mar_86_4, ylab="86 (score)" )
mar_86_4_effect$Means
