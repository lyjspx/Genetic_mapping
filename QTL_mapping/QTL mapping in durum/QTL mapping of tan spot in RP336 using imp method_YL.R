#Set working directory under Mac
setwd("/Users/xuehui.li/Desktop/XL/research/crops/durum_wheat/p6_genome wide nested association mapping of tan spot in durum/RP336")

#load packages "qtl"
RP336=read.cross("csvr", "/Users/desert/OneDrive - North Dakota University System/Research/Tan_Spot_Project", "Final_map_RP336.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "riself")

str(RP336)
nind(RP336)
nchr(RP336)
totmar(RP336)
nmar(RP336)
plot(RP336)
plotMap(RP336,main = "Genetic map of RP336")
plotMissing(RP336)
summaryMap(RP336)

write.table(summary.map(RP336),file = "temp.xls",row.names = TRUE,col.names = TRUE,sep = "\t")
#############################################
#Normality test, Kolmogorov-Smirnov (K-S) normality test and Shapiro-Wilk's test
shapiro.test(RP336.pheno$X86.124R1)
shapiro.test(RP336.pheno$X86.124R2)
shapiro.test(RP336.pheno$X86.124R3)
shapiro.test(RP336.pheno$X86.124_knock_out)
shapiro.test(RP336.pheno$X86.124_knock_out_rep2)

#############################################
#Homogeneity of variance test
#load car package
with(RP336.pheno,leveneTest(X86.124R1,X86.124R3))
with(RP336.pheno,leveneTest(X86.124R2,X86.124R3))
with(RP336.pheno,leveneTest(X86.124_knock_out,X86.124_knock_out_rep2))
#passed homogeneity test


##########################################################
#Spearman's correlations among the phenotypes
str(RP336$pheno)
RP336.pheno=RP336$pheno[complete.cases(RP336$pheno[,2:8]),]
str(RP336.pheno)
mean(RP336.pheno[,2])
min(RP336.pheno[,2])
max(RP336.pheno[,2])

mean(RP336.pheno[,3])
min(RP336.pheno[,3])
max(RP336.pheno[,3])

mean(RP336.pheno[,4])
min(RP336.pheno[,4])
max(RP336.pheno[,4])

mean(RP336.pheno[,5])
min(RP336.pheno[,5])
max(RP336.pheno[,5])

mean(RP336.pheno[,8])
min(RP336.pheno[,8])
max(RP336.pheno[,8])


cor(RP336.pheno[,3], RP336.pheno[,6], method="pearson")
cor.test(RP336.pheno[,5], RP336.pheno[,8], method="pearson")

cor(RP336.pheno[,2], RP336.pheno[,4], method="pearson")
cor.test(RP336.pheno[,2], RP336.pheno[,4], method="pearson")

cor(RP336.pheno[,3], RP336.pheno[,4], method="pearson")
cor.test(RP336.pheno[,3], RP336.pheno[,4], method="pearson")
boxplot(RP336.pheno[,2:4], ylab="Score")

hist(RP336$pheno$X86.124mean,breaks = 4, labels = T,ylim = c(0,120),xlim = c(0,5),right = T,main = "RP336 with 86-124",xlab = "lesion type")###left open, right close
hist(RP336$pheno$knock_out_mean,breaks = 4, labels = T,ylim = c(0,80),xlim = c(0,5),right = T,main = "RP336 with 86-124 ToxA knockout",xlab = "lesion type")

t.test(RP336$pheno$X86.124mean,RP336.pheno$knock_out_mean,alternative = "greater")

##########################################################
#################
#plot genetic map of RP336
tiff(filename="Fig. Genetic map of RP336.tiff", width=8, height=6,units='in',res=300)
plotMap(RP336, main="Genetic map of RP336" )
dev.off()

###################################################################
#single QTL genome scan
#calculate conditional probability if using "hk" method
RP336.prob <- calc.genoprob(RP336, step=1, error.prob=0.01)
str(RP336.prob)

#simulate missing genotype if using "imp" method
RP336=jittermap(RP336, amount=1e-6)
RP336.sim <- sim.geno(RP336, step=2, n.draws=128, err=0.001)
str(RP336.sim)

##1.1 single-QTL genome scan for 86
##scanone wiht method "imp"
out.imp.86R1 <- scanone(RP336.sim, method="imp", pheno.col=2)
summary(out.imp.86R1)
plot(out.imp.86R1)
operm.imp.86R1 <- scanone(RP336.sim, pheno.col=2, method="imp", n.perm=1000)

out.imp.86R2 <- scanone(RP336.sim, method="imp", pheno.col=3)
summary(out.imp.86R2)
plot(out.imp.86R2)
operm.imp.86R2 <- scanone(RP336.sim, pheno.col=3, method="imp", n.perm=1000)

out.imp.86R3 <- scanone(RP336.sim, method="imp", pheno.col=4)
summary(out.imp.86R3)
plot(out.imp.86R3)
operm.imp.86R3 <- scanone(RP336.sim, pheno.col=4, method="imp", n.perm=1000)

out.imp.86mean <- scanone(RP336.sim, method="imp", pheno.col=5)
summary(out.imp.86mean)
plot(out.imp.86mean)
operm.imp.86mean <- scanone(RP336.sim, pheno.col=5, method="imp", n.perm=1000) #3.27
summary(operm.imp.86mean)

out.imp.86.knockR1 <- scanone(RP336.sim, method="imp", pheno.col=6)
summary(out.imp.86.knockR1)
plot(out.imp.86.knockR1)
operm.imp.86knockR1 <- scanone(RP336.sim, pheno.col=6, method="imp", n.perm=1000)

out.imp.86.knockR2 <- scanone(RP336.sim, method="imp", pheno.col=7)
summary(out.imp.86.knockR2)
plot(out.imp.86.knockR2)
operm.imp.86knockR2 <- scanone(RP336.sim, pheno.col=7, method="imp", n.perm=1000)

out.imp.86.knockmean <- scanone(RP336.sim, method="imp", pheno.col=8)
summary(out.imp.86.knockmean )
plot(out.imp.86.knockmean)
operm.imp.86knockmean <- scanone(RP336.sim, pheno.col=8, method="imp", n.perm=1000) #3.34
summary(operm.imp.86knockmean)

######
#plot two QTL in one plot
plot(out.imp.86mean,out.imp.86.knockmean,col = c("gray0","blue"), 
     main="RP336, 86-124 (Race 2 ToxA only) and 86124 ToxA knockedout")
text(900,16, "Black: 86124 Blue: 86124 ToxA knockedout")


## 1.5-lod interval
# lodint(out.imp.86R2,chr="5A",drop=1.5)
# lodint(out.imp.86mean,chr="5A",drop=1.5)

#QTL effect
# max(out.imp.86)
# mar=find.marker(RP336.sim, chr="5A", pos=204.5)
# plotPXG(RP336.sim, pheno.col=3, marker=mar)
# effectplot(RP336.sim, pheno.col=3, mname1=mar)

#two-QTL genome scan
##1.2 two-QTL genome scan with "scantwo" for 86
out2.imp.86mean <- scantwo(RP336.sim, pheno.col=5, method="imp")
plot(out2.imp.86mean)

operm2.imp.86mean <- scantwo(RP336.sim, pheno.col=5, method="imp", n.perm=100)
str(operm2.imp.86mean)
summary(operm2.imp.86mean, alpha=0.01)
summary(operm2.imp.86mean, alpha=0.05)
summary(out2.imp.86mean, thresholds=c(6.77,5.28, 4.37, 5.58, 2.7), what="int") #no pairs

out2.imp.86knockmean <- scantwo(RP336.sim, pheno.col=8, method="imp")
plot(out2.imp.86knockmean)

operm2.imp.86knockmean <- scantwo(RP336.sim, pheno.col=8, method="imp", n.perm=100)
str(operm2.imp.86knockmean)
summary(operm2.imp.86knockmean, alpha=0.01)
summary(operm2.imp.86knockmean, alpha=0.05)
summary(out2.imp.86mean, thresholds=c(6.52, 4.85, 4.11, 5.27, 2.67), what="int") #no pair

##multiple QTL mapping and forward/backward selection using "imp" method


#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans
calc.penalties(operm2.imp.86mean, alpha=0.05)
calc.penalties(operm2.imp.86knockmean, alpha=0.05)
##1.4 stepwise selection for 86124mean
chr.86 <- c("5A")
pos.86 <- c(238)
qtl.86 <- makeqtl(RP336.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1
stepout.RP336.86mean=stepwiseqtl(RP336.sim, pheno.col=5, additive.only=FALSE, penalties=c(2.93, 4.34, 2.15), qtl = qtl.86, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP336.86mean

chr.86 <- c("2B","5A","5A","5B")
pos.86 <- c(14.216,86,238.408,140)
qtl.86 <- makeqtl(RP336.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1+Q2+Q3+Q4
out.fq.86mean <- fitqtl(RP336.sim, pheno.col=5, qtl=qtl.86, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.86mean

rqtl.86mean =refineqtl(RP336.sim, pheno.col=5, qtl=qtl.86, method="imp", formula=my.formula)
rqtl.86mean
par(mfrow=c(1,1))
plotLodProfile(rqtl.86mean, col=2, qtl.labels = FALSE, lwd=1)

plot(out.imp.86mean,chr = c("2B","5A","5B"),col=4, ylim=c(0,15), ylab="LOD", lwd=1,main="QTL mapping of 86-124 isolate on RP336")
plotLodProfile(rqtl.86mean, col=2, add=T, qtl.labels = TRUE, lwd=1)
abline(h=3.32, col=3, lty=2)
text(450,4,"1000 permutation at p=0.05, threshold = 3.32",col = 3)
legend(340,15,c("Standard Interval Mapping","Multiple QTL Method"),col=c(4,2),lwd=1)

mar.5A.86=find.marker(RP336.sim, chr="5A",pos=238)
mar.5A.86

find.marker(RP336.sim,chr = "5A", pos = 86)

lodint(rqtl.86mean,qtl.index=1, drop=1.5)
lodint(rqtl.86mean,qtl.index=2, drop=1.5)
lodint(rqtl.86mean,qtl.index=3, drop=1.5)
lodint(rqtl.86mean,qtl.index=4, drop=1.5)

out.fq.86 <- fitqtl(RP336.sim, pheno.col=2, qtl=rqtl.86, method="imp", formula=my.formula, get.ests=TRUE)
out.fq.86

tiff(filename="Fig. QTL effect for 86.tiff", width=6, height=6,units='in',res=300)
effectplot(RP336.sim, pheno.col=5, mname=mar.5A.86)
dev.off()

##################################################
##estimated QTL effects
mar_86_1=find.marker(RP336,chr="2B",pos=14.2)
mar_86_1
plotPXG(RP336,pheno.col=5,marker=mar_86_1)
mar_86_1_effect =effectplot(RP336,pheno.col=5,mname1=mar_86_1, ylab="86 (score)" )
mar_86_1_effect$Means

mar_86_2=find.marker(RP336,chr="5A",pos=86)
mar_86_2
plotPXG(RP336,pheno.col=5,marker=mar_86_2)
mar_86_2_effect =effectplot(RP336,pheno.col=5,mname1=mar_86_2, ylab="86 (score)" )
mar_86_2_effect$Means

mar_86_3=find.marker(RP336,chr="5A",pos=238.4)
mar_86_3
plotPXG(RP336,pheno.col=5,marker=mar_86_3)
mar_86_3_effect =effectplot(RP336,pheno.col=5,mname1=mar_86_3, ylab="86 (score)" )
mar_86_3_effect$Means

mar_86_4=find.marker(RP336,chr="5B",pos=140)
mar_86_4
x = plotPXG(RP336,pheno.col=5,marker=mar_86_4)
mar_86_4_effect =effectplot(RP336,pheno.col=5,mname1=mar_86_4, ylab="86 (score)" )
mar_86_4_effect$Means
#################################################



##1.4 stepwise selection for 86124 knock mean
chr.86 <- c("5A")
pos.86 <- c(236)
qtl.86<- makeqtl(RP336.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1
stepout.RP336.86knockmean=stepwiseqtl(RP336.sim, pheno.col=8, additive.only=FALSE, penalties=c(3.08, 4.29, 1.89), qtl = qtl.86, formula=y~Q1, max.qtl=6, verbose=FALSE, method="imp", keeptrace = TRUE)
stepout.RP336.86knockmean

chr.86 <- c("1B","3A","5A","7A")
pos.86 <- c(90.457,13.551,236.000,3.536)
qtl.86 <- makeqtl(RP336.sim, chr.86, pos.86)
qtl.86
my.formula=y~Q1+Q2+Q3+Q4
out.fq.86knockmean <- fitqtl(RP336.sim, pheno.col=8, qtl=qtl.86, method="imp", drop=TRUE, formula=my.formula, get.ests=TRUE)
out.fq.86knockmean

rqtl.86knockmean =refineqtl(RP336.sim, pheno.col=8, qtl=qtl.86, method="imp", formula=my.formula)
rqtl.86knockmean
par(mfrow=c(1,1))
plotLodProfile(rqtl.86knockmean, col=2, qtl.labels = FALSE, lwd=1)


plot(out.imp.86.knockmean,chr = c("1B","3A","5A","7A"),col=4, ylim=c(0,25), ylab="LOD", lwd=1,main="QTL mapping of 86-124 knock out on RP336")
plotLodProfile(rqtl.86knockmean, col=2, add=T, qtl.labels = TRUE, lwd=1)
abline(h=3.32, col=3, lty=2)
text(450,4,"1000 permutation at p=0.05, threshold = 3.32",col = 3)
legend(35,25,c("Standard Interval Mapping","Multiple QTL Method"),col=c(4,2),lwd=1)

dev.off()

mar.5A.86=find.marker(RP336.sim, chr="5A",pos=236)
mar.5A.86

bayesint(rqtl.86knockmean, qtl.index=1, prob=0.99)
bayesint(rqtl.86knockmean, qtl.index=2, prob=0.99)
bayesint(rqtl.86knockmean, qtl.index=3, prob=0.99)
bayesint(rqtl.86knockmean, qtl.index=4, prob=0.99)

lodint(rqtl.86knockmean,qtl.index=1, drop=1.5)
lodint(rqtl.86knockmean,qtl.index=2, drop=1.5)
lodint(rqtl.86knockmean,qtl.index=3, drop=1.5)
lodint(rqtl.86knockmean,qtl.index=4, drop=1.5)

find.marker(RP336.sim,chr="3A",pos = 13.6)
find.marker(RP336.sim,chr="5A",pos = 236)


tiff(filename="Fig. QTL effect for 86.tiff", width=6, height=6,units='in',res=300)
effectplot(RP336.sim, pheno.col=8, mname1=mar.5A.86)
dev.off()

##################################################
##estimated QTL effects
mar_86_knock_1=find.marker(RP336,chr="1B",pos=90.5)
mar_86_knock_1
plotPXG(RP336,pheno.col=8,marker=mar_86_knock_1)
mar_86_knock_1_effect =effectplot(RP336,pheno.col=8,mname1=mar_86_knock_1, ylab="86 (score)" )
mar_86_knock_1_effect$Means

mar_86_knock_2=find.marker(RP336,chr="3A",pos=13.6)
mar_86_knock_2
plotPXG(RP336,pheno.col=8,marker=mar_86_knock_2)
mar_86_knock_2_effect =effectplot(RP336,pheno.col=8,mname1=mar_86_knock_2, ylab="86 (score)" )
mar_86_knock_2_effect$Means

mar_86_knock_3=find.marker(RP336,chr="5A",pos=236)
mar_86_knock_3
plotPXG(RP336,pheno.col=8,marker=mar_86_knock_3)
mar_86_knock_3_effect =effectplot(RP336,pheno.col=8,mname1=mar_86_knock_3, ylab="86 (score)" )
mar_86_knock_3_effect$Means

mar_86_knock_4=find.marker(RP336,chr="7A",pos=3.5)
mar_86_knock_4
plotPXG(RP336,pheno.col=8,marker=mar_86_knock_4)
mar_86_knock_4_effect =effectplot(RP336,pheno.col=8,mname1=mar_86_knock_4, ylab="86 (score)" )
mar_86_knock_4_effect$Means
