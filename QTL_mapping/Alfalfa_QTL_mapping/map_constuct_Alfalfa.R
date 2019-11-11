
#setwd
currentwd = getwd()
setwd(currentwd)
setwd("F:/alfalfa_6_r")
#load ASMap and qtl package

Alfalfa = read.cross(format=c("csvr"),"F:/alfalfa_6_r","data_for_r.csv" , na.strings=c("-"),genotypes = c("A","H","B"),estimate.map = FALSE, crosstype = "bc")

str(Alfalfa)
nind(Alfalfa)
nchr(Alfalfa)
totmar(Alfalfa)
nmar(Alfalfa)
plot(Alfalfa)
plotMap(Alfalfa)
plotMissing(Alfalfa)

#remove individuals with too many missings
par(mfrow=c(1,1), las = 1)
plot(ntyped(Alfalfa,what=c("mar")))
Alfalfa = subset(Alfalfa, ind = (ntyped(Alfalfa)> 7000))

#remove markers with too many missings
plot(ntyped(Alfalfa, "mar"), xlab="Marker", ylab="No. typed samples",main="No. typed samples by markers")
num_bymar = ntyped(Alfalfa,"mar")
todrop = names(num_bymar[num_bymar < 180])
Alfalfa = drop.markers(Alfalfa,todrop)

str(Alfalfa)
nind(Alfalfa)
nchr(Alfalfa)
totmar(Alfalfa)

cg = comparegeno(Alfalfa)
str(cg)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes") 
rug(cg[lower.tri(cg)])
wh <- which(cg > 0.95, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh
Genotype = pull.geno(Alfalfa)
str(Genotype)
dim(Genotype)
Genotype[1:10,1:5]
table(Genotype[7,], Genotype[106,])

gc = genClones(Alfalfa, tol = 0.85)
gc$cgd
cgd <- gc$cgd
Alfalfa <- fixClones(Alfalfa, cgd, consensus = TRUE)

gt <- geno.table(Alfalfa)
str(gt)
dim(gt)
gt[1:5,]

todrop1 <- rownames(gt[gt$P.value < 1e-10,])
length(todrop1)
todrop2=rownames(gt[gt$P.value < 0.05/totmar(Alfalfa),])
length(todrop2)


##estimate minor allele frequency

MAF=rep(NA,13550)
for(i in 1:13550){
  MAF[i]=min()
}



