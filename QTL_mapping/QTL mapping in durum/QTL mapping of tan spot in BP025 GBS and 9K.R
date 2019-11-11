library("qtl")
BP025=read.cross("csvr", "Y:/Yuan/integrated map/Ben_PI41025", "BP025_9K_GBS_Rqtl.csv", na.strings=c("-"), 
                 estimate.map = FALSE, crosstype = "riself")

str(BP025)
nind(BP025)
nchr(BP025)
totmar(BP025)
nmar(BP025)
plot(BP025)
plotMap(BP025,main = "Genetic map of BP025")
plotMissing(BP025)
summaryMap(BP025)

BP025.prob <- calc.genoprob(BP025, step=1, error.prob=0.01)
str(BP025.prob)

#simulate missing genotype if using "imp" method
BP025=jittermap(BP025, amount=1e-6)
BP025.sim <- sim.geno(BP025, step=2, n.draws=128, err=0.001)
str(BP025.sim)

################################
#Scanone

# out.imp.BP025L13 = scanone(BP025.sim,method = "imp",pheno.col = 2)
# summary(out.imp.BP025L13)
# plot(out.imp.BP025L13)
# operm.imp.BP025L13 = scanone(BP025.sim,pheno.col = 2,method = "imp",n.perm = 1000)

out.imp.BP0253319 = scanone(BP025.sim,method = "imp",pheno.col = 3)
summary(out.imp.BP0253319)
plot(out.imp.BP0253319)
operm.imp.BP0253319 = scanone(BP025.sim,pheno.col = 3,method = "imp",n.perm = 1000)
summary(out.imp.BP0253319,perms = operm.imp.BP0253319,alpha = 0.05)

out.imp.BP025DW5 = scanone(BP025.sim,method = "imp",pheno.col = 4)
summary(out.imp.BP025DW5)
plot(out.imp.BP025DW5)
operm.imp.BP025DW5 = scanone(BP025.sim,pheno.col = 4,method = "imp",n.perm = 1000)
summary(out.imp.BP025DW5,perms = operm.imp.BP025DW5,alpha = 0.05)

out.imp.BP025ToxA = scanone(BP025.sim,method = "imp",pheno.col = 5)
summary(out.imp.BP025ToxA)
plot(out.imp.BP025ToxA)
operm.imp.BP025ToxA = scanone(BP025.sim,pheno.col = 5,method = "imp",n.perm = 1000)

out.imp.BP02586124 = scanone(BP025.sim,method = "imp",pheno.col = 6)
summary(out.imp.BP02586124)
plot(out.imp.BP02586124)
operm.imp.BP02586124 = scanone(BP025.sim,pheno.col = 6,method = "imp",n.perm = 1000)
summary(out.imp.BP02586124,perms = operm.imp.BP02586124,alpha = 0.05)

out.imp.BP025Pti2 = scanone(BP025.sim,method = "imp",pheno.col = 7)
summary(out.imp.BP025Pti2)
plot(out.imp.BP025Pti2)
operm.imp.BP025Pti2 = scanone(BP025.sim,pheno.col = 7,method = "imp",n.perm = 1000)
summary(out.imp.BP025Pti2,perms = operm.imp.BP025Pti2,alpha = 0.05)

out.imp.BP025ARB10 = scanone(BP025.sim,method = "imp",pheno.col = 8)
summary(out.imp.BP025ARB10)
plot(out.imp.BP025ARB10)
operm.imp.BP025ARB10 = scanone(BP025.sim,pheno.col = 8,method = "imp",n.perm = 1000)
summary(out.imp.BP025ARB10,perms = operm.imp.BP025ARB10,alpha = 0.05)

##two-QTL and multiple QTL
out2.imp.BP0253319 = scantwo(BP025.sim,pheno.col = 3,method = "imp")
plot(out2.imp.BP0253319)
operm2.imp.BP0253319 = scantwo(BP025.sim,pheno.col = 3,method = "imp",n.perm = 100)
str(operm2.imp.BP0253319)
summary(operm2.imp.BP0253319,alpha=0.05)
summary(out2.imp.BP0253319,perms = operm2.imp.BP0253319,what = "int")

out2.imp.BP025DW5 = scantwo(BP025.sim,pheno.col = 4,method = "imp")
plot(out2.imp.BP025DW5)
operm2.imp.BP025DW5 = scantwo(BP025.sim,pheno.col = 4,method = "imp",n.perm = 100)
str(operm2.imp.BP025DW5)
summary(operm2.imp.BP025DW5,alpha=0.05)
summary(out2.imp.BP025DW5,perms = operm2.imp.BP025DW5,what = "int")

out2.imp.BP02586124 = scantwo(BP025.sim,pheno.col = 6,method = "imp")
plot(out2.imp.BP02586124)
operm2.imp.BP02586124 = scantwo(BP025.sim,pheno.col = 6,method = "imp",n.perm = 100)
str(operm2.imp.BP02586124)
summary(operm2.imp.BP02586124,alpha=0.05)

out2.imp.BP025Pti2 = scantwo(BP025.sim,pheno.col = 7,method = "imp")
plot(out2.imp.BP025Pti2)
operm2.imp.BP025Pti2 = scantwo(BP025.sim,pheno.col = 7,method = "imp",n.perm = 100)
str(operm2.imp.BP025Pti2)
summary(operm2.imp.BP025Pti2,alpha=0.05)

out2.imp.BP025ARB10 = scantwo(BP025.sim,pheno.col = 8,method = "imp")
plot(out2.imp.BP025ARB10)
operm2.imp.BP025ARB10 = scantwo(BP025.sim,pheno.col = 8,method = "imp",n.perm = 100)
str(operm2.imp.BP025ARB10)
summary(operm2.imp.BP025ARB10,alpha=0.05)


#stepwise selection of multiple-QTL model
##calculate penalty threshold based on the results from two-dimentional two-QTL genome scans

####331-9
chr.3319 = c("5A")
pos.3319 = c(329)
qtl.3319 = makeqtl(BP025.sim,chr.3319,pos.3319)
my.formula = y~Q1
stepout.BP025.3319 = stepwiseqtl(BP025.sim,pheno.col = 3,additive.only = F,
                                 penalties = calc.penalties(operm2.imp.BP0253319, alpha = 0.05),
                                 qtl = qtl.3319,formula = my.formula,max.qtl=6,
                                 verbose=T, method="imp", keeptrace = TRUE)
rqtl.BP025.3319 = refineqtl(BP025.sim,pheno.col = 3,qtl = qtl.3319,method = "imp",formula = y~Q1)
rqtl.BP025.3319
plot(out.imp.BP0253319,chr = c("5A"),col = 4)
plotLodProfile(rqtl.BP025.3319,col = 2,add = T,qtl.labels = FALSE,lwd = 1)
lodint(rqtl.BP025.3319,qtl.index = 1,drop = 1.5)
fitqtl.3319 = fitqtl(BP025.sim,pheno.col = 3,qtl = rqtl.BP025.3319,method = "imp",
                     model=c("normal"), dropone=TRUE, get.ests=FALSE,
                     run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

####DW5
chr.DW5 = c("5A","5B")
pos.DW5 = c(322.72,5.76)
qtl.DW5 = makeqtl(BP025.sim,chr = chr.DW5,pos = pos.DW5)
my.formula = y~Q1+Q2
stepout.BP025.DW5 = stepwiseqtl(BP025.sim,pheno.col = 4,additive.only = F,
                                 penalties = calc.penalties(operm2.imp.BP025DW5, alpha = 0.05),
                                 qtl = qtl.DW5,formula = my.formula,max.qtl=6,
                                 verbose=T, method="imp", keeptrace = TRUE)
chr.DW5 = c("2B","4A","5A","5B")
pos.DW5 = c(398.52,0.00,322.72,2.00)
qtl.DW5 = makeqtl(BP025.sim,chr = chr.DW5,pos = pos.DW5)
my.formula = y ~ Q1 + Q2 + Q3 + Q4
rqtl.BO025.DW5 = refineqtl(BP025.sim,pheno.col = 4,qtl = qtl.DW5,method = "imp",formula = my.formula)
rqtl.BO025.DW5
plotLodProfile(rqtl.BO025.DW5,col=2,add = F,qtl.labels = F,lwd=1)
lodint(rqtl.BO025.DW5,qtl.index = 1,drop = 1.5)
lodint(rqtl.BO025.DW5,qtl.index = 2,drop = 1.5)
lodint(rqtl.BO025.DW5,qtl.index = 3,drop = 1.5)
lodint(rqtl.BO025.DW5,qtl.index = 4,drop = 1.5)
fitqtl.DW5 = fitqtl(BP025.sim,pheno.col = 4,qtl = rqtl.BO025.DW5,method = "imp",
                     model=c("normal"), dropone=TRUE, get.ests=FALSE,
                     run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
fitqtl.DW5$result.drop

####86124
chr.86124 = c("5A")
pos.86124 = c(320)
qtl.86124 = makeqtl(BP025.sim,chr = chr.86124,pos = pos.86124)
my.formula = y~Q1
stepout.BP025.86124 = stepwiseqtl(BP025.sim,pheno.col = 6,additive.only = F,
                                  penalties = calc.penalties(operm2.imp.BP02586124, alpha = 0.05),
                                  qtl = qtl.86124,formula = my.formula,max.qtl=6,
                                  verbose=T, method="imp", keeptrace = TRUE)

stepout.BP025.86124
chr.86124 = c("5A","5B","5B")
pos.86124 = c(319.56,139.67,240.0)
my.formula = y ~ Q1 + Q2 + Q3 
qtl.86124 = makeqtl(BP025.sim,chr = chr.86124,pos = pos.86124)
rqtl.BP025.86124 = refineqtl(BP025.sim,pheno.col = 6,qtl = qtl.86124,method = "imp",formula = my.formula)
plot(out.imp.BP02586124,chr = chr.86124,ylim = c(0,10))
plotLodProfile(rqtl.BP025.86124,col=2,add = T,qtl.labels = F,lwd=1)
lodint(rqtl.BP025.86124,qtl.index = 1,drop = 1.5)
lodint(rqtl.BP025.86124,qtl.index = 2,drop = 1.5)
lodint(rqtl.BP025.86124,qtl.index = 3,drop = 1.5)

fitqtl.86124 = fitqtl(BP025.sim,pheno.col = 6,qtl = rqtl.BP025.86124,method = "imp",
                    model=c("normal"), dropone=TRUE, get.ests=FALSE,
                    run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
fitqtl.86124$result.drop

#####Pti2
chr.Pti2 = c("5A")
pos.Pti2 = c(322.72)
qtl.Pti2 = makeqtl(BP025.sim, chr = chr.Pti2, pos = pos.Pti2)
my.formula = y~Q1
stepout.BP025.Pti2 = stepwiseqtl(BP025.sim,pheno.col = 7,additive.only = F,
                                 penalties = calc.penalties(operm2.imp.BP025Pti2, alpha = 0.05),
                                 qtl = qtl.Pti2,formula = my.formula,max.qtl=6,
                                 verbose=T, method="imp", keeptrace = TRUE)
stepout.BP025.Pti2
chr.Pti2 = c("1A","5A","5B","6B")
pos.Pti2 = c(75.924,32,140.187,273.137)
qtl.Pti2 = makeqtl(BP025.sim, chr = chr.Pti2, pos = pos.Pti2)
my.formula = y~Q1+Q2+Q3+Q4
rqtl.BP025.Pti2 = refineqtl(BP025.sim,pheno.col = 7,qtl = qtl.Pti2,method = "imp",formula = my.formula)
plot(out.imp.BP025Pti2,chr = chr.Pti2,ylim = c(0,10))
plotLodProfile(rqtl.BP025.Pti2,col=2,add = T,qtl.labels = F,lwd=1)
lodint(rqtl.BP025.Pti2,qtl.index = 1,drop = 1.5)
lodint(rqtl.BP025.Pti2,qtl.index = 2,drop = 1.5)
lodint(rqtl.BP025.Pti2,qtl.index = 3,drop = 1.5)
lodint(rqtl.BP025.Pti2,qtl.index = 4,drop = 1.5)
fitqtl.Pti2 = fitqtl(BP025.sim,pheno.col = 7,qtl = rqtl.BP025.Pti2,method = "imp",
       model=c("normal"), dropone=TRUE, get.ests=FALSE,
       run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
fitqtl.Pti2$result.drop
#bayesint(rqtl.BP025.Pti2,qtl.index = 1,prob = 0.99)

#####ARB10
chr.ARB10 = c("1A")
pos.ARB10 = c(78.0)
qtl.ARB10 = makeqtl(BP025.sim, chr = chr.ARB10, pos = pos.ARB10)
my.formula = y~Q1
stepout.BP025.ARB10 = stepwiseqtl(BP025.sim,pheno.col = 8,additive.only = F,
                                  penalties = calc.penalties(operm2.imp.BP025ARB10, alpha = 0.05),
                                  qtl = qtl.ARB10,formula = my.formula,max.qtl=6,
                                  verbose=T, method="imp", keeptrace = TRUE)

stepout.BP025.ARB10
chr.ARB10 = c("1A","4A","5A")
pos.ARB10 = c(75.924,10,319)
qtl.ARB10 = makeqtl(BP025.sim, chr = chr.ARB10, pos = pos.ARB10)
my.formula = y~Q1+Q2+Q3
rqtl.BP025.ARB10 = refineqtl(BP025.sim,pheno.col = 8,qtl = qtl.ARB10,method = "imp",formula = my.formula)
plot(out.imp.BP025ARB10,chr = chr.ARB10,ylim = c(0,10))
plotLodProfile(rqtl.BP025.ARB10,col=2,add = T,qtl.labels = F,lwd=1)
lodint(rqtl.BP025.ARB10,qtl.index = 1,drop = 1.5)
lodint(rqtl.BP025.ARB10,qtl.index = 2,drop = 1.5)
lodint(rqtl.BP025.ARB10,qtl.index = 3,drop = 1.5)
fitqtl.ARB10 = fitqtl(BP025.sim,pheno.col = 8, qtl = rqtl.BP025.ARB10,method = "imp",
                      model=c("normal"), dropone=TRUE, get.ests=FALSE,
                      run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
fitqtl.ARB10$result.drop
