setwd("Y:/Yuan/Breeding_Tanspot_Regional_Heritability_mapping")
list.files()
library("EMMREML")
library("rrBLUP")

markerData = read.table("Numeric_format.txt",header = T,na.strings = "NA",row.names = 1)

allPheno = read.table("phenotype_AYT12-16_ND12-quality_BLUPs.txt",header = T,na.strings = "NA",stringsAsFactors = F,row.names = 1)
sedVolume=allPheno[,7]

markerData=markerData[row.names(allPheno),]
validPheno=!is.na(sedVolume)
sedVolume=sedVolume[validPheno]
sedVolume=as.numeric(sedVolume)
mean(sedVolume)

markerData=markerData[validPheno,]



#################################################
#Missing value must be filled with NA, not na.
###Marker Imputation
imputed1 = A.mat(markerData, return.imputed = T)
imputed_marker1 = imputed1$imputed
imputed_marker1[1:5,1:5]

imputed2 = A.mat(markerData, return.imputed = T)
imputed_marker2 = imputed2$imputed
imputed_marker2[1:5,1:5]



find_major_allele_frequncy =  function(marker){
  allele_0 = length(which(marker == 0))
  allele_1 = length(which(marker == 1))
  return(max(allele_0,allele_1)/length(marker))
}

#######Phenotype HSEVR.a2
HSEVR.a2 = sedVolume
##HSEVR.a2 = pheno_data$HSEVR.a2[-which(is.na(pheno_data$HSEVR.a2))]

##imputed_marker1 = imputed_marker1[-which(is.na(pheno_data$HSEVR.a2)),]
##imputed_marker2 = imputed_marker2[-which(is.na(pheno_data$HSEVR.a2)),]

dim(imputed_marker1)
dim(imputed_marker2)
length(HSEVR.a2)

MAF_trait =  apply(imputed_marker1, 2,find_major_allele_frequncy)    

background_segement = imputed_marker1

kinship = (as.matrix(t(background_segement)) %*%  as.matrix(background_segement)/2/sum(MAF_trait*(1-MAF_trait)))

#X must be a matrix. It can be a matrix with only one column.
#Missing must be avoid in Y

#######################################
null_model = emmremlMultiKernel(y = HSEVR.a2, X = matrix(rep(mean(HSEVR.a2), length(HSEVR.a2)),ncol = 1), Klist = list(kinship), Zlist = list(background_segement))
null_model$loglik

#####Parallel computing
library(foreach)
library(doParallel)

MAF_trait_NA50 =  apply(imputed_marker2, 2,find_major_allele_frequncy)  

emmreml_seq = function(n){
  
  start_marker = n 
  
  end_marker = n+50 ###define window size
  
  test_segment = imputed_marker2[,(start_marker):(end_marker)]
  
  segment_GR = ((as.matrix(t(test_segment)) %*% as.matrix(test_segment)))/2/sum(MAF_trait_NA50[start_marker:end_marker]*(1-MAF_trait_NA50[start_marker:end_marker]))
  
  full_model = emmremlMultiKernel(y = HSEVR.a2, X = matrix(rep(mean(HSEVR.a2), length(HSEVR.a2)),ncol = 1), Klist = list(segment_GR,kinship), Zlist = list(test_segment, background_segement))
  
  result = c(colnames(imputed_marker2)[start_marker], colnames(imputed_marker2)[end_marker], full_model$loglik)
  
  return(result)
}



x=emmreml_seq(1)
x2=emmreml_seq(2000)


###############################
num_cores = 12 #Lab desktop = 12; abc cloud = 12; office desktop = 6
registerDoParallel(num_cores)
ptm = proc.time() 



test_contain_chr1A = foreach(m = seq(from = 1, to = 757, by = 25), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr1A)
write.table(test_contain_chr1A, file="test_contain_chr1A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr1B = foreach(m = seq(from = 758, to = 2084, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr1B)
write.table(test_contain_chr1B, file="test_contain_chr1B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr2A = foreach(m = seq(from = 2085, to = 2568, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr2A)
write.table(test_contain_chr2A, file="test_contain_chr2A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr2B = foreach(m = seq(from = 2569, to = 3679, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr2B)
write.table(test_contain_chr2B, file="test_contain_chr2B.txt", quote = F, col.names = F,row.names = F)


test_contain_chr3A = foreach(m = seq(from = 3680, to = 4518, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr3A)
write.table(test_contain_chr3A, file="test_contain_chr3A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr3B = foreach(m = seq(from = 4519, to = 5441, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr3B)
write.table(test_contain_chr3B, file="test_contain_chr3B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr4A = foreach(m = seq(from = 5442, to = 6552, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr4A)
write.table(test_contain_chr4A, file="test_contain_chr4A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr4B = foreach(m = seq(from = 6553, to = 7105, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr4B)
write.table(test_contain_chr4B, file="test_contain_chr4B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr5A = foreach(m = seq(from = 7106, to = 7857, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr5A)
write.table(test_contain_chr5A, file="test_contain_chr5A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr5B = foreach(m = seq(from = 7858, to = 8894, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr5B)
write.table(test_contain_chr5B, file="test_contain_chr5B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr6A = foreach(m = seq(from = 8895, to = 10028, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr6A)
write.table(test_contain_chr6A, file="test_contain_chr6A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr6B = foreach(m = seq(from = 10029, to = 11202, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr6B)
write.table(test_contain_chr6B, file="test_contain_chr6B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr7A = foreach(m = seq(from = 11203, to = 12158, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr7A)
write.table(test_contain_chr7A, file="test_contain_chr7A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr7B = foreach(m = seq(from = 12159, to = 13161, by = 10), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr7B)
write.table(test_contain_chr7B, file="test_contain_chr7B.txt", quote = F, col.names = F,row.names = F)

proc.time() - ptm ###Record multithreading running time
stopImplicitCluster()

########################
as.numeric(test_contain_chr1A[,3]) # this contains all loglikelihood results
null_model$loglik

D_chr1A=2*(as.numeric(test_contain_chr1A[,3])-null_model$loglik)

pvalue.chr1A = pchisq(2*(as.numeric(test_contain_chr1A[,3])-min(as.numeric(test_contain_chr1A[,3]))),df=1, lower.tail = FALSE)
pvalue.chr1A = pchisq(D_chr1A,df=1, lower.tail = FALSE)
pvalue.chr1A

qvalue = p.adjust(pvalue.chr1A, method = "fdr")
qvalue
plot(-log10(qvalue))

