####RHM
#Running environment setting
getwd()
setwd("/Users/xuehui.li/Desktop/XL/research/crops/common_wheat/spring wheat/project2_fhb in NDSU_HRSW/GWAS/RHM")
list.files()


setwd("/mnt/volum4/RHM_Spring_wheat/xuehui_editted")
######Begin
library(EMMREML)
library(rrBLUP)

####Spring Wheat genotype and phenotype
marker_data1 = read.table("FHB-nursery spring wheat lines_439 GBS SNPna10ABD_LDKNNimp.F01.txt", header = T, stringsAsFactors = F)
marker_data2 = read.table("FHB-nursery spring wheat lines_439 GBS SNPna50ABD_LDKNNimp.F01.txt", header = T, stringsAsFactors = F)
pheno_data = read.table("phenotype_estimated blups of HSEVR across all trials_439.txt", header = T, na.strings = "NA")

dim(marker_data1)
dim(marker_data2)
dim(pheno_data)
marker_data1[1:5,1:5]
marker_data2[1:5,1:5]
pheno_data[1:5,1:5]
rownames(pheno_data) = pheno_data[,1]
pheno_data = pheno_data[,-1]

#################################################
#Missing value must be filled with NA, not na.
###Marker Imputation
imputed1 = A.mat(marker_data1[,-1], return.imputed = T)
imputed_marker1 = imputed1$imputed
imputed_marker1[1:5,1:5]

imputed2 = A.mat(marker_data2[,-1], return.imputed = T)
imputed_marker2 = imputed2$imputed
imputed_marker2[1:5,1:5]

find_major_allele_frequncy =  function(marker){
  allele_0 = length(which(marker == 0))
  allele_1 = length(which(marker == 1))
  return(max(allele_0,allele_1)/length(marker))
}

#######Phenotype HSEVR.a2
HSEVR.a2 = pheno_data$HSEVR.a2
##HSEVR.a2 = pheno_data$HSEVR.a2[-which(is.na(pheno_data$HSEVR.a2))]

##imputed_marker1 = imputed_marker1[-which(is.na(pheno_data$HSEVR.a2)),]
##imputed_marker2 = imputed_marker2[-which(is.na(pheno_data$HSEVR.a2)),]

dim(imputed_marker1)
dim(imputed_marker2)
length(HSEVR.a2)

MAF_trait =  apply(imputed_marker1, 2,find_major_allele_frequncy)    

background_segement = imputed_marker1

kinship = (as.matrix(t(background_segement)) %*%  as.matrix(background_segement)/2/sum(MAF_trait*(1-MAF_trait)))

###########Single thread execution
#likelihood_container = c() ### Store likelihood generated from full model
#
# n = 1d
# 
# #load EMMREML package
# ptm = proc.time()
# 
# while(n < 4){
#   
#   start_marker = n 
#   
#   end_marker = n + 50
#   
#   test_segment = marker_HSEVR.a2[,(start_marker+1):(end_marker+1)]
#   
#   segment_GR = ((as.matrix(t(test_segment)) %*% as.matrix(test_segment)))/2/sum(MAF_trait[start_marker:end_marker]*(1-MAF_trait[start_marker:end_marker]))
#   
#   full_model = emmremlMultiKernel(y = HSEVR.a2, X = matrix(rep(1, length(HSEVR.a2)),ncol = 1), Klist = list(segment_GR,kinship), Zlist = list(test_segment, background_segement))
#   
#   likelihood_container = c(likelihood_container, full_model$loglik)
#   
#   print(n) ###Make progress visibile
#   
#   n = n + 1
# }
# 
# proc.time() - ptm

#X must be a matrix. It can be a matrix with only one column.
#Missing must be avoid in Y

#######################################
null_model = emmremlMultiKernel(y = HSEVR.a2, X = matrix(rep(mean(HSEVR.a2), length(HSEVR.a2)),ncol = 1), Klist = list(kinship), Zlist = list(background_segement))
null_model$loglik

#####Parallel computing
library(parallel)
library(foreach)
library(doParallel)

MAF_trait_NA50 =  apply(imputed_marker2, 2,find_major_allele_frequncy)  

emmreml_seq = function(n){
  
  start_marker = n 
  
  end_marker = n+200 ###define window size
  
  test_segment = imputed_marker2[,(start_marker):(end_marker)]
  
  segment_GR = ((as.matrix(t(test_segment)) %*% as.matrix(test_segment)))/2/sum(MAF_trait_NA50[start_marker:end_marker]*(1-MAF_trait_NA50[start_marker:end_marker]))
  
  full_model = emmremlMultiKernel(y = HSEVR.a2, X = matrix(rep(1, length(HSEVR.a2)),ncol = 1), Klist = list(segment_GR,kinship), Zlist = list(test_segment, background_segement))
  
  result = c(colnames(imputed_marker2)[start_marker], colnames(imputed_marker2)[end_marker], full_model$loglik)
  
  return(result)
}

###############################
num_cores = detectCores()-2 #Lab desktop = 12; abc cloud = 12; office desktop = 6
registerDoParallel(num_cores)
ptm = proc.time()

###the argument, by define window moving size
# 1	S1A
# 5295	S1B
# 12216	S1D
# 15013	S2A
# 19891	S2B
# 28510	S2D
# 31822	S3A
# 38243	S3B
# 47062	S3D
# 49779	S4A
# 54828	S4B
# 59585	S4D
# 60753	S5A
# 66467	S5B
# 73840	S5D
# 76266	S6A
# 81064	S6B
# 89618	S6D
# 92148	S7A
# 101392	S7B
# 108731	S7D

#Task assignment: Lab desktop chr 4,5,6,7; abc cloud chr 2,3; office desktop chr 1

test_contain_chr1A = foreach(m = seq(from = 1, to = 5294, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr1A)
write.table(test_contain_chr1A, file="test_contain_chr1A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr1B = foreach(m = seq(from = 5295, to = 12215, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr1B)
write.table(test_contain_chr1B, file="test_contain_chr1B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr1D = foreach(m = seq(from = 12216, to = 15012, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr1D)
write.table(test_contain_chr1D, file="test_contain_chr1D.txt", quote = F, col.names = F,row.names = F)


test_contain_chr2A = foreach(m = seq(from = 15013, to = 19890, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr2A)
write.table(test_contain_chr2A, file="test_contain_chr2A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr2B = foreach(m = seq(from = 19891, to = 28509, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr2B)
write.table(test_contain_chr2B, file="test_contain_chr2B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr2D = foreach(m = seq(from = 28510, to = 31821, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr2D)
write.table(test_contain_chr2D, file="test_contain_chr2D.txt", quote = F, col.names = F,row.names = F)

test_contain_chr3A = foreach(m = seq(from = 31822, to = 38242, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr3A)
write.table(test_contain_chr3A, file="test_contain_chr3A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr3B = foreach(m = seq(from = 38243, to = 47061, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr3B)
write.table(test_contain_chr3B, file="test_contain_chr3B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr3D = foreach(m = seq(from = 47062, to = 49778, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr3D)
write.table(test_contain_chr3D, file="test_contain_chr3D.txt", quote = F, col.names = F,row.names = F)


test_contain_chr4A = foreach(m = seq(from = 49779, to = 54827, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr4A)
write.table(test_contain_chr4A, file="test_contain_chr4A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr4B = foreach(m = seq(from = 54828, to = 59584, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr4B)
write.table(test_contain_chr4B, file="test_contain_chr4B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr4D = foreach(m = seq(from = 59585, to = 60752, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr4D)
write.table(test_contain_chr4D, file="test_contain_chr4D.txt", quote = F, col.names = F,row.names = F)

test_contain_chr5A = foreach(m = seq(from = 60753, to = 66466, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr5A)
write.table(test_contain_chr5A, file="test_contain_chr5A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr5B = foreach(m = seq(from = 66467, to = 73479, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr5B)
write.table(test_contain_chr5B, file="test_contain_chr5B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr5D = foreach(m = seq(from = 73840, to = 76265, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr5D)
write.table(test_contain_chr5D, file="test_contain_chr5D.txt", quote = F, col.names = F,row.names = F)

test_contain_chr6A = foreach(m = seq(from = 76266, to = 81063, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr6A)
write.table(test_contain_chr6A, file="test_contain_chr6A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr6B = foreach(m = seq(from = 81064, to = 89617, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr6B)
write.table(test_contain_chr6B, file="test_contain_chr6B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr6D = foreach(m = seq(from = 89618, to = 92147, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr6D)
write.table(test_contain_chr6D, file="test_contain_chr6D.txt", quote = F, col.names = F,row.names = F)

test_contain_chr7A = foreach(m = seq(from = 92148, to = 101391, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr7A)
write.table(test_contain_chr7A, file="test_contain_chr7A.txt", quote = F, col.names = F,row.names = F)

test_contain_chr7B = foreach(m = seq(from = 101392, to = 108730, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr7B)
write.table(test_contain_chr7B, file="test_contain_chr7B.txt", quote = F, col.names = F,row.names = F)

test_contain_chr7D = foreach(m = seq(from = 108731, to = 112250, by = 100), .combine =rbind, .errorhandling = "pass", .packages = "EMMREML", .inorder = T) %dopar%
  emmreml_seq(m)
length(test_contain_chr7D)
write.table(test_contain_chr7D, file="test_contain_chr7D.txt", quote = F, col.names = F,row.names = F)

proc.time() - ptm ###Record multithreading running time
stopImplicitCluster()

########################
test_contain_chr1A # this contains all loglikelihood results
null_model$loglik

D_chr1A=2*(test_contain_chr1A-null_model$loglik)

pvalue.chr1A = pchisq(2*(test_contain_chr1A-null_model$loglik),df=1, lower.tail = FALSE)
pvalue.chr1A = pchisq(D_chr1A,df=1, lower.tail = FALSE)
pvalue.chr1A

qvalue = p.adjust(pvalue, method = "fdr")
qvalue
plot(qvalue)



