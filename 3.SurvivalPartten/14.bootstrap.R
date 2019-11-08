#-------------------------------------
#' bootstrap cox回归模型的自证实
#-------------------------------------

##导入cox.os.all.gene OS生存数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.OS.all.gene.RData")

mvcoxres_1 <- function(mod){
  hr <- round(summary(mod)$conf.int[, 1], 2)
  lower <- round(summary(mod)$conf.int[, 3], 2)
  upper <- round(summary(mod)$conf.int[, 4], 2)
  p <- round(summary(mod)$coefficients[, 5], 3)
  res <- data.frame(hr,lower,upper, p, stringsAsFactors = F)
  res$p[res$p < 0.001] <- "<.001"
  colnames(res) <- c("HR","lower","upper","p-value")
  return(res)
}

library(bootStepAIC)
library(survival)
library(help="bootStepAIC")
?boot.stepAIC

## 所有临床字段的结果
coxFit <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.OS.all.gene)[-c(1:3)],collapse = "+"))), data=na.omit(COX.OS.all.gene)) #boot.stepAIC不允许存在NA
a1 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a2 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
a3 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a4 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
a5 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a6 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a7 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
save(a1,a2,a3,a4,a5,a6,a7,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/bootstrap1.RData")

## 不加Residual临床字段的结果
coxFit1 <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.OS.all.gene)[-c(1:3,10)],collapse = "+"))),data=na.omit(COX.OS.all.gene))
aaaa <- boot.stepAIC(coxFit1, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
aaaaa <- boot.stepAIC(coxFit1, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both") 
a6 <- boot.stepAIC(coxFit1, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a9 <- boot.stepAIC(coxFit1, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
save(aa,aaa,a7,a8,aaaa,aaaaa,a6,a9,a10,a11,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/bootstrap.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/bootstrap.RData")

library(bootStepAIC)
library(survival)
bootstrap.resultTable.a7 <- lapply(a7$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a11 <- lapply(a11$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.aa <- lapply(aa$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a10 <- lapply(a10$BootStepAIC,mvcoxres_1)
class(bootstrap.resultTable[[1]]) # data.frame
sj <- function(cox.resuts,rownam){
  ##cox.results
  ##rownam 
  result <- matrix(0,nrow=length(rownam),ncol=3,dimnames = list(rownam,c("mean HR","mean lower","mean upper")))
  for(j in 1:length(rownam)){
    temp <- list()
    for(i in 1:length(cox.resuts)){
      pos <- rownames(cox.resuts[[i]]) %in% rownam[j]
      if(any(pos)){
        temp[[i]] <- cox.resuts[[i]][pos,]
      }
    }
    tempMatrix <- do.call(rbind,temp)
    tempMatrix <- as.matrix(tempMatrix[,1:3])
    tempMatrix[tempMatrix >= 50] <- Inf ## 将太大的值设为Inf
    tempMatrix[tempMatrix <= -50] <- -Inf ## 将太大的值设为Inf
    tempMatrix[!is.finite(tempMatrix)] <- NA  ## 将Inf和-Inf改为NA
    tempRow <- colMeans(na.omit(tempMatrix))
    result[j,] <- tempRow
  }
  return(result)
}
rownam <- rownames(mvcoxres_1(coxFit))
boot.Mean.a7 <- sj(bootstrap.resultTable.a7,rownam) #a7
boot.Mean.a11 <- sj(bootstrap.resultTable.a11,rownam) #a11
boot.Mean.aa <- sj(bootstrap.resultTable.aa,rownam) #aa
boot.Mean.a10 <- sj(bootstrap.resultTable.a10,rownam) #a10
save(boot.Mean.a7,boot.Mean.a11,boot.Mean.aa,boot.Mean.a10,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/bootstrap.mean.HR.RData")

a1
# Covariates selected
# (%)
# age         100.0
# M_stage      97.4
# CASP8        92.3
# TNM_stage    90.4
# ARID1A       86.1
# SMAD2        86.0
# ANK1         61.7
# Tumour_site  48.5
# N_stage      43.9
# T_stage      41.2
# GRIN2A       40.7
# sex          39.0
# MSS_state    29.5
a2
# Covariates selected
# (%)
# age         100.0
# M_stage      98.5
# CASP8        92.6
# TNM_stage    89.4
# ARID1A       86.7
# SMAD2        85.8
# ANK1         61.5
# Tumour_site  50.1
# N_stage      43.5
# GRIN2A       41.9
# sex          40.0
# T_stage      39.9
# MSS_state    35.1
a3
# Covariates selected
# (%)
# age         99.9
# M_stage     98.6
# CASP8       91.9
# TNM_stage   87.4
# SMAD2       87.1
# ARID1A      85.9
# ANK1        63.0
# Tumour_site 46.5
# N_stage     44.9
# T_stage     42.5
# GRIN2A      39.2
# sex         36.8
# MSS_state   31.7
a4
# Covariates selected
# (%)
# age         100.0
# M_stage      98.2
# CASP8        91.0
# TNM_stage    88.7
# ARID1A       88.0
# SMAD2        85.8
# ANK1         64.9
# Tumour_site  48.5
# N_stage      45.7
# GRIN2A       41.1
# T_stage      39.3
# sex          36.1
# MSS_state    34.1
a5
# Covariates selected
# (%)
# age         100.0
# M_stage      98.2
# CASP8        91.8
# TNM_stage    89.6
# ARID1A       87.1
# SMAD2        83.3
# ANK1         65.0
# Tumour_site  48.3
# N_stage      46.3
# GRIN2A       41.9
# T_stage      39.5
# sex          39.0
# MSS_state    31.6
a6
# Covariates selected
# (%)
# age         100.0
# M_stage      97.4
# CASP8        92.6
# TNM_stage    88.4
# SMAD2        86.1
# ARID1A       85.4
# ANK1         60.1
# N_stage      46.8
# Tumour_site  46.4
# sex          40.1
# T_stage      39.0
# GRIN2A       38.8
# MSS_state    34.1
a7
# Covariates selected
# (%)
# age         100.0
# M_stage      98.3
# CASP8        92.3
# TNM_stage    89.2
# ARID1A       86.4
# SMAD2        86.3
# ANK1         60.5
# Tumour_site  51.3
# N_stage      43.6
# GRIN2A       39.9
# sex          39.5
# T_stage      39.4
# MSS_state    30.1
bootstrap.resultTable.a1 <- lapply(a1$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a2 <- lapply(a2$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a3 <- lapply(a3$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a4 <- lapply(a4$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a5 <- lapply(a5$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a6 <- lapply(a6$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a7 <- lapply(a7$BootStepAIC,mvcoxres_1)
sj <- function(cox.resuts,rownam){
  ##cox.results
  ##rownam 
  result <- matrix(0,nrow=length(rownam),ncol=3,dimnames = list(rownam,c("mean HR","mean lower","mean upper")))
  for(j in 1:length(rownam)){
    temp <- list()
    for(i in 1:length(cox.resuts)){
      pos <- rownames(cox.resuts[[i]]) %in% rownam[j]
      if(any(pos)){
        temp[[i]] <- cox.resuts[[i]][pos,]
      }
    }
    tempMatrix <- do.call(rbind,temp)
    tempMatrix <- as.matrix(tempMatrix[,1:3])
    tempMatrix[tempMatrix >= 50] <- Inf ## 将太大的值设为Inf
    tempMatrix[tempMatrix <= -50] <- -Inf ## 将太大的值设为Inf
    tempMatrix[!is.finite(tempMatrix)] <- NA  ## 将Inf和-Inf改为NA
    tempRow <- colMeans(na.omit(tempMatrix))
    result[j,] <- tempRow
  }
  return(result)
}
cox.resuts <- bootstrap.resultTable.a4
#result <- matrix(0,nrow=length(rownam),ncol=3,dimnames = list(rownam,c("mean HR","mean lower","mean upper")))
"N_stage"
"T_stage"
"ANK1ANK1 Subclonal"
"SMAD2SMAD2 Subclonal"
  temp <- list()
  for(i in 1:length(cox.resuts)){
    pos <- rownames(cox.resuts[[i]]) %in% "CASP8CASP8 Subclonal"
    if(any(pos)){
      temp[[i]] <- cox.resuts[[i]][pos,]
    }
  }
  tempMatrix <- do.call(rbind,temp)
  tempMatrix <- as.matrix(tempMatrix[,1:3])
  tempMatrix[tempMatrix >= 50] <- Inf ## 将太大的值设为Inf
  tempMatrix[tempMatrix <= -50] <- -Inf ## 将太大的值设为Inf
  tempMatrix[!is.finite(tempMatrix)] <- NA  ## 将Inf和-Inf改为NA
  tempRow <- colMeans(na.omit(tempMatrix))
  #result[j,] <- tempRow


boot.Mean.a1 <- sj(bootstrap.resultTable.a1,rownam) #a1
boot.Mean.a2 <- sj(bootstrap.resultTable.a2,rownam) #a2
boot.Mean.a3 <- sj(bootstrap.resultTable.a3,rownam) #a3
boot.Mean.a4 <- sj(bootstrap.resultTable.a4,rownam) #a4
boot.Mean.a5 <- sj(bootstrap.resultTable.a5,rownam) #a5
boot.Mean.a6 <- sj(bootstrap.resultTable.a6,rownam) #a6
boot.Mean.a7 <- sj(bootstrap.resultTable.a7,rownam) #a7

bootstrap.resultTable.a1[[1]]
bootstrap.resultTable.a1[[2]]






######################DFS bootstrap#########################
library(bootStepAIC)
library(survival)
coxFit <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.DFS.all.gene)[-c(1:3,10)],collapse = "+"))),data=na.omit(COX.DFS.all.gene))
a<- boot.stepAIC(coxFit, data = na.omit(COX.DFS.all.gene), B = 1000,direction = "backward")
aa <- boot.stepAIC(coxFit, data = na.omit(COX.DFS.all.gene), B = 1000,direction = "both") 
aaa <- boot.stepAIC(coxFit, data = na.omit(COX.DFS.all.gene), B = 1000,direction = "both")
aaaa <- boot.stepAIC(coxFit, data = na.omit(COX.DFS.all.gene), B = 1000,direction = "backward")
save(a,aa,aaa,aaaa,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.bootstrap.RData")
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.bootstrap.RData")
library(bootStepAIC)
bootstrap.resultTable.a <- lapply(a$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.aa <- lapply(aa$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.aaa <- lapply(aaa$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.aaaa <- lapply(aaaa$BootStepAIC,mvcoxres_1)
sj <- function(cox.resuts,rownam){
  ##cox.results
  ##rownam 
  result <- matrix(0,nrow=length(rownam),ncol=3,dimnames = list(rownam,c("mean HR","mean lower","mean upper")))
  for(j in 1:length(rownam)){
    temp <- list()
    for(i in 1:length(cox.resuts)){
      pos <- rownames(cox.resuts[[i]]) %in% rownam[j]
      if(any(pos)){
        temp[[i]] <- cox.resuts[[i]][pos,]
      }
    }
    tempMatrix <- do.call(rbind,temp)
    tempMatrix <- as.matrix(tempMatrix[,1:3])
    tempMatrix[tempMatrix>=7.261615e+24] <- Inf ## 将太大的值设为Inf
    tempMatrix[!is.finite(tempMatrix)] <- NA  ## 将Inf和-Inf改为NA
    tempRow <- colMeans(na.omit(tempMatrix))
    result[j,] <- tempRow
  }
  return(result)
}
rownam <- rownames(mvcoxres_1(coxFit))
boot.Mean.a <- sj(bootstrap.resultTable.a,rownam) #a
boot.Mean.aa <- sj(bootstrap.resultTable.aa,rownam) #aa
boot.Mean.aaa <- sj(bootstrap.resultTable.aaa,rownam) #aaa
boot.Mean.aaaa <- sj(bootstrap.resultTable.aaaa,rownam) #aaaa
save(boot.Mean.a,boot.Mean.aa,boot.Mean.aaa,boot.Mean.aaaa,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.bootstrap.mean.HR.RData")


