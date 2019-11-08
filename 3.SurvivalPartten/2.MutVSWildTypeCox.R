#######################################################################################
#driver基因突变和不突变情况下的OS和DFS 单因素和多因素COX分析
#######################################################################################
##################导入生存数据####################
#导入OS生存数据                                                                           
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.SingleMarker.RData");
#导入DFS生存数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.byDriverGene.SingleMarker.RData");

#处理清洗OS COX生存数据
{
  OS.significant <- c("AKT1","AXIN2","BCLAF1","CASP8","CNBD1","CD70","EPHA7","GRIN2A","SMAD2")
  pos <- match(OS.significant,names(OS.byDriverGene.SingleMarker))
  COX.OS.byDriverGene.SingleMarker <- vector(mode="list",length=length(OS.significant))
  for(i in 1:length(pos)){
    COX.OS.byDriverGene.SingleMarker[[i]] <- OS.byDriverGene.SingleMarker[[pos[i]]]
  }
  names(COX.OS.byDriverGene.SingleMarker) <- OS.significant
  head(COX.OS.byDriverGene.SingleMarker[[1]])
  COX.OS.byDriverGene.SingleMarker <- lapply(COX.OS.byDriverGene.SingleMarker,function(x){
    #temp <- x[,c(1,22,23,2,6,7,8,9,19,21,26)];
    #colnames(temp) <- c("Patient_ID","event","time","sex","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site","Mutation_status");
    temp <- x[,c(1,10,11,2,3,4,5,6,7,9,13,14)] #剔除亚型信息
    #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Subtype","Tumour_site","Mutation_status")
    #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Tumour_site","Mutation_status")
    #对协变量进行factor转换，控制levels的第一项为参考项,目前默认第一个level为参考项 
    temp$sex <- factor(temp$sex,levels=c("Female","Male"))
    temp$MSS_state[grep("MSI",temp$MSS_state)] <- "MSI"
    temp$MSS_state <- factor(temp$MSS_state,levels=c("MSI","MSS"))
    temp$T_stage[grep("T1|T2",temp$T_stage)] <- "T1/T2"
    temp$T_stage[grep("T3|T4",temp$T_stage)] <- "T3/T4"
    temp$T_stage <- factor(temp$T_stage,levels=c("T1/T2","T3/T4"))
    temp$N_stage[grep("N1|N2",temp$N_stage)] <- "N1/N2"
    temp$N_stage <- factor(temp$N_stage,levels=c("N0","N1/N2"))
    temp$M_stage <- factor(temp$M_stage,levels=c("M0","M1"))
    temp$TNM_stage[grep("Stage 1|Stage 2",temp$TNM_stage)] <- "Stage 1/Stage 2"
    temp$TNM_stage[grep("Stage 3|Stage 4",temp$TNM_stage)] <- "Stage 3/Stage 4"
    temp$TNM_stage <- factor(temp$TNM_stage,levels=c("Stage 1/Stage 2","Stage 3/Stage 4"))
    # temp$Residual[grep("R1|R2",temp$Residual)] <- "R1/R2"
    # temp$Residual <- factor(temp$Residual,levels=c("R0","R1/R2"))
    temp$age <- factor(temp$age,levels=c("<67.5",">=67.5"))
    #temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
    temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
    temp
  })
  names(COX.OS.byDriverGene.SingleMarker) <- OS.significant
  head(COX.OS.byDriverGene.SingleMarker[[1]])
}

## 进行OS 单多因素COX分析 得到结果COX.OS.result
COX.OS.result <- lapply(COX.OS.byDriverGene.SingleMarker,function(x){
	Cox.function(x);
})
names(COX.OS.result) <- OS.significant
save(COX.OS.result,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.OS.result.RData")

## 进行OS 多个基因组合多因素cox分析
{
  length(COX.OS.byDriverGene.SingleMarker)
  names(COX.OS.byDriverGene.SingleMarker)
  head(COX.OS.byDriverGene.SingleMarker[[1]])
  library(dplyr)
  temp <- full_join(COX.OS.byDriverGene.SingleMarker[[1]],COX.OS.byDriverGene.SingleMarker[[2]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site"))
  for(i in names(COX.OS.byDriverGene.SingleMarker)[c(-1,-2)]){
    temp <- full_join(temp,COX.OS.byDriverGene.SingleMarker[[i]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site"))
  }
  colnames(temp)[c(12,13,14,15,16,17,18,19,20)] <- OS.significant
  head(temp)
  COX.OS.byDriverGene.SingleMarker.allgenes <- temp
}
COX.OS.allgenes.result <- Cox.function(COX.OS.byDriverGene.SingleMarker.allgenes)
Cox.function(na.omit(COX.OS.byDriverGene.SingleMarker.allgenes))
save(COX.OS.allgenes.result,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.OS.allgenes.result.RData")

# 进行Test.proportional.hazards分析
{
  ## 只对基因突变克隆状态进行等比例风险分析
  uni.proportional.hazards <- list()
  for(i in names(COX.OS.byDriverGene.SingleMarker)){
    uni.proportional.hazards[[i]] <- Test.proportional.hazards(COX.OS.byDriverGene.SingleMarker[[i]][,c(1,2,3,12)],paste0("F:/",i,".pdf"))
  }
  ## 对所有临床字段进行等比例风险假设分析
  proportional.hazards <- list()
  for(i in names(COX.OS.byDriverGene.SingleMarker.allgenes)){
    proportional.hazards[[i]] <- Test.proportional.hazards(COX.OS.byDriverGene.SingleMarker.allgenes[[i]],paste0("F:/",i,".pdf"))
  }
  Test.proportional.hazards(COX.OS.byDriverGene.SingleMarker.allgenes,paste0("F:/",i,".pdf"))
}

## 进行OS多个基因组合逐步cox回归分析
{
  ##提取COX模型的HR,CI,p值
  mvcoxres <- function(mod){
    hr.ci <- paste0(
      round(summary(mod)$conf.int[, 1], 2), " (",
      round(summary(mod)$conf.int[, 3], 2), ", ",
      round(summary(mod)$conf.int[, 4], 2), ")"
    )
    p <- round(summary(mod)$coefficients[, 5], 3)
    res <- data.frame(hr.ci, p, stringsAsFactors = F)
    res$p[res$p < 0.001] <- "<.001"
    colnames(res) <- c("**HR (95% CI)**", "**p-value**")
    return(res)
  }
  library(MASS)
  library(survival)
  head(COX.OS.byDriverGene.SingleMarker.allgenes)
  coxmodel <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.OS.byDriverGene.SingleMarker.allgenes)[-c(1:3)],collapse = "+"))),
                    data=na.omit(COX.OS.byDriverGene.SingleMarker.allgenes)) #stepAIC 不允许生存数据有NA,考虑补缺失值
  summary(coxmodel)
  extractAIC(coxmodel)
  step.cox <- stepAIC(coxmodel,direction = "both",trace = F) #stepAIC 不允许生存数据有NA
  #step.cox.result <- summary(step.cox)
  stepwise.table <- mvcoxres(step.cox)
  save(stepwise.table,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.os.MutvsWT.stepwise.RData")
}
## 处理清洗DFS COX生存数据
DFS.significant <- c("AKT1","AXIN2","CUX1","KRAS")
pos <- match(DFS.significant,names(DFS.byDriverGene.SingleMarker))
COX.DFS.byDriverGene.SingleMarker <- vector(mode="list",length=length(DFS.significant))
for(i in 1:length(pos)){
   COX.DFS.byDriverGene.SingleMarker[[i]] <- DFS.byDriverGene.SingleMarker[[pos[i]]]
}
names(COX.DFS.byDriverGene.SingleMarker) <- DFS.significant

COX.DFS.byDriverGene.SingleMarker <- lapply(COX.DFS.byDriverGene.SingleMarker,function(x){
  #temp <- x[,c(1,22,23,2,6,7,8,9,19,21,26)];
  #colnames(temp) <- c("Patient_ID","event","time","sex","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site","Mutation_status");
  temp <- x[,c(1,12,13,2,3,4,5,6,7,9,17,18)] #剔除亚型信息
  colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site","Mutation_status")
  #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Tumour_site","Mutation_status")
  #对协变量进行factor转换，控制levels的第一项为参考项,目前默认第一个level为参考项 
  temp$sex <- factor(temp$sex,levels=c("Female","Male"))
  temp$MSS_state[grep("MSI",temp$MSS_state)] <- "MSI"
  temp$MSS_state <- factor(temp$MSS_state,levels=c("MSI","MSS"))
  temp$T_stage[grep("T1|T2",temp$T_stage)] <- "T1/T2"
  temp$T_stage[grep("T3|T4",temp$T_stage)] <- "T3/T4"
  temp$T_stage <- factor(temp$T_stage,levels=c("T1/T2","T3/T4"))
  temp$N_stage[grep("N1|N2",temp$N_stage)] <- "N1/N2"
  temp$N_stage <- factor(temp$N_stage,levels=c("N0","N1/N2"))
  temp$M_stage <- factor(temp$M_stage,levels=c("M0","M1"))
  temp$TNM_stage[grep("Stage 1|Stage 2",temp$TNM_stage)] <- "Stage 1/Stage 2"
  temp$TNM_stage[grep("Stage 3|Stage 4",temp$TNM_stage)] <- "Stage 3/Stage 4"
  temp$TNM_stage <- factor(temp$TNM_stage,levels=c("Stage 1/Stage 2","Stage 3/Stage 4"))
  # temp$Residual[grep("R1|R2",temp$Residual)] <- "R1/R2"
  # temp$Residual <- factor(temp$Residual,levels=c("R0","R1/R2"))
  med <- 67.5
  tempage <- ifelse(temp$age < 67.5,"<67.5",">=67.5")
  temp$age <- tempage
  temp$age <- factor(temp$age,levels=c("<67.5",">=67.5"))
  #temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
  temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
  temp
});
names(COX.DFS.byDriverGene.SingleMarker) <- DFS.significant;

#进行DFS COX分析 得到结果COX.DFS.result
COX.DFS.result <- lapply(COX.DFS.byDriverGene.SingleMarker,function(x){
	Cox.function(x);
})
names(COX.DFS.result) <- DFS.significant;
save(COX.DFS.result,file="D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.DFS.result.RData");

## 进行DFS 多个基因组合多因素cox分析
{
  length(COX.DFS.byDriverGene.SingleMarker)
  names(COX.DFS.byDriverGene.SingleMarker)
  head(COX.DFS.byDriverGene.SingleMarker[[1]])
  library(dplyr)
  temp <- full_join(COX.DFS.byDriverGene.SingleMarker[[1]],COX.DFS.byDriverGene.SingleMarker[[2]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site"))
  for(i in names(COX.DFS.byDriverGene.SingleMarker)[c(-1,-2)]){
    temp <- full_join(temp,COX.DFS.byDriverGene.SingleMarker[[i]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site"))
  }
  colnames(temp)[c(12,13,14,15)] <- DFS.significant
  head(temp)
  COX.DFS.byDriverGene.SingleMarker.allgenes <- temp
}
COX.DFS.allgenes.result <- Cox.function(COX.DFS.byDriverGene.SingleMarker.allgenes)
Cox.function(na.omit(COX.DFS.byDriverGene.SingleMarker.allgenes))
save(COX.DFS.allgenes.result,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.DFS.allgenes.result.RData")

# 进行Test.proportional.hazards分析
{
  ## 只对基因突变克隆状态进行等比例风险分析
  uni.proportional.hazards <- list()
  for(i in names(COX.DFS.allgenes.result)){
    uni.proportional.hazards[[i]] <- Test.proportional.hazards(COX.DFS.allgenes.result[[i]][,c(1,2,3,12)],paste0("F:/",i,".pdf"))
  }
  ## 对所有临床字段进行等比例风险假设分析
  proportional.hazards <- list()
  for(i in names(COX.DFS.byDriverGene.SingleMarker.allgenes)){
    proportional.hazards[[i]] <- Test.proportional.hazards(COX.DFS.byDriverGene.SingleMarker.allgenes[[i]],paste0("F:/",i,".pdf"))
  }
  Test.proportional.hazards(COX.DFS.byDriverGene.SingleMarker.allgenes,paste0("D:/",i,".pdf"))
}

## 进行OS多个基因组合逐步cox回归分析
{
  ##提取COX模型的HR,CI,p值
  mvcoxres <- function(mod){
    hr.ci <- paste0(
      round(summary(mod)$conf.int[, 1], 2), " (",
      round(summary(mod)$conf.int[, 3], 2), ", ",
      round(summary(mod)$conf.int[, 4], 2), ")"
    )
    p <- round(summary(mod)$coefficients[, 5], 3)
    res <- data.frame(hr.ci, p, stringsAsFactors = F)
    res$p[res$p < 0.001] <- "<.001"
    colnames(res) <- c("**HR (95% CI)**", "**p-value**")
    return(res)
  }
  library(MASS)
  library(survival)
  head(COX.DFS.byDriverGene.SingleMarker.allgenes)
  coxmodel <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.DFS.byDriverGene.SingleMarker.allgenes)[-c(1:3)],collapse = "+"))),
                    data=na.omit(COX.DFS.byDriverGene.SingleMarker.allgenes)) #stepAIC 不允许生存数据有NA,考虑补缺失值
  summary(coxmodel)
  extractAIC(coxmodel)
  step.cox <- stepAIC(coxmodel,direction = "both",trace = F) #stepAIC 不允许生存数据有NA
  #step.cox.result <- summary(step.cox)
  stepwise.table <- mvcoxres(step.cox)
  save(stepwise.table,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.DFS.MutvsWT.stepwise.RData")
}


