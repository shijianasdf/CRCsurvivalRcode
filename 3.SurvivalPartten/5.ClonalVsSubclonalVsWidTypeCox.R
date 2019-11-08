####################################################################
#依据Driver基因3种状态的单因素COX和多因素COX
#####################################################################
##########导入OS COX分析所需数据#########
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.ClonalSubclonal.RData");
##########导入DFS COX分析所需数据#########
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.byDriverGene.ClonalSubclonal.RData");

#处理清洗OS COX生存数据
{
  OS.significant <- c("ANK1","ARID1A","CASP8","GRIN2A","SMAD2")
  pos <- match(OS.significant,names(OS.byDriverGene.ClonalSubclonal))
  COX.OS.byDriverGene.ClonalSubclonal <- vector(mode="list",length=length(OS.significant))
  for(i in 1:length(pos)){
    COX.OS.byDriverGene.ClonalSubclonal[[i]] <- OS.byDriverGene.ClonalSubclonal[[pos[i]]]
  }
  names(COX.OS.byDriverGene.ClonalSubclonal) <- OS.significant
  head(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  levels(COX.OS.byDriverGene.ClonalSubclonal[[1]]$sample.label)
  ## 对临床字段进行处理
  COX.OS.byDriverGene.ClonalSubclonal <- lapply(COX.OS.byDriverGene.ClonalSubclonal,function(x){
    #temp <- x[,c(1,10,11,2,3,4,5,6,7,8,9,12,13,14)] #所有临床变量
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
    #temp$Residual[grep("R1|R2",temp$Residual)] <- "R1/R2"
    #temp$Residual <- factor(temp$Residual,levels=c("R0","R1/R2"))
    temp$age <- factor(temp$age,levels=c("<67.5",">=67.5"))
    #temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
    temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
    temp
  })
  head(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  levels(COX.OS.byDriverGene.ClonalSubclonal[[1]]$TNM_stage)
  names(COX.OS.byDriverGene.ClonalSubclonal) <- OS.significant
  head(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  library(mice)
  md.pattern(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  library("VIM") #观察NA分布，发现Residual tumour非常多,所以删除Residual tumour
  dev.new();aggr(COX.OS.byDriverGene.ClonalSubclonal[[1]],prop=FALSE,numbers=TRUE)
  dim(na.omit(COX.OS.byDriverGene.ClonalSubclonal[[1]]))
}
# 进行Test.proportional.hazards分析
{
  ## 只对基因突变克隆状态进行等比例风险分析
  uni.proportional.hazards <- list()
  for(i in names(COX.OS.byDriverGene.ClonalSubclonal)){
    uni.proportional.hazards[[i]] <- Test.proportional.hazards(COX.OS.byDriverGene.ClonalSubclonal[[i]][,c(1,2,3,12)],paste0("F:/",i,".pdf"))
  }
  ## 对所有临床字段进行等比例风险假设分析
  proportional.hazards <- list()
  for(i in names(COX.OS.byDriverGene.ClonalSubclonal)){
    proportional.hazards[[i]] <- Test.proportional.hazards(COX.OS.byDriverGene.ClonalSubclonal[[i]],paste0("F:/",i,".pdf"))
  }
}
# 对每个基因单独进行OS单多因素COX分析 得到结果COX.OS.ClonalSubclonal.result
{
  COX.OS.ClonalSubclonal.result <- lapply(COX.OS.byDriverGene.ClonalSubclonal,function(x){
    Cox.function(x)
  })
  names(COX.OS.ClonalSubclonal.result) <- OS.significant
  save(COX.OS.ClonalSubclonal.result,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.OS.ClonalSubclonal.result.RData")
}

# 将所有基因合在一起做单多因素COX
{
  library(dplyr)
  COX.OS.all.gene <- COX.OS.byDriverGene.ClonalSubclonal[[1]] %>% 
    full_join(COX.OS.byDriverGene.ClonalSubclonal[[2]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site")) %>% 
    full_join(COX.OS.byDriverGene.ClonalSubclonal[[3]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site")) %>% 
    full_join(COX.OS.byDriverGene.ClonalSubclonal[[4]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site")) %>% 
    full_join(COX.OS.byDriverGene.ClonalSubclonal[[5]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site"))
  colnames(COX.OS.all.gene)[c(12,13,14,15,16)] <- OS.significant
  save(COX.OS.all.gene,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.OS.all.gene.RData")
  ## 对变量进行等比例风险回归假设进行检验
  Test.proportional.hazards(na.omit(COX.OS.all.gene),paste0("D:/","test",".pdf"))
  ## 单多因素COX分析
  Cox.function(COX.OS.all.gene) #不过滤掉NA 使用封装的Cox.function不应该过滤NA
  Cox.function(na.omit(COX.OS.all.gene)) #过滤掉NA，这是错的
  dim(COX.OS.all.gene)
  ## 5个基因的互斥共发生检验
  occur.mutual.matrix <- COX.OS.all.gene[,c(12,13,14,15,16)]
}

# 5个基因混在一起进行OS逐步COX回归分析 
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
  #pos <- complete.cases(COX.OS.all.gene[,c(2,3)])
  head(COX.OS.all.gene)
  coxmodel <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.OS.all.gene)[-c(1:3)],collapse = "+"))),
                    data=na.omit(COX.OS.all.gene)) #stepAIC 不允许生存数据有NA
  summary(coxmodel)
  extractAIC(coxmodel)
  step.cox <- stepAIC(coxmodel,direction = "both",trace = F) #stepAIC 不允许生存数据有NA
  #step.cox.result <- summary(step.cox)
  stepwise.table <- mvcoxres(step.cox)
  
}



## 处理清洗DFS COX生存数据
DFS.significant <- c("ARID1A","KRAS");
pos <- match(DFS.significant,names(DFS.byDriverGene.ClonalSubclonal));
COX.DFS.byDriverGene.ClonalSubclonal <- vector(mode="list",length=length(DFS.significant));
for(i in 1:length(pos)){
   COX.DFS.byDriverGene.ClonalSubclonal[[i]] <- DFS.byDriverGene.ClonalSubclonal[[pos[i]]];
}
names(COX.DFS.byDriverGene.ClonalSubclonal) <- DFS.significant;
head(COX.DFS.byDriverGene.ClonalSubclonal[[1]])
COX.DFS.byDriverGene.ClonalSubclonal <- lapply(COX.DFS.byDriverGene.ClonalSubclonal,function(x){
	 temp <- x[,c(1,9,10,2,3,4,5,6,7,8,11,12)];
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
	 #temp$Residual[grep("R1|R2",temp$Residual)] <- "R1/R2"
	 #temp$Residual <- factor(temp$Residual,levels=c("R0","R1/R2"))
	 temp$age <- factor(temp$age,levels=c("<67.5",">=67.5"))
	 #temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
	 temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
	 temp;
});
names(COX.DFS.byDriverGene.ClonalSubclonal) <- DFS.significant
head(COX.DFS.byDriverGene.ClonalSubclonal[[1]])

# 将所有基因合在一起做单多因素COX
{
  library(dplyr)
  COX.DFS.all.gene <- COX.DFS.byDriverGene.ClonalSubclonal[[1]] %>% 
    full_join(COX.DFS.byDriverGene.ClonalSubclonal[[2]],by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site")) 
  colnames(COX.DFS.all.gene)[c(12,13)] <- DFS.significant
  head(COX.DFS.all.gene)
  ## 对变量进行等比例风险回归假设进行检验
  Test.proportional.hazards(na.omit(COX.DFS.all.gene),paste0("F:/","test",".pdf"))
  ## 单多因素COX分析
  Cox.function(COX.DFS.all.gene) #不过滤掉NA
  Cox.function(COX.DFS.all.gene[,-13]) #不过滤掉NA
  #Cox.function(na.omit(COX.DFS.all.gene)) #过滤掉NA
  dim(COX.DFS.all.gene)
}

# 2个基因混在一起进行DFS逐步COX回归分析 
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
  #pos <- complete.cases(COX.OS.all.gene[,c(2,3)])
  head(COX.DFS.all.gene)
  coxmodel <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.DFS.all.gene)[-c(1:3)],collapse = "+"))),
                    data=na.omit(COX.DFS.all.gene)) #stepAIC 不允许生存数据有NA
  summary(coxmodel)
  extractAIC(coxmodel)
  step.cox <- stepAIC(coxmodel,direction = "both",trace = F) #stepAIC 不允许生存数据有NA
  #step.cox.result <- summary(step.cox)
  stepwise.table <- mvcoxres(step.cox)
}
save(COX.DFS.ClonalSubclonal.result,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/2.SurvivalKM/COX.DFS.ClonalSubclonal.result.RData");




