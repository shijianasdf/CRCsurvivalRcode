##################################################################
#          克隆亚克隆3种状态的OS，DFS分层分析
###################################################################
#############导入OS 生存数据##################
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.ClonalSubclonal.RData");
#############导入DFS 生存数据#################
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.byDriverGene.ClonalSubclonal.RData");

##基于微卫星不稳定性状态分2层 MSS/MSI
##对临床字段进行处理
strata.OS.byDriverGene.ClonalSubclonal <- lapply(OS.byDriverGene.ClonalSubclonal,function(x){
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
#OS.significant <- c("ANK1","ARID1A","GRIN2A","CASP8","SMAD2");
TS.OS.survival.data <- lapply(strata.OS.byDriverGene.ClonalSubclonal,function(x){
  split(x,x$Tumour_site)
})
names(TS.OS.survival.data) <- names(OS.byDriverGene.ClonalSubclonal)
#AJCC.OS.survival.data <- AJCC.OS.survival.data[match(OS.significant,names(AJCC.OS.survival.data))];
#names(AJCC.OS.survival.data) <- c("ANK1","ARID1A","GRIN2A","CASP8","SMAD2");
#KM生存曲线如果有1组样本数太少会报错，一般至少为5
TS_medianOS_survivalRate <- list()
for(i in names(TS.OS.survival.data)){
  print(i);
  for(j in 1:length(TS.OS.survival.data[[i]])){
    inputFilePath <- paste0(file.path("F:/Rsources/Project/StrataByTS",i,"TumourSite",unique(TS.OS.survival.data[[i]][[j]]$Tumour_site)),"/");
    #要有三类样本
    if(length(table(TS.OS.survival.data[[i]][[j]]$sample.label)) >= 2 & min(table(TS.OS.survival.data[[i]][[j]]$sample.label)) >= 5){		    
      TS_medianOS_survivalRate[[i]] <- plot.surv(clinical.data = TS.OS.survival.data[[i]][[j]],upper.time=NULL, unit.xlabel = "month",
                                                  risk.table = TRUE, pval = TRUE, conf.int = F, surv.median.line="hv", 
                                                  main = i,survival.event = "Overall Survival",inputFilePath=inputFilePath,picture.name = i);
    }else{
      next;
    }
  }	
}
#只有"ANK1","ARID1A"，"CASP8"在stage_2中显著


##基于AJCC_stage分两层
#DFS.significant <- c("ARID1A","KRAS");
strata.DFS.byDriverGene.ClonalSubclonal <- lapply(DFS.byDriverGene.ClonalSubclonal,function(x){
  #temp <- x[,c(1,10,11,2,3,4,5,6,7,8,9,12,13,14)] #所有临床变量
  temp <- x[,c(1,9,10,2,3,4,5,6,7,8,11,12)] #剔除亚型信息
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
TS.DFS.survival.data <- lapply(strata.DFS.byDriverGene.ClonalSubclonal,function(x){
  split(x,x$Tumour_site)
})
names(TS.DFS.survival.data) <- names(DFS.byDriverGene.ClonalSubclonal)
#AJCC.DFS.survival.data <- AJCC.DFS.survival.data[match(DFS.significant,names(AJCC.DFS.survival.data))];
#names(AJCC.DFS.survival.data) <- DFS.significant;
#KM生存曲线如果有1组样本数太少会报错，一般至少为5
TS_medianDFS_survivalRate <- list()
for(i in names(TS.DFS.survival.data)){
  print(i);
  for(j in 1:length(TS.DFS.survival.data[[i]])){
    inputFilePath <- paste0(file.path("F:/Rsources/Project/StrataByTS",i,"TumourSite",unique(TS.DFS.survival.data[[i]][[j]]$Tumour_site)),"/");
    #要有两类样本并且每类样本个数都要大于5
    if(length(table(TS.DFS.survival.data[[i]][[j]]$sample.label)) >= 2 & min(table(TS.DFS.survival.data[[i]][[j]]$sample.label)) > 5){		    
      TS_medianDFS_survivalRate[[i]] <- plot.surv(clinical.data = TS.DFS.survival.data[[i]][[j]],upper.time=NULL, unit.xlabel = "month",
                                                   risk.table = TRUE, pval = TRUE, conf.int = F, surv.median.line="hv", 
                                                   main = i,survival.event = "Progress Free Survival",inputFilePath=inputFilePath,picture.name = i);
    }else{
      next;
    }
  }	
}

