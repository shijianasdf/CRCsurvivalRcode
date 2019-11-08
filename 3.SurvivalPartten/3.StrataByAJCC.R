#################################################################################################
#################依据AJCC_stage分层分析，分析依据driver基因突变不突变状态下的KM##################
#################################################################################################
#导入OS生存数据
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.SingleMarker.RData");
#导入DFS生存数据
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.byDriverGene.SingleMarker.RData");
##基于AJCC_stage分4层## stage1:95  stage2:196  stage3:153  stage4:73
AJCC.OS.survival.data <- lapply(OS.byDriverGene.SingleMarker,function(x){
	split(x,x$AJCC_PATHOLOGIC_TUMOR_STAGE);
});
names(AJCC.OS.survival.data) <- names(OS.byDriverGene.SingleMarker);

#KM生存曲线如果有1组样本数太少会报错，一般至少为5
for(i in names(AJCC.OS.survival.data)){
  print(i);
	for(j in 1:length(AJCC.OS.survival.data[[i]])){
	    inputFilePath <- paste0(file.path("D:/R/Project/StrataByAJCC",i,"AJCC_stage",unique(AJCC.OS.survival.data[[i]][[j]]$AJCC_PATHOLOGIC_TUMOR_STAGE)),"/");
		#要有两类样本并且每类样本个数都要大于5
		if(length(table(AJCC.OS.survival.data[[i]][[j]]$sample.label)) == 2 & min(table(AJCC.OS.survival.data[[i]][[j]]$sample.label)) > 5){		    
			plot.surv(clinical.data = AJCC.OS.survival.data[[i]][[j]],upper.time=NULL, unit.xlabel = "month",
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE, surv.median.line="hv", 
                      main = NULL,survival.event = "Overall Survival",inputFilePath=inputFilePath,picture.name = i);
		}else{
		    next;
		}
	}	
}
#得到87个基因,其中AKT1在stage2中显著,AXIN2在stage2中显著，BCLAF1在stage2中显著，CASP8在stage2中显著，EPHA7在stage2中显著

##基于AJCC_stage分4层## stage1:95  stage2:196  stage3:153  stage4:73
AJCC.DFS.survival.data <- lapply(DFS.byDriverGene.SingleMarker,function(x){
	split(x,x$AJCC_PATHOLOGIC_TUMOR_STAGE);
});
names(AJCC.DFS.survival.data ) <- names(DFS.byDriverGene.SingleMarker);

#KM生存曲线如果有1组样本数太少会报错，一般至少为5
for(i in names(AJCC.DFS.survival.data)){
    print(i);
	for(j in 1:length(AJCC.DFS.survival.data[[i]])){
	    inputFilePath <- paste0(file.path("D:/R/Project/StrataByAJCC",i,"AJCC_stage",unique(AJCC.OS.survival.data[[i]][[j]]$AJCC_PATHOLOGIC_TUMOR_STAGE)),"/");
		#要有两类样本并且每类样本个数都要大于5
		if(length(table(AJCC.DFS.survival.data[[i]][[j]]$sample.label)) == 2 & min(table(AJCC.DFS.survival.data[[i]][[j]]$sample.label)) > 5){		    
			plot.surv(clinical.data = AJCC.DFS.survival.data[[i]][[j]],upper.time=NULL, unit.xlabel = "month",
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE, surv.median.line="hv", 
                      main = NULL,survival.event = "Progress Free Survival",inputFilePath=inputFilePath,picture.name = i);
		}else{
		    next;
		}
	}	
}
#得到87个基因，其中AXIN2在stage3中显著，CUX1在stage2中显著

