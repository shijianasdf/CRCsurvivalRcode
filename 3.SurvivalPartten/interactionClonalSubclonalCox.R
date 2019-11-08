##########导入OS COX分析所需数据#########
load("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/2.SurvivalKM/OS.byDriverGene.ClonalSubclonal.RData");
##########导入DFS COX分析所需数据#########
load("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/2.SurvivalKM/DFS.byDriverGene.ClonalSubclonal.RData");
COX.OS.byDriverGene.ClonalSubclonal <- lapply(OS.byDriverGene.ClonalSubclonal,function(x){
	 temp <- x[,c(1,22,23,2,6,7,8,9,19,21,26)];
	 colnames(temp) <- c("Patient_ID","event","time","sex","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site","Mutation_status");
	 #年龄离散化
	 med <- median(temp$age);
	 tempage <- vector(mode="logical",length=length(temp$age));
	 tempage[which(temp$age >= med)] <- paste0(">=",med);
	 tempage[which(temp$age < med)] <- paste0("<",med);
	 temp$age <- tempage;
	 #对协变量进行factor转换，控制levels的第一项为参考项,目前默认第一个level为参考项 
	 temp;
});
tempList <- lapply(COX.OS.byDriverGene.ClonalSubclonal,function(x){
	 as.character(x$Mutation_status);
});
temp <- as.data.frame(do.call(cbind,tempList))

#画互斥图
library(pheatmap)
tt <- t(cbind(c(1,1,0,0,1,0,0,1,0,1,1,1,0,0,1,1,0),c(1,1,0,0,1,0,0,1,0,1,1,1,0,0,1,1,0)))
pheatmap(tt)

#interaction test
interaction_table <- cbind.data.frame(COX.OS.byDriverGene.ClonalSubclonal[[1]][,-11],temp)
interaction_table$APC <- factor(interaction_table$APC,levels=c("APC WT","APC Clonal","APC Subclonal"))
interaction_table$TP53 <- factor(interaction_table$TP53,levels=c("TP53 WT","TP53 Clonal","TP53 Subclonal"))
library(survival);
coxmodel <- coxph(Surv(time, event) ~ sex + T_stage + N_stage + M_stage + AJCC_stage + age + Tumour_site + APC*TP53, data = interaction_table);
tmp <- summary(coxmodel);


    
    #提取p值，保留两位有效数字
    p.value <- round(tmp$coefficients[ ,5], digits = 4);
    p.value[which(p.value < 0.0001)] <- "<0.0001"; 
    
    #提取beta值，这里得到的coefficients为矩阵，且有多行（每行对应一个协变量）
    #beta <- round(tmp$coefficients[ ,1], digits = 4);
    
    #提取风险比
    HR <- round(tmp$coefficients[ ,2], digits = 4);
    
    #提取95%置信区间上下界
    HR.confint.lower <- round(tmp$conf.int[ ,"lower .95"], digits = 4);
    HR.confint.upper <- round(tmp$conf.int[ ,"upper .95"], digits = 4);
    
    #合并风险比HR和置信区间为一个内容
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
	
	variate <- rownames(tmp$coefficients);
    
    #整合输出内容为data.frame
    multiv.result <- as.data.frame(cbind(variate, HR, p.value));
    colnames(multiv.result) <- c("variate", "multiv HR (95% CI for HR)", "multiv p value");