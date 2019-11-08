#------------------------------------------------------------------
#' 
#'    预后模型的比较(比较引进克隆状态和普通变异的预后模型)
#' 
#------------------------------------------------------------------
########导入OS COX分析所需数据#########
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.ClonalSubclonal.RData")
########导入OS生存数据##########                                                                          
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.SingleMarker.RData")

##处理清洗OS COX生存数据
#处理清洗OS.SubclonalClonal cox数据
{
  OS.SubclonalClonal.significant <- c("ANK1","ARID1A","CASP8","GRIN2A","SMAD2")
  pos <- match(OS.SubclonalClonal.significant,names(OS.byDriverGene.ClonalSubclonal))
  COX.OS.byDriverGene.ClonalSubclonal <- vector(mode="list",length=length(OS.SubclonalClonal.significant))
  for(i in 1:length(pos)){
    COX.OS.byDriverGene.ClonalSubclonal[[i]] <- OS.byDriverGene.ClonalSubclonal[[pos[i]]]
  }
  names(COX.OS.byDriverGene.ClonalSubclonal) <- OS.SubclonalClonal.significant
  head(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  levels(COX.OS.byDriverGene.ClonalSubclonal[[1]]$sample.label)
  ## 对临床字段进行处理
  COX.OS.byDriverGene.ClonalSubclonal <- lapply(COX.OS.byDriverGene.ClonalSubclonal,function(x){
    #temp <- x[,c(1,10,11,2,3,4,5,6,7,8,9,12,13,14)] #所有临床变量
    temp <- x[,c(1,10,11,2,3,4,5,6,7,9,13,14)] #剔除亚型信息
    #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Subtype","Tumour_site","Mutation_status")
    #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Tumour_site","Mutation_status")
    temp$time <- as.numeric(temp$time)
    temp$event <- as.numeric(temp$event)
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
  names(COX.OS.byDriverGene.ClonalSubclonal) <- OS.SubclonalClonal.significant;
  head(COX.OS.byDriverGene.ClonalSubclonal[[1]])
}
#处理清洗OS.Mut vs WT  cox数据
{
  OS.MutvsWT.significant <- c("AKT1","AXIN2","BCLAF1","CASP8","CNBD1","CD70","EPHA7","GRIN2A","SMAD2","ANK1","ARID1A")
  OS.MutvsWT.SubclonalClonal.significant <- intersect(OS.MutvsWT.significant,OS.SubclonalClonal.significant)
  pos <- match(OS.MutvsWT.SubclonalClonal.significant,names(OS.byDriverGene.SingleMarker))
  COX.OS.byDriverGene.SingleMarker <- vector(mode="list",length=length(OS.MutvsWT.SubclonalClonal.significant))
  for(i in 1:length(pos)){
    COX.OS.byDriverGene.SingleMarker[[i]] <- OS.byDriverGene.SingleMarker[[pos[i]]]
  }
  names(COX.OS.byDriverGene.SingleMarker) <- OS.MutvsWT.SubclonalClonal.significant
  head(COX.OS.byDriverGene.SingleMarker[[1]])
  COX.OS.byDriverGene.SingleMarker <- lapply(COX.OS.byDriverGene.SingleMarker,function(x){
    #temp <- x[,c(1,22,23,2,6,7,8,9,19,21,26)];
    #colnames(temp) <- c("Patient_ID","event","time","sex","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site","Mutation_status");
    temp <- x[,c(1,10,11,2,3,4,5,6,7,9,13,14)] #剔除亚型信息
    #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Subtype","Tumour_site","Mutation_status")
    #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Tumour_site","Mutation_status")
    temp$time <- as.numeric(temp$time)
    temp$event <- as.numeric(temp$event)
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
  names(COX.OS.byDriverGene.SingleMarker) <- OS.MutvsWT.SubclonalClonal.significant
  head(COX.OS.byDriverGene.SingleMarker[[1]])
}
## 融合上面两个临床数据框，为后面模型比较做准备OS.MutvsWT.SubclonalClonal.overlap
{
  library(dplyr)
  OS.MutvsWT.SubclonalClonal.significant #交集基因
  OS.MutvsWT.SubclonalClonal.overlap <- list()
  for(i in OS.MutvsWT.SubclonalClonal.significant){
    OS.MutvsWT.SubclonalClonal.overlap[[i]] <- full_join(x=COX.OS.byDriverGene.SingleMarker[[i]],
                                                         y=COX.OS.byDriverGene.ClonalSubclonal[[i]],
                                                         by=c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site"))
    colnames(OS.MutvsWT.SubclonalClonal.overlap[[i]])[c(12,13)] <- c("MutvsWT","SubclonalClonal") 
  }
  head(OS.MutvsWT.SubclonalClonal.overlap$CASP8)
  table(OS.MutvsWT.SubclonalClonal.overlap$CASP8$SubclonalClonal)
}


##---------------------------------------------------------------------------------------------
#  calculate the C-index for a cox model
#'@description 计算生存模型的C-index或者两个模型比较的p值
#'@param clinical.data: 临床变量信息，至少包含model中的变量信息
#'@param time: 数值型向量，患者对应的生存时间
#'@param event: 患者对应的生存状态，通常0=alive, 1=dead
#'@param models：列表，每一个元素是一个字符型向量，包含一个model所有变量对应的列名,至少包含两个元素
#'@param diff: 逻辑值，diff = FALSE不计算两个模型比较的p值，只返回对应的C-index值；diff = TRUE返回模型对应的C-index值及模型比较的p值
#'@return 返回一个data.frame，包含C-index值，置信区间，比较的p值（当diff = T时）
cIndex <- function(clinical.data, time, event, models, diff = TRUE){
  options(stringsAsFactors = FALSE)
  # load survival, Hmisc package
  suppressPackageStartupMessages(require(survival))
  suppressPackageStartupMessages(require(Hmisc))
  ##-------------------------------
  # CoxphObject：create a coxph object
  CoxphObject <- function(data, time, event, variables){
    if (length(variables) == 1) {
      formula <- as.formula(paste("Surv(time, event)~", variables))
    } else {
      formula <- as.formula(paste("Surv(time, event)~", paste(variables, collapse=" + ")))
    }      
    coxph.object <-coxph(formula, data = data)
    
    return(coxph.object)
  }
  ##-------------------------------
  # 输出的Cindex和置信区间
  c.ci <- function(rcorrobj){
    CIndex <- round(rcorrobj['C Index'], digits = 4)
    se     <- rcorrobj['S.D.']/2
    Lower <- round(CIndex - 1.96*se, digits = 4)
    Upper <- round(CIndex + 1.96*se, digits = 4)
    result <- c(CIndex, Lower, Upper)
    names(result) <- c("C-Index", "Lower", "Upper")
    
    return(result)
  }
  # 计算每个模型的C-index值和置信区间
  coxph.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})
  pred.models.coxph <- lapply(coxph.list, function(x){predict(x, type = "risk")})
  models.result <- lapply(pred.models.coxph, function(x){
              #rcorr.cens(-x, Surv(time = time, event = event))
              1-rcorr.cens(x, Surv(time = time, event = event))
    })
  models.filter.result <- lapply(models.result, function(x){c.ci(rcorrobj = x)})
  # 是否进行C-index1的比较
  if (diff == FALSE) {
    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval)
    colnames(result) <- c('C-index', '95%CI')
    return(result)
  }else{
    # 计算比较的p值，都是和第一个模型比较
    compare.cindex <- lapply(pred.models.coxph[-1], function(x){rcorrp.cens(pred.models.coxph[[1]], x, Surv(time = time, event = event))})
    p.value <- c("-", unlist(lapply(compare.cindex, function(x)(round(2*(1 - pnorm(abs(x[['Dxy']] / x[['S.D.']]))), digits=4))))) 
    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval, p.value)
    colnames(result) <- c('C-index', '95%CI', 'p_value')
    return(result)
  }
}
tempList <- list()
for(i in names(OS.MutvsWT.SubclonalClonal.overlap))
   tempList[[i]] <- OS.MutvsWT.SubclonalClonal.overlap[[i]][which(!is.na(OS.MutvsWT.SubclonalClonal.overlap[[i]]$SubclonalClonal)),]
table(tempList[[3]]$SubclonalClonal,useNA = "ifany")
head(tempList[[3]])
names(tempList)

for(i in 1:5){
  tt <- cIndex(clinical.data=na.omit(tempList[[i]]),
         time=na.omit(tempList[[i]])$time,
         event=na.omit(tempList[[i]])$event,
         models=list(c("sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site","MutvsWT"),c("sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site","SubclonalClonal")))
  print(tt)
}

for(i in 1:5){
  tt <- cIndex(clinical.data=na.omit(tempList[[i]]),
               time=na.omit(tempList[[i]])$time,
               event=na.omit(tempList[[i]])$event,
               models=list(c("MutvsWT"),c("SubclonalClonal")))
  print(tt)
}

