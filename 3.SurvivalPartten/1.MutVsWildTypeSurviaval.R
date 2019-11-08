#############################################################
#     比较driver基因突变和不突变情况的生存情况              #
#                 包括KM生存曲线                            #
#############################################################
#导入临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData"); #colorectal.clinic.data
#导入克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData");#colorectal.select.table 
#导入driver gene
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData"); 
#ColorectalDriverGenes
##提取需要的临床字段c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29)
{
  colnames(colorectal.clinic.data)
  head(colorectal.clinic.data,1)
  head(colorectal.clinic.data[,c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29,30)])
  colorectal.clinic.data <- colorectal.clinic.data[,c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29,30)]
  head(colorectal.clinic.data)
}
#对突变数据进行清洗
colnames(colorectal.select.table)[3] <- "Gene_Symbol";
colnames(colorectal.select.table)[1] <- "Patient_ID";
#对临床数据进行清洗
colnames(colorectal.clinic.data)[1] <- "Patient_ID";

#函数一 grouping.bymutation.singleMarker
#' @description 依据单个基因是否突变对样本进行分类
#' @param markers marker基因的名字
#' @param clinical.data  干净的临床数据，每行代表一个样本，每列代表样本的不同临床信息
#' @param mutation.data  干净的突变数据，每行代表一个突变，每列代表突变的相关信息
#' @return 返回的是带有标签的临床数据，最后一列sample.label为分类标签，是通过这个函数产生的，必须为factor
grouping.bymutation.singleMarker <- function(markers, 
                                             clinical.data, 
                                             mutation.data
)
{
  #STEP 1:处理样本，保留在突变数据和临床数据中共同出现的样本
  common.sample <- intersect( unique(mutation.data$Patient_ID), clinical.data$Patient_ID ) #取样本的交集
  clinical.data <- clinical.data[clinical.data$Patient_ID %in% common.sample, ]
  mutation.data <- mutation.data[mutation.data$Patient_ID %in% common.sample, ] 
  #STEP 2:处理markers，要求markers必须有效
  markers <- intersect(markers ,unique(mutation.data$Gene_Symbol))
  if(length(markers) == 0)
  {
    #stop("所选则的marker并未出现在突变谱中！");
  	sample.label <- rep(paste(markers,"WT"),length(common.sample))
  	return(cbind.data.frame(clinical.data,sample.label));
  }  
  #STEP 3:对样本进行分类，根据是否突变将样本分为突变型和野生型
  label <- c(paste(markers, "WT",sep = " "),paste(markers, "mut",sep = " "))
  #提取突变样本id
  mut.sample <- unique(mutation.data[mutation.data$Gene_Symbol %in% markers,"Patient_ID"]) 
  clinical.data$'sample.label' <- factor(label[as.numeric(clinical.data$Patient_ID %in% mut.sample)+1], levels = label)  #给样本打标签sample.label
  return(clinical.data)
}
#生成97个基因依据突变状态分类的临床数据
Survival.byDriverGene.SingleMarker <- lapply(ColorectalDriverGenes,function(x){
    grouping.bymutation.singleMarker(x,colorectal.clinic.data,colorectal.select.table)
})
names(Survival.byDriverGene.SingleMarker) <- ColorectalDriverGenes
#存储OS和DFS生存数据
save(Survival.byDriverGene.SingleMarker,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/Survival.byDriverGene.SingleMarker.RData");
lapply(Survival.byDriverGene.SingleMarker,function(x){
	x
});

#######################################################################################################
#进行KM曲线OS分析，对临床数据表头名字进行清洗，满足函数接口
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/Survival.byDriverGene.SingleMarker.RData")
#进行KM曲线OS分析，对临床数据表头名字进行清洗，满足函数接口
head(Survival.byDriverGene.SingleMarker[[1]])
OS.byDriverGene.SingleMarker <- lapply(Survival.byDriverGene.SingleMarker,function(x){
  x <- x[,c(-12,-13,-14,-16)] #对临床数据进行筛选
  colnames(x) <- c("Patient_ID","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","Residual","age","event","time","Subtype","Tumour_site","sample.label")
  #事件转换为数字类型
  x$event <- as.numeric(x$event)
  #年龄依据中值离散化处理(没问题)
  #med <- median(x$age)
  med <- 67.5
  tempage <- vector(mode="character",length=length(x$age))
  tempage[which(x$age >= med)] <- paste0(">=",med)
  tempage[which(x$age < med)] <- paste0("<",med)
  x$age <- tempage
  return(x) #如果不设置，默认返回最后一行
})
names(OS.byDriverGene.SingleMarker) <- ColorectalDriverGenes
head(OS.byDriverGene.SingleMarker[[1]])
#存储OS生存数据
save(OS.byDriverGene.SingleMarker,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.SingleMarker.RData")
#测试
plot.surv(clinical.data = OS.byDriverGene.SingleMarker[[27]],upper.time=NULL, unit.xlabel = "month",
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                      main = NULL,survival.event = "Overall Survival",inputFilePath="F:/Rsources/Project/DriverMutation.KMplot/",picture.name = names(OS.byDriverGene.SingleMarker)[27]);
#OS KM分析
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.SingleMarker.RData")
medianOS_survivalRate <- list() #中位生存期和5年生存率
for(i in names(OS.byDriverGene.SingleMarker)){
  print(i);
  if(length(table(OS.byDriverGene.SingleMarker[[i]]$sample.label)) == 2 & min(table(OS.byDriverGene.SingleMarker[[i]]$sample.label)) > 5)
  medianOS_survivalRate[[i]] <- plot.surv(clinical.data = OS.byDriverGene.SingleMarker[[i]],upper.time=NULL, unit.xlabel = "month",
                      risk.table = TRUE, pval = TRUE, conf.int = F,surv.median.line="hv",
                      main = i,survival.event = "Overall Survival",inputFilePath="D:/Rsources/Project/DriverMutation.KMplot/",picture.name = i);
}
length(medianOS_survivalRate) # 9
save(medianOS_survivalRate,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/MutvsWT.medianOS_survivalRate.RData")
#生存显著的有AKT1,AXIN2,BCLAF1,CASP8,CD70,CNBD1,EPHA7,GRIN2A,SMAD2
###############################################################################################################

###############################################################################################################
load(file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/Survival.byDriverGene.SingleMarker.RData")
#进行KM曲线DFS分析，对临床数据表头名字进行清洗，满足函数接口
DFS.byDriverGene.SingleMarker <- lapply(Survival.byDriverGene.SingleMarker,function(x){
	colnames(x)[12] <- "event";
  colnames(x)[13]	<- "time";
  x[,12] <- as.numeric(x[,12]);
	x[,13] <- as.numeric(x[,13]);
	return(x); #如果不设置，默认返回最后一行
});
names(DFS.byDriverGene.SingleMarker) <- ColorectalDriverGenes
#存储DFS生存数据
save(DFS.byDriverGene.SingleMarker,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.byDriverGene.SingleMarker.RData")
#DFS KM分析
medianDFS_survivalRate <- list()
for(i in names(DFS.byDriverGene.SingleMarker)){
  print(i);
  if(length(table(DFS.byDriverGene.SingleMarker[[i]]$sample.label)) == 2 & min(table(DFS.byDriverGene.SingleMarker[[i]]$sample.label)) > 5)
    medianDFS_survivalRate[[i]] <- plot.surv(clinical.data = DFS.byDriverGene.SingleMarker[[i]],upper.time=NULL, unit.xlabel = "month",
                                             risk.table = TRUE, pval = TRUE, conf.int = F,surv.median.line="hv",
                                             main = i,survival.event = "Progress Free Survival",inputFilePath="F:/Rsources/Project/DFS.DriverMutation.KMplot/",picture.name = i);
}

#生存显著的AXIN2,CUX1,KRAS(至少在10个样本中突变)
