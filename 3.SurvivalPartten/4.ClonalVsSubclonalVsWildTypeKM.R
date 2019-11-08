##########################################################################
#               DriverGene  3种状态的KM曲线分析
##########################################################################

########导入生存数据,依据克隆状态对临床数据进行分类##############
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.survival.data.list.RData")
#CRC.survival.data.list

#进行KM曲线OS分析，对临床数据表头名字进行清洗，满足函数接口
head(CRC.survival.data.list[[1]])
OS.byDriverGene.ClonalSubclonal <- lapply(CRC.survival.data.list,function(x){
  x <- x[,c(-12,-13,-14,-16)] #对临床数据进行筛选
  colnames(x) <- c("Patient_ID","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","Residual","age","event","time","Subtype","Tumour_site","sample.label")
  #事件转换为数字类型
  x$event <- as.numeric(x$event)
  #年龄依据中值离散化处理(没问题)
  #med <- median(x$age)
  # med <- 67.5
  # tempage <- vector(mode="character",length=length(x$age))
  # tempage[which(x$age >= med)] <- paste0(">=",med)
  # tempage[which(x$age < med)] <- paste0("<",med)
  # x$age <- tempage
	return(x) #如果不设置，默认返回最后一行
})
names(OS.byDriverGene.ClonalSubclonal) <- names(CRC.survival.data.list)
head(OS.byDriverGene.ClonalSubclonal[[1]])
class(OS.byDriverGene.ClonalSubclonal[[1]]$event)
class(OS.byDriverGene.ClonalSubclonal[[1]]$time)
#存储三种状态分类的OS生存数据
save(OS.byDriverGene.ClonalSubclonal,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.ClonalSubclonal.RData");
#OS KM分析,筛选亚克隆和克隆突变数>5的样本,如果选10就一个显著的都没有，画出KM曲线并且返回中位生存期CI和5年生存率
medianOS_survivalRate <- list() #中位生存期和5年生存率
for(i in names(OS.byDriverGene.ClonalSubclonal)){
  print(i);
  #if(length(table(OS.byDriverGene.ClonalSubclonal[[i]]$sample.label)) == 3 & min(table(OS.byDriverGene.ClonalSubclonal[[i]]$sample.label)) > 5)
  if(length(table(OS.byDriverGene.ClonalSubclonal[[i]]$sample.label)) >= 2 & min(table(OS.byDriverGene.ClonalSubclonal[[i]]$sample.label)) > 5)
    medianOS_survivalRate[[i]] <- plot.surv(clinical.data = OS.byDriverGene.ClonalSubclonal[[i]],upper.time=NULL, unit.xlabel = "month",
                                        risk.table = TRUE, pval = TRUE, conf.int = FALSE,surv.median.line="hv",
                                        main = i,survival.event = "Overall Survival",inputFilePath="D:/R/Project/DriverMutation.KMplot/",picture.name = i);
}
table(OS.byDriverGene.ClonalSubclonal[["ARID1A"]]$sample.label)
table(OS.byDriverGene.ClonalSubclonal[["ANK1"]]$sample.label)
table(OS.byDriverGene.ClonalSubclonal[["SMAD2"]]$sample.label)
table(OS.byDriverGene.ClonalSubclonal[["CASP8"]]$sample.label)
table(OS.byDriverGene.ClonalSubclonal[["GRIN2A"]]$sample.label)
table(OS.byDriverGene.ClonalSubclonal[["ARID1A"]]$sample.label)

#生存显著的有ANK1,ARID1A,GRIN2A,CASP8,SMAD2
#进行KM曲线DFS分析，对临床数据表头名字进行清洗，满足函数接口
DFS.byDriverGene.ClonalSubclonal <- lapply(CRC.survival.data.list,function(x){
  x <- x[,c(-8,-10,-11,-14,-15,-16)] #对临床数据进行筛选
  colnames(x) <- c("Patient_ID","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","event","time","Tumour_site","sample.label")
  #事件转换为数字类型
  x$event <- as.numeric(x$event)
  #年龄依据中值离散化处理(没问题)
  # med <- median(x$age)
  med <- 67.5
  tempage <- vector(mode="character",length=length(x$age))
  tempage[which(x$age >= med)] <- paste0(">=",med)
  tempage[which(x$age < med)] <- paste0("<",med)
  x$age <- tempage
  return(x) #如果不设置，默认返回最后一行
});
names(DFS.byDriverGene.ClonalSubclonal) <- names(CRC.survival.data.list);
#存储三种状态分类的OS生存数据
save(DFS.byDriverGene.ClonalSubclonal,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.byDriverGene.ClonalSubclonal.RData");
#DFS KM分析,筛选亚克隆和克隆突变数>5的样本
medianDFS_survivalRate <- list() #DFS中位生存期和5年生存率
for(i in names(DFS.byDriverGene.ClonalSubclonal)){
   print(i);
   if(length(table(DFS.byDriverGene.ClonalSubclonal[[i]]$sample.label)) >= 2 & min(table(DFS.byDriverGene.ClonalSubclonal[[i]]$sample.label)) > 5)
     medianDFS_survivalRate[[i]] <- plot.surv(clinical.data = DFS.byDriverGene.ClonalSubclonal[[i]],upper.time=NULL, unit.xlabel = "month",risk.table = TRUE, pval = TRUE, conf.int = F,surv.median.line="hv",
                  main = i,survival.event = "Progress Free Survival",inputFilePath="F:/Rsources/Project/DriverMutation.KMplot/",picture.name = i);
}

#DFS生存显著的有ARID1A,KRAS,其中KRAS单独比不显著，整体显著
######################################
