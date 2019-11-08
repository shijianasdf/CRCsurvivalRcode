#---------------------------------------
#     评估克隆性和临床变量的关系
#---------------------------------------
#####加载依据克隆状态分类的临床数据，用于克隆性和临床变量的关系#####
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.clonal.subclonal.data.list.RData");
##对CRC.clonal.subclonal.data.list进行处理，挑选要用的临床变量来分析
Correlation.table.list <- lapply(CRC.clonal.subclonal.data.list,function(x){
	temp <- x[,c(2,3,4,5,6,7,9,17,18)]
	colnames(temp) <- c("sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site",colnames(x)[18])
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
	temp$AJCC_stage[grep("Stage 1|Stage 2",temp$AJCC_stage)] <- "Stage 1/Stage 2"
	temp$AJCC_stage[grep("Stage 3|Stage 4",temp$AJCC_stage)] <- "Stage 3/Stage 4"
	temp$AJCC_stage <- factor(temp$AJCC_stage,levels=c("Stage 1/Stage 2","Stage 3/Stage 4"))
	tempage <- vector(mode="logical",length=length(temp$age))
	tempage[which(temp$age >= 67.5)] <- paste0(">=",67.5)
	tempage[which(temp$age < 67.5)] <- paste0("<",67.5)
	temp$age <- tempage
	temp$age <- factor(temp$age,levels=c("<67.5",">=67.5"))
	#temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
	temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
	temp
})
names(Correlation.table.list) <- names(CRC.clonal.subclonal.data.list)
head(Correlation.table.list[[1]])

##过滤掉克隆和亚克隆突变数小于5的基因以及没有突变的基因
length(Correlation.table.list)
min( table(Correlation.table.list[[1]][,9]) )
tt <- lapply(Correlation.table.list,function(x){
  if( min( table(x[,9]) ) >= 5 & length(table(x[,9])) >= 2){ #过滤没有突变的基因以及亚克隆突变数或者克隆突变数小于5的基因
    x
  }else{
    NULL
  }
})
names(tt) <- names(Correlation.table.list)
tt2 <- tt[which(!sapply(tt, is.null))] #去除列表中的NULL
length(tt2)

lapply(tt2,function(x){
   table(x[,9])
})
Correlation.table.list <- tt2

##建立chisq.list用于储存卡方检验和fisher'test需要的数据
Correlation.table.list <- lapply(Correlation.table.list,function(x){
	tablename <- colnames(x)
	chisq.fisher.list <- list()
	for(i in tablename[1:8]){
	   chisq.fisher.list[[paste(tablename[9],i,sep="&")]] <- unclass(table(x[,i],x[,tablename[9]]))
	}
	chisq.fisher.list
})
#names(Correlation.table.list) <- names(CRC.clonal.subclonal.data.list)
length(Correlation.table.list)
Correlation.table.list[[1]]
Filter.Correlation.table.list <- Correlation.table.list
Filter.Correlation.table.list[[1]] # list里面嵌套着list
Filter.Correlation.table.list <- unlist(Filter.Correlation.table.list,recursive=F,use.names=T)

##过滤掉不存在变异的基因RBM10, BCOR, DCAF12L2, ZMYM3以及只有clonal突变的基因AXIN1, POLE, SMAD3, NTHL1
# filter.genes <- c("RBM10", "BCOR", "DCAF12L2", "ZMYM3","AXIN1", "POLE", "SMAD3", "NTHL1") #过滤掉不存在变异的基因以及只有克隆突变的基因
# filter.genes <- c("RBM10", "BCOR", "DCAF12L2", "ZMYM3")   #过滤掉不存在变异的基因，但保留只有clonal突变的基因AXIN1, POLE, SMAD3, NTHL1
# filter.genes.positions <- match(filter.genes,names(Correlation.table.list))
# Filter.Correlation.table.list <- Correlation.table.list[-filter.genes.positions]
# Filter.Correlation.table.list <- unlist(Filter.Correlation.table.list,recursive=F,use.names=T)

#存储基因克隆状态和临床变量的列联表信息
save(Filter.Correlation.table.list,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/Filter.Correlation.table.list.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/Filter.Correlation.table.list.RData")
Filter.Correlation.table.list[grep("TP53",names(Filter.Correlation.table.list))]
Filter.Correlation.table.list[grep("BRAF",names(Filter.Correlation.table.list))] #5
Filter.Correlation.table.list[grep("BCL9L",names(Filter.Correlation.table.list))] #7
Filter.Correlation.table.list[grep("CASP8",names(Filter.Correlation.table.list))]



##非矫正，矫正，和fisher检验的计算，所有p值存在result中##
noncorrect.chisq.result <- lapply(Filter.Correlation.table.list, function(x){chisq.test(x, correct = F)})
correct.chisq.result <- lapply(Filter.Correlation.table.list, function(x){chisq.test(x, correct = T)})
fisher.result <- lapply(Filter.Correlation.table.list, function(x){fisher.test(x, workspace=100000000)})
fisher.p.result <-  lapply(fisher.result, function(x){x$p.value})
correct.p.result <- lapply(correct.chisq.result, function(x){x$p.value})
noncorrect.p.result <- lapply(noncorrect.chisq.result, function(x){x$p.value})
result <- cbind(names(noncorrect.chisq.result), do.call(rbind.data.frame, correct.p.result), do.call(rbind.data.frame, noncorrect.p.result), do.call(rbind.data.frame, fisher.p.result));
colnames(result) <- c("name", "correct.p", "noncorrect.p", "fisher.p");
save(result,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CloneAndVariables.RData");

#提取显著的基因
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CloneAndVariables.RData")
head(result)

result[grep("TP53",result$name),]
format(result[grep("TP53",result$name),4],scientific = F,digits = 2)
result[grep("BCL9L",result$name),]
format(result[grep("BCL9L",result$name),4],scientific = F,digits = 2)
result[grep("BRAF",result$name),]
format(result[grep("BRAF",result$name),4],scientific = F,digits = 2)

sigResult <- result[result$fisher.p <= 0.05,]
library(stringr)
?str_split
unique(str_split(as.character(sigResult$name),"\\.",simplify = T)[,1]) 
head(sigResult)
#微卫星不稳定性状态相关的基因
MSS_genes <- unique( str_split( as.character(sigResult[grep( "MSS_state",as.character(sigResult$name) ),]$name),"\\.",simplify = T)[,1] )
#肿瘤位置相关的基因
Tumour_site_genes <- unique( str_split(as.character(sigResult[grep( "Tumour_site",as.character(sigResult$name) ),]$name),"\\.",simplify = T)[,1] )
#AJCC_stage相关的基因
AJCC_stage_genes <- unique( str_split(as.character(sigResult[grep( "AJCC_stage",as.character(sigResult$name) ),]$name),"\\.",simplify = T)[,1] )
#sex相关的基因
Sex_genes <- unique( str_split(as.character(sigResult[grep( "sex",as.character(sigResult$name) ),]$name),"\\.",simplify = T)[,1] )
#N_stage相关的基因
N_stage_genes <- unique( str_split(as.character(sigResult[grep( "N_stage",as.character(sigResult$name) ),]$name),"\\.",simplify = T)[,1] )
#T_stage相关的基因
T_stage_genes <- unique( str_split(as.character(sigResult[grep( "T_stage",as.character(sigResult$name) ),]$name),"\\.",simplify = T)[,1] )
#M_stage相关的基因
M_stage_genes <- unique( str_split(as.character(sigResult[grep( "M_stage",as.character(sigResult$name) ),]$name),"\\.",simplify = T)[,1] )
#age相关的基因
Age_genes <- unique( as.character(sigResult[grep( "age",as.character(sigResult$name) ),]$name) )[grep( "stage",unique( as.character(sigResult[grep( "age",as.character(sigResult$name) ),]$name) ),invert=T )]

# 生成可视化热图的矩阵
z <- strsplit(as.character(result[,1]), split = "&");
dim(result)
p.result <- matrix(NA, ncol= 43, nrow = 8, dimnames = list(unique(do.call(rbind,lapply(z, function(x){x[2]}))), unique(do.call(rbind,lapply(z, function(x){x[1]})))));
for(i in rownames(p.result))
{
	for(j in colnames(p.result))
	{
		p.result[i,j] <- result[which(result$name == paste(j, i, sep = "&")),4] #使用的是fisher.test的结果
	}
}
colnames(p.result) <- sapply(strsplit(colnames(p.result),"\\."),function(x){unique(x)})

#统计多少
apply(p.result,1,function(x){
  length(which(x <= 0.05))/length(x)
})
# sex   MSS_state     T_stage     N_stage     M_stage  AJCC_stage         age Tumour_site 
# 0.09677419  0.69892473  0.04301075  0.21505376  0.05376344  0.27956989  0.05376344  0.58064516

#-log2转换
test.p.result <- -log2(p.result)
library(pheatmap)
opar <- par(no.readonly = TRUE)

#par(mai=c(2,0,2,0)); #调整画布上下左右大小
pdf("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/2.刻画克隆性和临床变量的关系/heatmap6.pdf",width=24,height=4)
pheatmap(test.p.result, 
		     display_numbers = matrix(ifelse(test.p.result > -log2(0.05), "*", ""), nrow(test.p.result)),	
		     cellwidth = 10, cellheight = 7,
		     treeheight_row=0, treeheight_col=0,
		     #border_color = NA,
		     color = colorRampPalette(c("RoyalBlue", "white","red"),bias=10)(20),fontsize=6)
dev.off()
par(opar)

###################logist回归判断克隆状态和转移(远端和淋巴结转移)以及疾病stage的关系###########################
## 选择响应变量为淋巴结转移进行多因素logistic回归分析
head(CRC.clonal.subclonal.data.list[[1]])
dim(CRC.clonal.subclonal.data.list[[1]])
## 对临床数据进行处理，N_stage为响应变量
logist.clinical.N.data <- lapply(CRC.clonal.subclonal.data.list,function(x){
                                        #筛选协变量进行后续logistic分析
                                        temp <- x[,c(2,3,4,5,6,7,9,17,18)]
                                        colnames(temp) <- c("sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site",colnames(x)[18])
                                        #处理AJCC_NODES_PATHOLOGIC_PN作为响应变量
                                        temp$N_stage[temp$N_stage == "N0"] <- 0 #没转移
                                        temp$N_stage[temp$N_stage == "N1" | temp$N_stage == "N2" | temp$N_stage == "N3"] <- 1#转移
                                        temp$N_stage <- as.numeric(temp$N_stage)
                                        #年龄依据中值离散化处理(没问题)
                                        # med <- median(temp$age);
                                        # tempage <- vector(mode="logical",length=length(temp$age));
                                        # tempage[which(temp$age >= med)] <- paste0(">=",med);
                                        # tempage[which(temp$age < med)] <- paste0("<",med);
                                        # temp$age <- tempage;
                                        # temp;
                                        #对协变量进行factor转换，控制levels的第一项为参考项,目前默认第一个level为参考项 
                                        temp$sex <- factor(temp$sex,levels=c("Female","Male"))
                                        temp$MSS_state[grep("MSI",temp$MSS_state)] <- "MSI"
                                        temp$MSS_state <- factor(temp$MSS_state,levels=c("MSI","MSS"))
                                        temp$T_stage[grep("T1|T2",temp$T_stage)] <- "T1/T2"
                                        temp$T_stage[grep("T3|T4",temp$T_stage)] <- "T3/T4"
                                        temp$T_stage <- factor(temp$T_stage,levels=c("T1/T2","T3/T4"))
                                        #temp$N_stage[grep("N1|N2",temp$N_stage)] <- "N1/N2"
                                        #temp$N_stage <- factor(temp$N_stage,levels=c("N0","N1/N2"))
                                        temp$M_stage <- factor(temp$M_stage,levels=c("M0","M1"))
                                        temp$AJCC_stage[grep("Stage 1|Stage 2",temp$AJCC_stage)] <- "Stage 1/Stage 2"
                                        temp$AJCC_stage[grep("Stage 3|Stage 4",temp$AJCC_stage)] <- "Stage 3/Stage 4"
                                        temp$AJCC_stage <- factor(temp$AJCC_stage,levels=c("Stage 1/Stage 2","Stage 3/Stage 4"))
                                        #temp$Residual[grep("R1|R2",temp$Residual)] <- "R1/R2"
                                        #temp$Residual <- factor(temp$Residual,levels=c("R0","R1/R2"))
                                        tempage <- vector(mode="logical",length=length(temp$age))
                                        tempage[which(temp$age >= 67.5)] <- paste0(">=",67.5)
                                        tempage[which(temp$age < 67.5)] <- paste0("<",67.5)
                                        temp$age <- tempage
                                        temp$age <- factor(temp$age,levels=c("<67.5",">=67.5"))
                                        #temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
                                        temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
                                        temp
                                    })
names(logist.clinical.N.data) <- names(CRC.clonal.subclonal.data.list)
##过滤掉没有变异的基因以及只有clone突变的基因进行后续分析，否则程序报错(因为突变状态之前已经处理为因子且设为3个等级)
filter.genes <- c("RBM10","BCOR","DCAF12L2","ZMYM3","AXIN1","POLE","SMAD3","NTHL1")
filter.genes <- c("RBM10", "BCOR", "DCAF12L2", "ZMYM3")
filter.genes.positions <- match(filter.genes,names(logist.clinical.N.data))
logist.clinical.N.data <- logist.clinical.N.data[-filter.genes.positions]
head(logist.clinical.N.data[[1]])
## 对临床数据进行logistic分析，以N_stage为响应变量
logist_N_result <- lapply(logist.clinical.N.data,function(x){
  print(colnames(x)[9])
  My_logistic(x,"N_stage",step=T)
})
sink(file="F:/logistic.stepwise.txt")
logist_N_result 
sink()
#结果:BCLAF1显著

logist_N_result_1 <- lapply(logist.clinical.N.data,function(x){
  print(colnames(x)[9])
  My_logistic(x,"N_stage",step=F)
})
sink(file="F:/logistic.txt")
logist_N_result_1 
sink()
#结果:BCLAF1显著
## 对临床数据进行处理，M_stage为响应变量
factor(CRC.clonal.subclonal.data.list[[1]]$AJCC_PATHOLOGIC_TUMOR_STAGE)
logist.clinical.M.data <- lapply(CRC.clonal.subclonal.data.list,function(x){
  #筛选协变量进行后续logistic分析
  temp <- x[,c(2,3,4,5,6,7,9,17,18)]
  colnames(temp) <- c("sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site",colnames(x)[18])
  #处理AJCC_NODES_PATHOLOGIC_PN作为响应变量
  temp$M_stage[temp$M_stage == "M0"] <- 0 #没转移
  temp$M_stage[temp$M_stage == "M1" | temp$N_stage == "M2" | temp$N_stage == "M3"] <- 1#转移
  temp$M_stage <- as.numeric(temp$M_stage)
  #年龄依据中值离散化处理(没问题)
  # med <- median(temp$age);
  # tempage <- vector(mode="logical",length=length(temp$age));
  # tempage[which(temp$age >= med)] <- paste0(">=",med);
  # tempage[which(temp$age < med)] <- paste0("<",med);
  # temp$age <- tempage;
  # temp;
  #对协变量进行factor转换，控制levels的第一项为参考项,目前默认第一个level为参考项 
  temp$sex <- factor(temp$sex,levels=c("Female","Male"))
  temp$MSS_state[grep("MSI",temp$MSS_state)] <- "MSI"
  temp$MSS_state <- factor(temp$MSS_state,levels=c("MSI","MSS"))
  temp$T_stage[grep("T1|T2",temp$T_stage)] <- "T1/T2"
  temp$T_stage[grep("T3|T4",temp$T_stage)] <- "T3/T4"
  temp$T_stage <- factor(temp$T_stage,levels=c("T1/T2","T3/T4"))
  temp$N_stage[grep("N1|N2",temp$N_stage)] <- "N1/N2"
  temp$N_stage <- factor(temp$N_stage,levels=c("N0","N1/N2"))
  #temp$M_stage <- factor(temp$M_stage,levels=c("M0","M1"))
  temp$AJCC_stage[grep("Stage 1|Stage 2",temp$AJCC_stage)] <- "Stage 1/Stage 2"
  temp$AJCC_stage[grep("Stage 3|Stage 4",temp$AJCC_stage)] <- "Stage 3/Stage 4"
  temp$AJCC_stage <- factor(temp$AJCC_stage,levels=c("Stage 1/Stage 2","Stage 3/Stage 4"))
  #temp$Residual[grep("R1|R2",temp$Residual)] <- "R1/R2"
  #temp$Residual <- factor(temp$Residual,levels=c("R0","R1/R2"))
  tempage <- vector(mode="logical",length=length(temp$age))
  tempage[which(temp$age >= 67.5)] <- paste0(">=",67.5)
  tempage[which(temp$age < 67.5)] <- paste0("<",67.5)
  temp$age <- tempage
  temp$age <- factor(temp$age,levels=c("<67.5",">=67.5"))
  #temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
  temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
  temp
})
names(logist.clinical.M.data) <- names(CRC.clonal.subclonal.data.list)
head(logist.clinical.M.data[[1]])
##过滤掉没有变异的基因以及只有clone突变的基因进行后续分析，否则程序报错(因为突变状态之前已经处理为因子且设为3个等级)
filter.genes <- c("RBM10", "BCOR", "DCAF12L2", "ZMYM3","AXIN1", "POLE", "SMAD3", "NTHL1")
filter.genes <- c("RBM10", "BCOR", "DCAF12L2", "ZMYM3")
filter.genes.positions <- match(filter.genes,names(logist.clinical.M.data))
logist.clinical.M.data <- logist.clinical.M.data[-filter.genes.positions]
head(logist.clinical.M.data[[1]])
dim(logist.clinical.M.data[[1]])
## 对临床数据进行logistic分析，以M_stage为响应变量
logist_M_result <- lapply(logist.clinical.M.data,function(x){
  print(colnames(x)[9])
  My_logistic(na.omit(x),"M_stage",step=T)
})
sink(file="F:/logistic_M_stepwise.txt")
logist_M_result
sink()

logist_M_result_1 <- lapply(logist.clinical.M.data,function(x){
  print(colnames(x)[9])
  My_logistic(x,"M_stage",step=F)
})
sink(file="F:/logistic_M.txt")
logist_M_result_1 
sink()
single.logistic(logist.clinical.M.data[[1]],"AJCC_METASTASIS_PATHOLOGIC_PM")
logist.clinical.M.data[[1]]$AJCC_PATHOLOGIC_TUMOR_STAGE

#########################评估驱动基因CCF值与临床指标的关系###############################
## 读入每个基因打上CCF值标签的临床数据
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.clinical.CCF.data.list.RData")

## 对临床数据表头进行清洗
head(CRC.clinical.CCF.data.list$APC)
logist.byDriverGene.CCF <- lapply(CRC.clinical.CCF.data.list,function(x){
  x <- x[,c(2,3,4,5,6,7,9,17,18)] #对临床数据进行筛选
  colnames(x) <- c("sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site",colnames(x)[9])
  tempage <- vector(mode="logical",length=length(x$age))
  tempage[which(x$age >= 67.5)] <- paste0(">=",67.5)
  tempage[which(x$age < 67.5)] <- paste0("<",67.5)
  x$age <- tempage
  x$age <- factor(x$age,levels=c("<67.5",">=67.5"))
  #对协变量进行factor转换，控制levels的第一项为参考项,目前默认第一个level为参考项 
  x$sex <- factor(x$sex,levels=c("Female","Male"))
  x$MSS_state[grep("MSI",x$MSS_state)] <- "MSI"
  x$MSS_state <- factor(x$MSS_state,levels=c("MSI","MSS"))
  x$T_stage[grep("T1|T2",x$T_stage)] <- "T1/T2"
  x$T_stage[grep("T3|T4",x$T_stage)] <- "T3/T4"
  x$T_stage <- factor(x$T_stage,levels=c("T1/T2","T3/T4"))
  x$N_stage[grep("N1|N2",x$N_stage)] <- "N1/N2"
  x$N_stage <- factor(x$N_stage,levels=c("N0","N1/N2"))
  x$M_stage <- factor(x$M_stage,levels=c("M0","M1"))
  x$AJCC_stage[grep("Stage 1|Stage 2",x$AJCC_stage)] <- "Stage 1/Stage 2"
  x$AJCC_stage[grep("Stage 3|Stage 4",x$AJCC_stage)] <- "Stage 3/Stage 4"
  x$AJCC_stage <- factor(x$AJCC_stage,levels=c("Stage 1/Stage 2","Stage 3/Stage 4"))
  #temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
  x$Tumour_site <- factor(x$Tumour_site,levels=c("right-sided","left-sided"))
  return(x) #如果不设置，默认返回最后一行
})
## 判断驱动基因CCF值和淋巴结转移的关系
head(logist.byDriverGene.CCF[[1]])
logist.N.byDriverGene.CCF <- lapply(logist.byDriverGene.CCF,function(x){
    x$N_stage <- ifelse(x$N_stage == "N0",0,1)
    x
})
head(logist.N.byDriverGene.CCF[[1]])
logist.N.byDriverGene.CCF.result <- lapply(logist.N.byDriverGene.CCF,function(x){
  My_logistic(na.omit(x),y="N_stage",step=F)
})
names(logist.N.byDriverGene.CCF.result) <- names(logist.N.byDriverGene.CCF)
sink(file="F:/CCF.logistic.N.txt")
logist.N.byDriverGene.CCF.result
sink()

## 判断驱动基因CCF值和远端转移的关系
logist.M.byDriverGene.CCF <- lapply(logist.byDriverGene.CCF,function(x){
  x$M_stage <- ifelse(x$M_stage == "M0",0,1)
  x
})
head(logist.M.byDriverGene.CCF[[1]])
logist.M.byDriverGene.CCF.result <- lapply(logist.M.byDriverGene.CCF,function(x){
  My_logistic(na.omit(x),y="M_stage",step=F)
})
names(logist.M.byDriverGene.CCF.result) <- names(logist.M.byDriverGene.CCF)
sink(file="F:/CCF.logistic.M.txt")
logist.M.byDriverGene.CCF.result
sink()






