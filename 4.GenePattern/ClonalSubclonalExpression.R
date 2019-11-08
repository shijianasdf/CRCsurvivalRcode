#'-----------------------------------------
#'        对表达数据进行处理和分析
#'-----------------------------------------
#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
#加载结直肠癌driver gene
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData")
#加载结直肠癌临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData")
#加载TCGA结肠癌和直肠癌表达数据
READ.ReadCount.matrix <- read.csv(file="D:/Rsources/Project/TCGA-READ/TCGAREAD-Counts.csv",stringsAsFactors = F)
COAD.ReadCount.matrix <- read.csv(file="D:/Rsources/Project/TCGA-COAD/TCGACOAD-Counts.csv",stringsAsFactors = F)
#加载ensemble gene注释信息
mart_export <- read.table("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/mart_export.txt",
                          sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
#加载直肠和结肠表达矩阵(表达数据行名不唯一，也就是基因名字相同的情况出现)
{
  READ.ReadCount.matrix[1:4,1:4]
  COAD.ReadCount.matrix[1:4,1:4]
  dim(READ.ReadCount.matrix) #56734   178
  dim(COAD.ReadCount.matrix) #56734   522
  intersect(substr(colnames(READ.ReadCount.matrix),1,12),substr(colnames(COAD.ReadCount.matrix),1,12))
  length(READ.ReadCount.matrix[,1])-length(unique(READ.ReadCount.matrix[,1])) #直肠癌有132个基因名字是重复的
  length(COAD.ReadCount.matrix[,1])-length(unique(COAD.ReadCount.matrix[,1])) #结肠癌有132个基因名字是重复的
  length(colnames(READ.ReadCount.matrix)) - length(unique(colnames(READ.ReadCount.matrix))) #列名字是唯一的
  length(colnames(COAD.ReadCount.matrix)) - length(unique(colnames(COAD.ReadCount.matrix))) #列名字是唯一的
  tail(sort(table(READ.ReadCount.matrix[,1])))
  table(sort(table(READ.ReadCount.matrix[,1])))
  #1     2     3     4     5     7 
  #56495    97     3     1     5     1
  tail(sort(table(COAD.ReadCount.matrix[,1])))
  table(sort(table(COAD.ReadCount.matrix[,1])))
  #1     2     3     4     5     7 
  #56495    97     3     1     5     1
}
{
  #过滤两个表达矩阵，重复基因名字随便选一个,算了选择第一个吧
  {
    test <- by(READ.ReadCount.matrix,READ.ReadCount.matrix[,1],function(x){x[1,]})
    class(test) # "by"
    test <- do.call(rbind.data.frame,test)
    dim(test) #56602   178
    test[1:6,1:6]
    READ.ReadCount.matrix <- test
    dim(READ.ReadCount.matrix)
    READ.ReadCount.matrix[1:6,1:6]
    
    test <- by(COAD.ReadCount.matrix,COAD.ReadCount.matrix[,1],function(x){x[1,]})
    class(test) # "by"
    test <- do.call(rbind.data.frame,test)
    dim(test) #56602   522
    test[1:6,1:6]
    COAD.ReadCount.matrix <- test
    dim(COAD.ReadCount.matrix) #56602   522
    COAD.ReadCount.matrix[1:6,1:6]
  }
  #合并两个表达矩阵
  {
    library(dplyr)
    COADREAD.ReadCount.matrix <- full_join(COAD.ReadCount.matrix,READ.ReadCount.matrix,by="X")
    dim(COADREAD.ReadCount.matrix) #56602 699
    COADREAD.ReadCount.matrix[1:8,1:4]
    rownames(COADREAD.ReadCount.matrix) <- COADREAD.ReadCount.matrix[,1]
    COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[,-1]
    COADREAD.ReadCount.matrix[1:8,1:4]
  }
  #过滤正常样本，留下肿瘤样本，并去除技术重复的样本
  {
    group_list <- ifelse(as.numeric(substr(colnames(COADREAD.ReadCount.matrix),14,15)) < 10,'tumor','normal') #样本标签
    table(group_list,useNA = "ifany")
    # group_list
    # normal  tumor 
    # 51    647
    pos <- group_list %in% "tumor"
    table(pos)
    COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[,pos]
    dim(COADREAD.ReadCount.matrix)
    #colnames(COADREAD.ReadCount.matrix) <- substr(colnames(COADREAD.ReadCount.matrix),1,12)
    COADREAD.ReadCount.matrix[1:8,1:4]
    duplicated(substr(colnames(COADREAD.ReadCount.matrix),1,12)) #病人会有技术重复
    COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[,!duplicated(substr(colnames(COADREAD.ReadCount.matrix),1,12))] #去除技术重复的样本
    colnames(COADREAD.ReadCount.matrix) <- gsub("\\.","-",colnames(COADREAD.ReadCount.matrix))
  }
  #筛选表达矩阵的样本,和我们克隆突变数据以及临床数据一致
  {
    intersect.patients <- intersect(colorectal.clinic.data$PATIENT_ID,colorectal.select.table$patient)
    length(intersect.patients) #536个交集病人
    length(unique(intersect.patients)) #536个交集病人
    pos <- substr(colnames(COADREAD.ReadCount.matrix),1,12) %in% intersect.patients
    intersect(substr(colnames(COADREAD.ReadCount.matrix),1,12),intersect.patients) #532个病人有交集,奇怪为什么只有532个病人的表达数据
    COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[,pos]
    dim(COADREAD.ReadCount.matrix) #56602   532
    head(COADREAD.ReadCount.matrix)
  }
  #ENSG id 转换,并保存表达矩阵
  {
    #下载ensemble gene注释信息
    mart_export <- read.table("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/mart_export.txt",
                              sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
    head(mart_export)
    tail(sort(table(mart_export$Gene.stable.ID))) #ENSEM无重复
    tail(sort(table(mart_export$Gene.name))) #基因symbol有重复，多个ENSEM对应同一个gene symbol
    table(sort(table(mart_export$Gene.name))) 
    pos <- match(rownames(COADREAD.ReadCount.matrix),mart_export$Gene.stable.ID) #默认match到第一个，后面有就不要了
    which(is.na(pos)) #无NA
    tail(mart_export[pos,])
    mart_export <- mart_export[pos,]
    sort(table(mart_export$Gene.type))
    pos1 <- which(mart_export$Gene.type == "protein_coding")
    mart_export <- mart_export[pos1,]
    COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[pos1,]
    COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[!duplicated(mart_export$Gene.name),] #去除重复的基因
    rownames(COADREAD.ReadCount.matrix) <- mart_export$Gene.name[!duplicated(mart_export$Gene.name)]
    dim(COADREAD.ReadCount.matrix) #19634   532
    COADREAD.ReadCount.matrix[1:4,1:4]
    #储存表达矩阵
    save(COADREAD.ReadCount.matrix,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/colorectal.expression.data.RData")
  }
  
}
#处理TCGA表达数据
temp.COADREAD.ReadCount.matrix <- dataProcessing(COAD.ReadCount.matrix,READ.ReadCount.matrix,mart_export)
save(temp.COADREAD.ReadCount.matrix,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/colorectal.expression.ENSG.data.RData")
##突变相对于正常样本算差异
{
  #导入结直肠癌依照克隆状态分类的临床数据
  load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.survival.data.list.RData")
  #导入结直肠癌表达数据
  #load(file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/colorectal.expression.data.RData")
  load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/colorectal.expression.ENSG.data.RData")
  sig <- c("CASP8","SMAD2","ANK1")
  pos <- match(sig,names(CRC.survival.data.list))
  CRC.clonalSubclonal.data.list <- CRC.survival.data.list[pos]
  
  #生成计算差异所需要的表达矩阵和样本标签
  source("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/CloneSurvivalCode/0.Codelibrary/DiffExpression.R")
  source("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/CloneSurvivalCode/0.Codelibrary/ExpressionScore.R")
  source("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/CloneSurvivalCode/0.Codelibrary/ScientificPlot.R")
  #列表里面嵌入列表,每个列表对应一个基因的表达矩阵,样本基于突变状态的分组信息,以及对应表达矩阵样本的克隆状态信息
  tempList <- list()
  for(i in names(CRC.clonalSubclonal.data.list)){
    tempList[[i]] <- makeDiffData(temp.COADREAD.ReadCount.matrix,CRC.clonalSubclonal.data.list[[i]])
  }
  save(tempList,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/exprSet.diffLabel.clinical.RData")
  #储存过滤后的最终表达矩阵，
  #save(temp.COADREAD.ReadCount.matrix,temp.clinical,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP7_DEG.RData")
  #三种算法差异分析
  {
    DESeq2_DEG <- DESeq2.function(tempList[[1]]$exprSet,tempList[[1]]$group_list,"MUT","WT")
    EdgeR_DEG <- EdgeR.function(tempList[[1]]$exprSet,tempList[[1]]$group_list)
    limma_DEG <- limma.voom(tempList[[1]]$exprSet,tempList[[1]]$group_list,"MUT-WT")
    DESeq2_DEG <- IDmap(DESeq2_DEG,mart_export)
    EdgeR_DEG <- IDmap(EdgeR_DEG,mart_export)
    limma_DEG <- IDmap(limma_DEG,mart_export)
    save(DESeq2_DEG,EdgeR_DEG,limma_DEG,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP8_DEG.RData")
    
    DESeq2_DEG_SMAD2 <- DESeq2.function(tempList[[2]]$exprSet,tempList[[2]]$group_list,"MUT","WT")
    EdgeR_DEG_SMAD2 <- EdgeR.function(tempList[[2]]$exprSet,tempList[[2]]$group_list)
    limma_DEG_SMAD2 <- limma.voom(tempList[[2]]$exprSet,tempList[[2]]$group_list,"MUT-WT")
    DESeq2_DEG_SMAD2 <- IDmap(DESeq2_DEG_SMAD2,mart_export)
    EdgeR_DEG_SMAD2 <- IDmap(EdgeR_DEG_SMAD2,mart_export)
    limma_DEG_SMAD2 <- IDmap(limma_DEG_SMAD2,mart_export)
    save(DESeq2_DEG_SMAD2,EdgeR_DEG_SMAD2,limma_DEG_SMAD2,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/SMAD2_DEG.RData")
    
    
    DESeq2_DEG_ANK1 <- DESeq2.function(tempList[[3]]$exprSet,tempList[[3]]$group_list,"MUT","WT")
    EdgeR_DEG_ANK1 <- EdgeR.function(tempList[[3]]$exprSet,tempList[[3]]$group_list)
    limma_DEG_ANK1 <- limma.voom(tempList[[3]]$exprSet,tempList[[3]]$group_list,"MUT-WT")
    DESeq2_DEG_ANK1 <- IDmap(DESeq2_DEG_ANK1,mart_export)
    EdgeR_DEG_ANK1 <- IDmap(EdgeR_DEG_ANK1,mart_export)
    limma_DEG_ANK1 <- IDmap(limma_DEG_ANK1,mart_export)
    save(DESeq2_DEG_ANK1,EdgeR_DEG_ANK1,limma_DEG_ANK1,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ANK1_DEG.RData")
    
  }
  
}
#火山图
{
  DESeq2_DEG$lab <- ifelse(abs(DESeq2_DEG$log2FoldChange) > 1 & DESeq2_DEG$padj < 0.05,ifelse(DESeq2_DEG$log2FoldChange > 1,"UP","DOWN"),"NO")
  volcano.Deseq2.plot(DESeq2_DEG,"CASP8","F:/Rsources/Project/CASP8.Deseq2.volcano.pdf")
  EdgeR_DEG$lab <- ifelse(abs(EdgeR_DEG$logFC) > 1 & EdgeR_DEG$FDR < 0.05,ifelse(EdgeR_DEG$logFC > 1,"UP","DOWN"),"NO")
  volcano.EdgeR.plot(EdgeR_DEG,"CASP8","F:/Rsources/Project/CASP8.EdgeR.volcano.pdf") #预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/
  limma_DEG$lab <- ifelse(abs(limma_DEG$logFC) > 1 & limma_DEG$adj.P.Val < 0.05,ifelse(limma_DEG$logFC > 1,"UP","DOWN"),"NO")
  volcano.limma.voom.plot(limma_DEG,"CASP8","F:/Rsources/Project/CASP8.limma.volcano.pdf") #预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/
  save(DESeq2_DEG,EdgeR_DEG,limma_DEG,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP7_DEG.RData")
  
  DESeq2_DEG_SMAD2$lab <- ifelse(abs(DESeq2_DEG_SMAD2$log2FoldChange) > 1 & DESeq2_DEG_SMAD2$padj < 0.05,ifelse(DESeq2_DEG_SMAD2$log2FoldChange > 1,"UP","DOWN"),"NO")
  volcano.Deseq2.plot(DESeq2_DEG_SMAD2,"SMAD2","F:/Rsources/Project/SMAD2.Deseq2.volcano.pdf")
  EdgeR_DEG_SMAD2$lab <- ifelse(abs(EdgeR_DEG_SMAD2$logFC) > 1 & EdgeR_DEG_SMAD2$FDR < 0.05,ifelse(EdgeR_DEG_SMAD2$logFC > 1,"UP","DOWN"),"NO")
  volcano.EdgeR.plot(EdgeR_DEG_SMAD2,"SMAD2","F:/Rsources/Project/SMAD2.EdgeR.volcano.pdf") #预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/
  limma_DEG_SMAD2$lab <- ifelse(abs(limma_DEG_SMAD2$logFC) > 1 & limma_DEG_SMAD2$adj.P.Val < 0.05,ifelse(limma_DEG_SMAD2$logFC > 1,"UP","DOWN"),"NO")
  volcano.limma.voom.plot(limma_DEG_SMAD2,"SMAD2","F:/Rsources/Project/SMAD2.limma.volcano.pdf") #预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/
  save(DESeq2_DEG_SMAD2,EdgeR_DEG_SMAD2,limma_DEG_SMAD2,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/SMAD2_DEG.RData")
  
  DESeq2_DEG_ANK1$lab <- ifelse(abs(DESeq2_DEG_ANK1$log2FoldChange) > 1 & DESeq2_DEG_ANK1$padj < 0.05,ifelse(DESeq2_DEG_ANK1$log2FoldChange > 1,"UP","DOWN"),"NO")
  volcano.Deseq2.plot(DESeq2_DEG_ANK1,"ANK1","F:/Rsources/Project/ANK1.Deseq2.volcano.pdf")
  EdgeR_DEG_ANK1$lab <- ifelse(abs(EdgeR_DEG_ANK1$logFC) > 1 & EdgeR_DEG_ANK1$FDR < 0.05,ifelse(EdgeR_DEG_ANK1$logFC > 1,"UP","DOWN"),"NO")
  volcano.EdgeR.plot(EdgeR_DEG_ANK1,"ANK1","F:/Rsources/Project/ANK1.EdgeR.volcano.pdf") #预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/
  limma_DEG_ANK1$lab <- ifelse(abs(limma_DEG_ANK1$logFC) > 1 & limma_DEG_ANK1$adj.P.Val < 0.05,ifelse(limma_DEG_ANK1$logFC > 1,"UP","DOWN"),"NO")
  volcano.limma.voom.plot(limma_DEG_ANK1,"ANK1","F:/Rsources/Project/ANK1.limma.volcano.pdf") #预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/
  save(DESeq2_DEG_ANK1,EdgeR_DEG_ANK1,limma_DEG_SMAD2,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ANK1_DEG.RData")

}
#上调韦恩图
{
  GroupA <- rownames(DESeq2_DEG)[DESeq2_DEG$lab == "UP"]
  GroupB <- rownames(EdgeR_DEG)[EdgeR_DEG$lab == "UP"]
  GroupC <- rownames(limma_DEG)[limma_DEG$lab == "UP"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/CASP8.up.venn.tiff")
  
  GroupA <- rownames(DESeq2_DEG_SMAD2)[DESeq2_DEG_SMAD2$lab == "UP"]
  GroupB <- rownames(EdgeR_DEG_SMAD2)[EdgeR_DEG_SMAD2$lab == "UP"]
  GroupC <- rownames(limma_DEG_SMAD2)[limma_DEG_SMAD2$lab == "UP"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/SMAD2.up.venn.tiff")
  
  GroupA <- rownames(DESeq2_DEG_ANK1)[DESeq2_DEG_ANK1$lab == "UP"]
  GroupB <- rownames(EdgeR_DEG_ANK1)[EdgeR_DEG_ANK1$lab == "UP"]
  GroupC <- rownames(limma_DEG_ANK1)[limma_DEG_ANK1$lab == "UP"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/ANK1.up.venn.tiff")
  
  library(VennDiagram)
  library(grDevices)
  GroupA <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "UP"]
  GroupB <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "UP"]
  input <- list(GroupA,GroupB)
  temp <- venn.diagram(x=input,
               fill=c("skyblue", "pink1"),
               alpha=c(0.5,0.5),cex=2,
               cat.fontface=4,
               category = c("Deseq2", "EdgeR"),filename = NULL)
  pdf(file="D:/ANK1.up.venn.pdf")
  grid.draw(temp)
  dev.off()
  
  
  
  library(VennDiagram)
  library(grDevices)
  GroupA <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "DOWN"]
  GroupB <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "DOWN"]
  input <- list(GroupA,GroupB)
  temp <- venn.diagram(x=input,
                       fill=c("skyblue", "pink1"),
                       alpha=c(0.5,0.5),cex=2,
                       cat.fontface=4,
                       category = c("Deseq2", "EdgeR"),filename = NULL)
  pdf(file="D:/ANK1.down.venn.pdf")
  grid.draw(temp)
  dev.off()
  
}
{
  GroupA <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "UP"]
  GroupB <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "UP"]
  GroupC <- limma_DEG$ENSG[limma_DEG$lab == "UP"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/CASP8.up.venn.tiff")
  
  GroupA <- DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab == "UP"]
  GroupB <- EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab == "UP"]
  GroupC <- limma_DEG_SMAD2$ENSG[limma_DEG_SMAD2$lab == "UP"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/SMAD2.up.venn.tiff")
  
  GroupA <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "UP"]
  GroupB <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "UP"]
  GroupC <- limma_DEG_ANK1$ENSG[limma_DEG_ANK1$lab == "UP"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/ANK1.up.venn.tiff")
}
#下调韦恩图
{
  GroupA1 <- rownames(DESeq2_DEG)[DESeq2_DEG$lab == "DOWN"]
  GroupB1 <- rownames(EdgeR_DEG)[EdgeR_DEG$lab == "DOWN"]
  GroupC1 <- rownames(limma_DEG)[limma_DEG$lab == "DOWN"]
  input1 <- list(GroupA1,GroupB1,GroupC1)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/CASP8.down.venn.tiff")
  
  GroupA <- rownames(DESeq2_DEG_SMAD2)[DESeq2_DEG_SMAD2$lab == "DOWN"]
  GroupB <- rownames(EdgeR_DEG_SMAD2)[EdgeR_DEG_SMAD2$lab == "DOWN"]
  GroupC <- rownames(limma_DEG_SMAD2)[limma_DEG_SMAD2$lab == "DOWN"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/SMAD2.down.venn.tiff")
  
  GroupA <- rownames(DESeq2_DEG_ANK1)[DESeq2_DEG_ANK1$lab == "DOWN"]
  GroupB <- rownames(EdgeR_DEG_ANK1)[EdgeR_DEG_ANK1$lab == "DOWN"]
  GroupC <- rownames(limma_DEG_ANK1)[limma_DEG_ANK1$lab == "DOWN"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/ANK1.down.venn.tiff")
}
{
  GroupA1 <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "DOWN"]
  GroupB1 <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "DOWN"]
  GroupC1 <- limma_DEG$ENSG[limma_DEG$lab == "DOWN"]
  input1 <- list(GroupA1,GroupB1,GroupC1)
  Venn.plot(input1,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/CASP8.down.venn.tiff")
  
  GroupA <- DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab == "DOWN"]
  GroupB <- EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab == "DOWN"]
  GroupC <- limma_DEG_SMAD2$ENSG[limma_DEG_SMAD2$lab == "DOWN"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/SMAD2.down.venn.tiff")
  
  GroupA <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "DOWN"]
  GroupB <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "DOWN"]
  GroupC <- limma_DEG_ANK1$ENSG[limma_DEG_ANK1$lab == "DOWN"]
  input <- list(GroupA,GroupB,GroupC)
  Venn.plot(input,fill=c("skyblue", "pink1", "mediumorchid"),
            category = c("Deseq2", "EdgeR", "Limma-Voom"),
            filename="F:/ANK1.down.venn.tiff")
}
#比较不同克隆状态下差异基因的表达活性的比较
{
  table(tempList[[1]]$clinical.data$sample.label)
  #X <- log2(tempList[[1]]$exprSet[,!grepl("WT",tempList[[1]]$clinical.data$sample.label)]+1) #read count标化
  #去除正常样本
  X <- tempList[[1]]$exprSet[,!grepl("WT",tempList[[1]]$clinical.data$sample.label)] 
  label <- tempList[[1]]$clinical.data$sample.label[!grepl("WT",tempList[[1]]$clinical.data$sample.label)]
  GroupA <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "UP"]
  GroupB <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "UP"]
  GroupC <- limma_DEG$ENSG[limma_DEG$lab == "UP"]
  GroupA1 <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "DOWN"]
  GroupB1 <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "DOWN"]
  GroupC1 <- limma_DEG$ENSG[limma_DEG$lab == "DOWN"]
  UP_genes <- intersect(GroupA,GroupB)  #intersect(intersect(GroupA,GroupB),GroupC)  #上调基因
  Down_genes <- intersect(GroupA1,GroupB1) #intersect(intersect(GroupA1,GroupB1),GroupC1)  #下调基因
  genesets <- list(up=UP_genes,down=Down_genes)
  
  #计算表达活性得分
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson") 
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.up.plot,"D:/Rsources/Project/tempReuslt/Casp8UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.down.plot,"D:/Rsources/Project/tempReuslt/Casp8DOWNgenes.score.pdf")

  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.up.plot,"D:/Rsources/Project/tempReuslt/Casp8UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.down.plot,"D:/Rsources/Project/tempReuslt/Casp8DOWNgenes.ssgsea.score.pdf")
  
  X <- tempList[[2]]$exprSet[,!grepl("WT",tempList[[2]]$clinical.data$sample.label)] 
  label <- tempList[[2]]$clinical.data$sample.label[!grepl("WT",tempList[[2]]$clinical.data$sample.label)]
  GroupA <- rownames(DESeq2_DEG_SMAD2)[DESeq2_DEG_SMAD2$lab == "UP"]
  GroupB <- rownames(EdgeR_DEG_SMAD2)[EdgeR_DEG_SMAD2$lab == "UP"]
  GroupC <- rownames(limma_DEG_SMAD2)[limma_DEG_SMAD2$lab == "UP"]
  GroupA1 <- rownames(DESeq2_DEG_SMAD2)[DESeq2_DEG_SMAD2$lab == "DOWN"]
  GroupB1 <- rownames(EdgeR_DEG_SMAD2)[EdgeR_DEG_SMAD2$lab == "DOWN"]
  GroupC1 <- rownames(limma_DEG_SMAD2)[limma_DEG_SMAD2$lab == "DOWN"]
  UP_genes <- intersect(intersect(GroupA,GroupB),GroupC)  #上调基因
  Down_genes <- intersect(intersect(GroupA1,GroupB1),GroupC1)  #下调基因
  genesets <- list(up=UP_genes,down=Down_genes)
  #计算表达活性得分
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson") 
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/SMAD2UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/SMAD2DOWNgenes.score.pdf")
  
  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/SMAD2UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/SMAD2DOWNgenes.ssgsea.score.pdf")
  
  
  X <- tempList[[3]]$exprSet[,!grepl("WT",tempList[[3]]$clinical.data$sample.label)] #read count标化
  label <- tempList[[3]]$clinical.data$sample.label[!grepl("WT",tempList[[3]]$clinical.data$sample.label)]
  GroupA <- rownames(DESeq2_DEG_ANK1)[DESeq2_DEG_ANK1$lab == "UP"]
  GroupB <- rownames(EdgeR_DEG_ANK1)[EdgeR_DEG_ANK1$lab == "UP"]
  GroupC <- rownames(limma_DEG_ANK1)[limma_DEG_ANK1$lab == "UP"]
  GroupA1 <- rownames(DESeq2_DEG_ANK1)[DESeq2_DEG_ANK1$lab == "DOWN"]
  GroupB1 <- rownames(EdgeR_DEG_ANK1)[EdgeR_DEG_ANK1$lab == "DOWN"]
  GroupC1 <- rownames(limma_DEG_ANK1)[limma_DEG_ANK1$lab == "DOWN"]
  UP_genes <- intersect(intersect(GroupA,GroupB),GroupC)  #上调基因
  Down_genes <- intersect(intersect(GroupA1,GroupB1),GroupC1)  #下调基因
  genesets <- list(up=UP_genes,down=Down_genes)
  #计算表达活性得分
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson") 
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/ANK1UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/ANK1DOWNgenes.score.pdf")
  
  library(reshape2)
  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/ANK1UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/ANK1DOWNgenes.ssgsea.score.pdf")
  
  pheatmap::pheatmap(es.max1, scale = "row",
                     clustering_method = "complete",
                     annotation_col = annotation_df[,c(4,5)],
                     clustering_distance_rows = "correlation",
                     show_rownames = F,
                     file = paste0("Results/2015-11-08-GSVA-biocart-primary.pdf"),
                     width = 10, height = 8)
  
    
}
{
  table(tempList[[1]]$clinical.data$sample.label)
  #X <- log2(tempList[[1]]$exprSet[,!grepl("WT",tempList[[1]]$clinical.data$sample.label)]+1) #read count标化
  X <- tempList[[1]]$exprSet[,!grepl("WT",tempList[[1]]$clinical.data$sample.label)]
  label <- tempList[[1]]$clinical.data$sample.label[!grepl("WT",tempList[[1]]$clinical.data$sample.label)]
  GroupA <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "UP"]
  GroupB <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "UP"]
  GroupC <- limma_DEG$ENSG[limma_DEG$lab == "UP"]
  GroupA1 <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "DOWN"]
  GroupB1 <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "DOWN"]
  GroupC1 <- limma_DEG$ENSG[limma_DEG$lab == "DOWN"]
  #UP_genes <- intersect(intersect(GroupA,GroupB),GroupC)  #上调基因
  #Down_genes <- intersect(intersect(GroupA1,GroupB1),GroupC1)  #下调基因
  UP_genes <- intersect(GroupA,GroupB)  #上调基因
  Down_genes <- intersect(GroupA1,GroupB1)  #下调基因
  genesets <- list(up=UP_genes,down=Down_genes)
  
  #计算表达活性得分
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson") 
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/Casp8UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/Casp8DOWNgenes.score.pdf")
  
  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/Casp8UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/Casp8DOWNgenes.ssgsea.score.pdf")
  
  X <- tempList[[2]]$exprSet[,!grepl("WT",tempList[[2]]$clinical.data$sample.label)] 
  label <- tempList[[2]]$clinical.data$sample.label[!grepl("WT",tempList[[2]]$clinical.data$sample.label)]
  GroupA <- DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab == "UP"]
  GroupB <- EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab == "UP"]
  GroupC <- limma_DEG_SMAD2$ENSG[limma_DEG_SMAD2$lab == "UP"]
  GroupA1 <- DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab == "DOWN"]
  GroupB1 <- EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab == "DOWN"]
  GroupC1 <- limma_DEG_SMAD2$ENSG[limma_DEG_SMAD2$lab == "DOWN"]
  #UP_genes <- intersect(intersect(GroupA,GroupB),GroupC)  #上调基因
  #Down_genes <- intersect(intersect(GroupA1,GroupB1),GroupC1)  #下调基因
  UP_genes <- intersect(GroupA,GroupB)  #上调基因
  Down_genes <- intersect(GroupA1,GroupB1)  #下调基因
  genesets <- list(up=UP_genes,down=Down_genes)
  #计算表达活性得分
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson") 
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/SMAD2UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/SMAD2DOWNgenes.score.pdf")
  
  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/SMAD2UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/SMAD2DOWNgenes.ssgsea.score.pdf")
  
  
  X <- tempList[[3]]$exprSet[,!grepl("WT",tempList[[3]]$clinical.data$sample.label)] #read count标化
  label <- tempList[[3]]$clinical.data$sample.label[!grepl("WT",tempList[[3]]$clinical.data$sample.label)]
  GroupA <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "UP"]
  GroupB <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "UP"]
  GroupC <- limma_DEG_ANK1$ENSG[limma_DEG_ANK1$lab == "UP"]
  GroupA1 <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "DOWN"]
  GroupB1 <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "DOWN"]
  GroupC1 <- limma_DEG_ANK1$ENSG[limma_DEG_ANK1$lab == "DOWN"]
  #UP_genes <- intersect(intersect(GroupA,GroupB),GroupC)  #上调基因
  #Down_genes <- intersect(intersect(GroupA1,GroupB1),GroupC1)  #下调基因
  UP_genes <- intersect(GroupA,GroupB)  #上调基因
  Down_genes <- intersect(GroupA1,GroupB1)  #下调基因
  genesets <- list(up=UP_genes,down=Down_genes)
  #计算表达活性得分
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson") 
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/ANK1UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/ANK1DOWNgenes.score.pdf")
  
  library(reshape2)
  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.up.plot,"F:/Rsources/Project/tempReuslt/ANK1UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.down.plot,"F:/Rsources/Project/tempReuslt/ANK1DOWNgenes.ssgsea.score.pdf")

}




