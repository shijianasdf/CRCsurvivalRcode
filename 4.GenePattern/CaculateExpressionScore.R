#-------------------------------------------------------------------------
#       计算基因集合表达活性得分(WT，clonal，subclonal)
#-------------------------------------------------------------------------
##导入表达矩阵
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/colorectal.expression.ENSG.data.RData")
##加载差异表达基因
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP8_DEG.RData")
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/SMAD2_DEG.RData")
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ANK1_DEG.RData")
#加载结直肠癌依照克隆状态分类的临床数据
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.survival.data.list.RData")

##提取过滤后的表达矩阵
sig <- c("CASP8","SMAD2","ANK1")
pos <- match(sig,names(CRC.survival.data.list))
CRC.clonalSubclonal.data.list <- CRC.survival.data.list[pos]
tempList <- list()
for(i in names(CRC.clonalSubclonal.data.list)){
  tempList[[i]] <- makeDiffData(temp.COADREAD.ReadCount.matrix,CRC.clonalSubclonal.data.list[[i]])
}
save(tempList,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/exprSet.diffLabel.clinical.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/exprSet.diffLabel.clinical.RData")

{
  ##筛选CASP8差异表达基因分别为DESeq2和EdgeR top20,top30,top50等等
  DESeq2_DEG$lab <- ifelse(abs(DESeq2_DEG$log2FoldChange) > 1 & DESeq2_DEG$padj < 0.05,ifelse(DESeq2_DEG$log2FoldChange > 1,"UP","DOWN"),"NO")
  EdgeR_DEG$lab <- ifelse(abs(EdgeR_DEG$logFC) > 1 & EdgeR_DEG$FDR < 0.05,ifelse(EdgeR_DEG$logFC > 1,"UP","DOWN"),"NO")
  GroupA <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "UP"]
  GroupB <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "UP"]
  Up_genes <- intersect(GroupA,GroupB) #318个上调基因
  pos1 <- match(Up_genes,DESeq2_DEG$ENSG[DESeq2_DEG$lab=="UP"])
  pos2 <- match(Up_genes,EdgeR_DEG$ENSG[EdgeR_DEG$lab=="UP"])
  Up_genes <- Up_genes[which((pos1 <= 10) & (pos2 <= 10))] #筛选在DESeq2和EdgeR中都显著差异并且排秩在上调中都靠前的基因
  
  GroupA1 <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "DOWN"]
  GroupB1 <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "DOWN"]
  Down_genes <- intersect(GroupA1,GroupB1)  #下调基因
  pos1 <- match(Down_genes,DESeq2_DEG$ENSG[DESeq2_DEG$lab=="DOWN"])
  pos2 <- match(Down_genes,EdgeR_DEG$ENSG[EdgeR_DEG$lab=="DOWN"])
  Down_genes <- Down_genes[which((pos1 <= 10) & (pos2 <= 10))] #筛选在DESeq2和EdgeR中都显著差异并且排秩在下调中都靠前的基因
  
  #计算表达活性得分
  X <- tempList[[1]]$exprSet
  #去除正常样本
  X <- tempList[[1]]$exprSet[,!grepl("WT",tempList[[1]]$clinical.data$sample.label)]
  label <- tempList[[1]]$clinical.data$sample.label[!grepl("WT",tempList[[1]]$clinical.data$sample.label)]
  
  genesets <- list(up=Up_genes,down=Down_genes)
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson")
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  colnames(es.max) <- as.character(tempList[[1]]$clinical.data$sample.label)
  colnames(es.max1) <- as.character(tempList[[1]]$clinical.data$sample.label)

  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.up.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/Casp8UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.down.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/Casp8DOWNgenes.score.pdf")
  
  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssGSEA
  box.jitter.ggpubr.plot(dat.up.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/Casp8UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssGSEA
  box.jitter.ggpubr.plot(dat.down.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/Casp8DOWNgenes.ssgsea.score.pdf")
}
{
  ##筛选SMAD2差异表达基因分别为DESeq2和EdgeR top20,top30,top50等等
  DESeq2_DEG_SMAD2$lab <- ifelse(abs(DESeq2_DEG_SMAD2$log2FoldChange) > 1 & DESeq2_DEG_SMAD2$padj < 0.05,ifelse(DESeq2_DEG_SMAD2$log2FoldChange > 1,"UP","DOWN"),"NO")
  EdgeR_DEG_SMAD2$lab <- ifelse(abs(EdgeR_DEG_SMAD2$logFC) > 1 & EdgeR_DEG_SMAD2$FDR < 0.05,ifelse(EdgeR_DEG_SMAD2$logFC > 1,"UP","DOWN"),"NO")
  GroupA <- DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab == "UP"]
  GroupB <- EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab == "UP"]
  Up_genes <- intersect(GroupA,GroupB) 
  pos1 <- match(Up_genes,DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab=="UP"])
  pos2 <- match(Up_genes,EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab=="UP"])
  Up_genes <- Up_genes[which((pos1 <= 30) & (pos2 <= 30))] #筛选在DESeq2和EdgeR中都显著差异并且排秩在上调中都靠前的基因
  
  GroupA1 <- DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab == "DOWN"]
  GroupB1 <- EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab == "DOWN"]
  Down_genes <- intersect(GroupA1,GroupB1)  #下调基因
  pos1 <- match(Down_genes,DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab=="DOWN"])
  pos2 <- match(Down_genes,EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab=="DOWN"])
  Down_genes <- Down_genes[which((pos1 <= 30) & (pos2 <= 30))] #筛选在DESeq2和EdgeR中都显著差异并且排秩在下调中都靠前的基因
  
  #计算表达活性得分
  X <- tempList[[2]]$exprSet
  #去除正常样本
  X <- tempList[[2]]$exprSet[,!grepl("WT",tempList[[2]]$clinical.data$sample.label)] 
  label <- tempList[[2]]$clinical.data$sample.label[!grepl("WT",tempList[[2]]$clinical.data$sample.label)]
  genesets <- list(up=Up_genes,down=Down_genes)
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson")
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  
  colnames(es.max) <- as.character(tempList[[2]]$clinical.data$sample.label)
  colnames(es.max1) <- as.character(tempList[[2]]$clinical.data$sample.label)
  
  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.up.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/SMAD2UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.down.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/SMAD2DOWNgenes.score.pdf")
  
  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.up.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/SMAD2UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.down.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/SMAD2DOWNgenes.ssgsea.score.pdf")
}

{
  ##筛选ANK1差异表达基因分别为DESeq2和EdgeR top20,top30,top50等等
  DESeq2_DEG_ANK1$lab <- ifelse(abs(DESeq2_DEG_ANK1$log2FoldChange) > 1 & DESeq2_DEG_ANK1$padj < 0.05,ifelse(DESeq2_DEG_ANK1$log2FoldChange > 1,"UP","DOWN"),"NO")
  EdgeR_DEG_ANK1$lab <- ifelse(abs(EdgeR_DEG_ANK1$logFC) > 1 & EdgeR_DEG_ANK1$FDR < 0.05,ifelse(EdgeR_DEG_ANK1$logFC > 1,"UP","DOWN"),"NO")
  GroupA <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "UP"]
  GroupB <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "UP"]
  Up_genes <- intersect(GroupA,GroupB) #318个上调基因
  pos1 <- match(Up_genes,DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab=="UP"])
  pos2 <- match(Up_genes,EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab=="UP"])
  Up_genes <- Up_genes[which((pos1 <= 50) & (pos2 <= 50))] #筛选在DESeq2和EdgeR中都显著差异并且排秩在上调中都靠前的基因
  
  GroupA1 <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "DOWN"]
  GroupB1 <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "DOWN"]
  Down_genes <- intersect(GroupA1,GroupB1)  #下调基因
  pos1 <- match(Down_genes,DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab=="DOWN"])
  pos2 <- match(Down_genes,EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab=="DOWN"])
  Down_genes <- Down_genes[which((pos1 <= 50) & (pos2 <= 50))] #筛选在DESeq2和EdgeR中都显著差异并且排秩在下调中都靠前的基因
  
  #计算表达活性得分
  X <- tempList[[3]]$exprSet
  #去除正常样本
  X <- tempList[[3]]$exprSet[,!grepl("WT",tempList[[3]]$clinical.data$sample.label)] #read count标化
  label <- tempList[[3]]$clinical.data$sample.label[!grepl("WT",tempList[[3]]$clinical.data$sample.label)]
  genesets <- list(up=Up_genes,down=Down_genes)
  es.max <- CaculateExprSCore(as.matrix(X),genesets,"gsva","Poisson") 
  es.max1 <- CaculateExprSCore(as.matrix(X),genesets,"ssgsea","Poisson")
  colnames(es.max) <- as.character(label)
  colnames(es.max1) <- as.character(label)
  
  colnames(es.max) <- as.character(tempList[[3]]$clinical.data$sample.label)
  colnames(es.max1) <- as.character(tempList[[3]]$clinical.data$sample.label)
  
  library(reshape2)
  dat.plot <- melt(es.max)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.up.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/ANK1UPgenes.score.pdf")
  #下调基因不同克隆状态下的表达活性差异gsva
  box.jitter.ggpubr.plot(dat.down.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/ANK1DOWNgenes.score.pdf")
  
  dat.plot <- melt(es.max1)
  dat.up.plot <- dat.plot[dat.plot$Var1=="up",]
  dat.down.plot <- dat.plot[dat.plot$Var1=="down",]
  #上调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.up.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/ANK1UPgenes.ssgsea.score.pdf")
  #下调基因不同克隆状态下的表达活性差异ssgsea
  box.jitter.ggpubr.plot(dat.down.plot,"D:/Rsources/Project/tempReuslt/SubclonalVSclonalVSwildtype/ANK1DOWNgenes.ssgsea.score.pdf")
}

##导入富集的GO和KEGG结果
files <- list.files(path="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/GoKegg",full.names = T)
files.csv <- files[grep(pattern=".csv",files)]
files.names <- substring(basename(files.csv),1,nchar(basename(files.csv))-4)
GoKeggEnrichList <- list()
for(i in 1:length(files.csv)){
  GoKeggEnrichList[[i]] <- read.csv(files.csv[i])
  GoKeggEnrichList[[i]] <- GoKeggEnrichList[[i]][,-1]
}
names(GoKeggEnrichList) <- files.names

##根据GO id和KEGG id提取基因集合

##计算GO和KEGG在表达矩阵的表达活性得分

