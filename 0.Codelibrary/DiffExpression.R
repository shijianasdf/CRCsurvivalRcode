#-----------------------------
#差异分析modual
#-----------------------------

makeDiffData <- function(exprSet,
                         clinical.data){
  #'@description:对表达数据基因和样本进行过滤,筛选至少在50%样本中有表达的样本,
  #'             提取有临床数据的表达样本，同时对处理好的表达数据样本打标签，用于后续差异分析    
  #'@param exprSet: 表达谱数据
  #'@param clinical.data: 临床数据
  temp.clinical <- clinical.data[,c("Patient_ID","sample.label")]
  #取表达矩阵和临床数据的病人交集
  intersect.patients <- intersect(substring(colnames(exprSet),1,12),temp.clinical$Patient_ID)
  temp.exprSet <- exprSet[,substring(colnames(exprSet),1,12) %in% intersect.patients]
  temp.clinical <- temp.clinical[temp.clinical$Patient_ID %in% intersect.patients,]
  pos <- match(substring(colnames(temp.exprSet),1,12),temp.clinical$Patient_ID)
  #生成和表达矩阵一致的样本标签
  temp.clinical <- temp.clinical[pos,]
  #给样本打标签,用于后续差异分析
  group_list <- ifelse(grepl("WT",temp.clinical$sample.label),"WT","MUT")
  #过滤表达矩阵的基因，每个基因至少在50%的样本有表达才留下
  pos <- rowSums(temp.exprSet > 0) >= ncol(temp.exprSet)/2 
  temp.exprSet <- temp.exprSet[pos,]
  return(list(exprSet=temp.exprSet,clinical.data=temp.clinical,group_list=group_list))
}


DESeq2.function <- function(ExpreSsionMatrix,
                            group_list,
                            case,
                            control){
  #' @param ExpreSsionMatrix: 表达矩阵
  #' @param group_list: 分组标签
  #' @param case： 分组标签中的处理组
  #' @param control： 分组标签中的对照组
  suppressMessages(library(DESeq2))
  (colData <- data.frame(row.names=colnames(ExpreSsionMatrix), group_list=group_list))
  #产生DESeq对象
  dds <- DESeqDataSetFromMatrix(countData = ExpreSsionMatrix,
                                colData = colData,
                                design = ~ group_list)
  #计算差异
  dds <- DESeq(dds)
  #替换异常值
  #dds <- replaceOutliersWithTrimmedMean(dds)
  #提取结果case相对于control的差异表达
  res <- results(dds, contrast=c("group_list",case,control))  
  resOrdered <- res[order(res$padj),]
  DEG <- as.data.frame(resOrdered) #转换成data.frame
  DESeq2_DEG <- na.omit(DEG)
  return(DESeq2_DEG)
}


EdgeR.function <- function(exprSet,
                           group_list){
  #' @param ExpreSsionMatrix: 表达矩阵
  #' @param group_list: 分组标签
  suppressMessages(library(edgeR))
  dge <- DGEList(counts=exprSet,group=factor(group_list))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  design
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit,  contrast=c(1,-1))
  nrDEG <- topTags(lrt, n=nrow(exprSet))
  nrDEG <- as.data.frame(nrDEG)
  edgeR_DEG <- nrDEG 
  return(edgeR_DEG)
}


limma.voom <- function(exprSet,
                       group_list,
                       case.control){
  #' @param ExpreSsionMatrix: 表达矩阵
  #' @param group_list: 分组标签,与表达矩阵的列名对应上
  #' @param case.control： "MUT-WT",MUT相对于WT算差异
  suppressMessages(library(limma))
  suppressMessages(library(edgeR))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  v <- voom(dge,design,normalize="quantile")
  fit <- lmFit(v, design)
  cont.matrix <- makeContrasts(contrasts=c(case.control),levels = design) #'trt-untrt'
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  tempOutput <- topTable(fit2, coef=case.control, n=Inf) #'trt-untrt'
  DEG_limma_voom <- na.omit(tempOutput)
  return(DEG_limma_voom)
}


  






























