#--------------------------
#   TCGA表达数据处理
#--------------------------

dataProcessing <- function(exprSet1,exprSet2,mart_export){
  #去除重复ENSG id号
  pos <- !duplicated(exprSet1[,1]) 
  exprSet1 <- exprSet1[pos,] 
  pos1 <- !duplicated(exprSet2[,1]) 
  exprSet2 <- exprSet2[pos1,]
  
  #合并表达矩阵
  library(dplyr) 
  COADREAD.ReadCount.matrix <- full_join(exprSet1,exprSet2,by="X") #overlap "X"
  rownames(COADREAD.ReadCount.matrix) <- COADREAD.ReadCount.matrix[,1]
  COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[,-1]
  
  #去除癌旁组织样本和血液样本和技术重复的样本
  group_list <- ifelse(as.numeric(substr(colnames(COADREAD.ReadCount.matrix),14,15)) < 10,'tumor','normal')
  pos <- group_list %in% "tumor"
  COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[,pos] 
  COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[,!duplicated(substr(colnames(COADREAD.ReadCount.matrix),1,12))] 
  colnames(COADREAD.ReadCount.matrix) <- gsub("\\.","-",colnames(COADREAD.ReadCount.matrix))
  
  #取临床样本和克隆突变样本交集的病人
  intersect.patients <- intersect(colorectal.clinic.data$PATIENT_ID,colorectal.select.table$patient)
  pos <- substr(colnames(COADREAD.ReadCount.matrix),1,12) %in% intersect.patients
  COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[,pos]
  
  #ENSG注释信息，筛选蛋白质编码基因
  mart_export <- mart_export[mart_export$Gene.type=="protein_coding",]
  pos <- rownames(COADREAD.ReadCount.matrix) %in% mart_export$Gene.stable.ID 
  COADREAD.ReadCount.matrix <- COADREAD.ReadCount.matrix[pos,]
  return(COADREAD.ReadCount.matrix)
}

IDmap <- function(DEG,mart_export){
  #' @description 差异表达基因id转换
  library(dplyr)
  DEG$ENSG <- rownames(DEG)
  temp <- inner_join(DEG,mart_export,by=c("ENSG" = "Gene.stable.ID"))
  return(temp)
}








