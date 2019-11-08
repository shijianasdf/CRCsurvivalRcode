#下载TCGA-COAD read-count
library(TCGAbiolinks)
{
  query <- GDCquery(project = "TCGA-COAD", 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  GDCdownload(query, method = "api", files.per.chunk = 100,directory = "D:/Rsources/Project/TCGA-COAD/ReadCount")
  expdat <- GDCprepare(query = query,directory = "F:/Rsources/Project/TCGA-COAD/ReadCount")
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste("D:/Rsources/Project/TCGA-COAD/","TCGACOAD","Counts.csv",sep = "-"))
}
{
  query <- GDCquery(project = "TCGA-READ", 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  GDCdownload(query, method = "api", files.per.chunk = 100,directory = "D:/Rsources/Project/TCGA-READ/ReadCount")
  expdat <- GDCprepare(query = query,directory = "D:/Rsources/Project/TCGA-READ/ReadCount")
  count_matrix=assay(expdat)
  colData(expdat)
  rowRanges(expdat)
  write.csv(count_matrix,file = paste("D:/Rsources/Project/TCGA-READ/","TCGAREAD","Counts.csv",sep = "-"))
}