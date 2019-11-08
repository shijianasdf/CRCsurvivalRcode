#-------------------------------------
#计算表达活性得分
#-------------------------------------
CaculateExprSCore <- function(X,
                              genesets,
                              method=c("gsva", "ssgsea", "zscore", "plage"),
                              kcdf=c("Gaussian", "Poisson","none")){
  #' @param X:表达矩阵
  #' @param genesets: 一个列表，储存基因集合
  #' @param method: "gsva" "ssGSEA" "Z-SCORE" "plage"
  #' @param kcdf: "Poisson" for read count "Gaussian" for microarray and log(TPM) log(cpm) log(RPKM)
  library(GSVA)
  es.max <- gsva(as.matrix(X),genesets,method=method, #表达矩阵和基因结合
                 mx.diff=FALSE, verbose=FALSE,kcdf=kcdf,
                 parallel.sz=1)
  return(es.max)
}




