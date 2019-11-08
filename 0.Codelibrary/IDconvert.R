#----------------------------------------
#              id转换
#----------------------------------------
#keytypes(org.Hs.eg.db)
IDconvert <- function(inputGenes,
                      from,
                      to=c("SYMBOL","ENTREZID")){
  #' @param inputGenes: 输入基因集合
  #' @param from: 该基因集合的类型
  #' @param  to: 待转换的基因类型
  library(org.Hs.eg.db)
  result <- select(org.Hs.eg.db,keys=inputGenes,keytype = from,columns = to)
  return(result)
}

GSM3674752    GSM3674752           tumoral
GSM3674753    GSM3674753             polyp
GSM3674754    GSM3674754             polyp
GSM3674755    GSM3674755           tumoral
GSM3674756    GSM3674756             polyp
GSM3674757    GSM3674757           tumoral