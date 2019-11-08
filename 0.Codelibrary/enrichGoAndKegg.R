#-------------------------
#     GO和KEGG功能注释
#-------------------------
enrichGoAndKEGG <- function(geneset,
                            geneAll){
  #' @param geneset: 基因集合，是gene Entrez id
  #' @param geneAll: 背景基因集合,一般是表达谱中的基因
  #' @return 上调和下调基因富集KEGG和GO的结果
  library(DOSE)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  go_enrich_results <- lapply(c('BP','MF','CC'),function(ont){
                                cat(paste('Now process',ont))
                                enrichGO(gene = geneset, #输入的基因集合必须是gene ID
                                                universe  = geneAll, #表达谱中所有的基因
                                                OrgDb = org.Hs.eg.db,
                                                ont = ont,
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                qvalueCutoff  = 0.05,
                                                readable  = TRUE)
                        })
  cat('KEGG')
  kegg_enrich_results <- enrichKEGG(gene = geneset, #输入的基因集合必须是gene ID
                                    organism = 'hsa',
                                    universe = geneAll, #表达谱中所有的基因
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff =0.05)
  enrich_results <- list(go_enrich_results=go_enrich_results,kegg_enrich_results=kegg_enrich_results)
  return(enrich_results)
}

enrichPlot <- function(cpClassList,
                       Go.filePath.csv,
                       Go.filePath.pdf,
                       KEGG.filePath.csv,
                       KEGG.filepath.pdf){
  #' @param cpClassList: enrichGoAndKEGG返回结果
  #' @param filePath: 文件和图片输出地址
  library(clusterProfiler)
  #输出KEGG富集结果文件和图片
  pdf(KEGG.filepath.pdf,height = 12,width = 20)
  print(dotplot(cpClassList[[2]],showCategory=50))
  dev.off()
  KEGG_results<- data.frame(cpClassList[[2]])
  write.csv(KEGG_results,KEGG.filePath.csv)
  
  #输出GO富集结果文件和图片
  GoList <- cpClassList[[1]]
  out.csvs <- file.path(dirname(Go.filePath.csv),paste(c('BP','MF','CC'),basename(Go.filePath.csv),sep = "_"))
  out.pdfs <- file.path(dirname(Go.filePath.pdf),paste(c('BP','MF','CC'),basename(Go.filePath.pdf),sep = "_"))
  for(i in 1:3){
    pdf(out.pdfs[i],height = 12,width = 20)
    print(dotplot(GoList[[i]],showCategory=50))
    dev.off()
    Go_results<- data.frame(GoList[[i]])
    write.csv(Go_results,out.csvs[i])
  }
}

#基因集合富集癌症hallmark函数
enrichHallMark <- function(inputgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path,pdf.path){
  #' @param  inputgenes 一个向量 输入基因集合
  #' @param  universe 一个向量 背景基因集合，表达谱有表达的基因
  #' @param  TERM2GENE 一个数据框 通路基因集合，第一列代表通路的名字，第二列表示该通路所具有的的基因名字
  library(clusterProfiler)
  enrich.result.ob <- enricher(inputgenes, pvalueCutoff = p.cut, pAdjustMethod = "BH", universe,
                               minGSSize = min, maxGSSize = max, qvalueCutoff = q.cut, TERM2GENE=TERM2GENE)
  enrich.result <- data.frame(enrich.result.ob)
  write.csv(enrich.result,file=csv.path,row.names = F)
  pdf(pdf.path)
  print(dotplot(enrich.result.ob))
  dev.off()
  return(enrich.result.ob)
}

























library(biomaRt)
listMarts()
ensembl=useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
listFilters(ensembl)[1:5,]
listAttributes(ensembl)[,1]
getBM(attributes = c('entrezgene','hgnc_symbol'), 
      filters = 'go', 
      values = 'GO:0046148', 
      mart = ensembl)

go=c("GO:0051330","GO:0000080","GO:0000114","GO:0000082")
getBM(attributes= "hgnc_symbol",
      filters="go",
      values="GO:0034368", mart=ensembl)

library(help="GO.db")
