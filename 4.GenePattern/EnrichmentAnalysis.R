#-------------------------------------------------
#对上调和下调基因分别进行GO和KEGG富集分析
#-------------------------------------------------

#加载3个基因的差异表达基因
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP7_DEG.RData")
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/SMAD2_DEG.RData")
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ANK1_DEG.RData")

#提取CAPS8上调和下调基因
GroupA <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "UP"]
GroupB <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "UP"]
GroupC <- limma_DEG$ENSG[limma_DEG$lab == "UP"]
GroupA1 <- DESeq2_DEG$ENSG[DESeq2_DEG$lab == "DOWN"]
GroupB1 <- EdgeR_DEG$ENSG[EdgeR_DEG$lab == "DOWN"]
GroupC1 <- limma_DEG$ENSG[limma_DEG$lab == "DOWN"]
Up_genes <- intersect(GroupA,GroupB)
Down_genes <- intersect(GroupA1,GroupB1)

#id转换
Up_genes_id <- IDconvert(Up_genes,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #富集基因集合
Down_genes_id <- IDconvert(Down_genes,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #富集基因集合
DESeq2_DEG_id <- IDconvert(DESeq2_DEG$ENSG,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #背景基因集合

#GO KEGG富集分析
CASP8_up_enrichGoKegg_result <- enrichGoAndKEGG(na.omit(Up_genes_id$ENTREZID),na.omit(DESeq2_DEG_id$ENTREZID))
CASP8_dowm_enrichGoKegg_result <- enrichGoAndKEGG(na.omit(Down_genes_id$ENTREZID),na.omit(DESeq2_DEG_id$ENTREZID))
#names(CASP8_up_enrichGoKegg_result) #[1] "go_enrich_results"   "kegg_enrich_results"
#存储Go KEGG富集分析结果
enrichPlot(CASP8_up_enrichGoKegg_result,"F:/CAPS8.Up.Go.csv","F:/CAPS8.Up.Go.pdf","F:/CASP8.Up.KEGG.csv","F:/CASP8.Up.KEGG.pdf")
enrichPlot(CASP8_dowm_enrichGoKegg_result,"F:/CAPS8.Down.Go.csv","F:/CAPS8.Down.Go.pdf","F:/CASP8.Down.KEGG.csv","F:/CASP8.Down.KEGG.pdf")


#提取SMAD2上调和下调基因
GroupA <- DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab == "UP"]
GroupB <- EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab == "UP"]
GroupC <- limma_DEG_SMAD2$ENSG[limma_DEG_SMAD2$lab == "UP"]
GroupA1 <- DESeq2_DEG_SMAD2$ENSG[DESeq2_DEG_SMAD2$lab == "DOWN"]
GroupB1 <- EdgeR_DEG_SMAD2$ENSG[EdgeR_DEG_SMAD2$lab == "DOWN"]
GroupC1 <- limma_DEG_SMAD2$ENSG[limma_DEG_SMAD2$lab == "DOWN"]
Up_genes <- intersect(GroupA,GroupB)
Down_genes <- intersect(GroupA1,GroupB1)

#id转换
Up_genes_id <- IDconvert(Up_genes,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #富集基因集合
Down_genes_id <- IDconvert(Down_genes,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #富集基因集合
DESeq2_DEG_SMAD2_id <- IDconvert(DESeq2_DEG_SMAD2$ENSG,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #背景基因集合

#GO KEGG富集分析
SMAD2_up_enrichGoKegg_result <- enrichGoAndKEGG(na.omit(Up_genes_id$ENTREZID),na.omit(DESeq2_DEG_SMAD2_id$ENTREZID))
SMAD2_dowm_enrichGoKegg_result <- enrichGoAndKEGG(na.omit(Down_genes_id$ENTREZID),na.omit(DESeq2_DEG_SMAD2_id$ENTREZID))
#names(CASP8_up_enrichGoKegg_result) #[1] "go_enrich_results"   "kegg_enrich_results"
#存储Go KEGG富集分析结果
enrichPlot(SMAD2_up_enrichGoKegg_result,"F:/SMAD2.Up.Go.csv","F:/SMAD2.Up.Go.pdf","F:/SMAD2.Up.KEGG.csv","F:/SMAD2.Up.KEGG.pdf")
enrichPlot(SMAD2_dowm_enrichGoKegg_result,"F:/SMAD2.Down.Go.csv","F:/SMAD2.Down.Go.pdf","F:/SMAD2.Down.KEGG.csv","F:/SMAD2.Down.KEGG.pdf")


#提取ANK1上调和下调基因
GroupA <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "UP"]
GroupB <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "UP"]
GroupC <- limma_DEG_ANK1$ENSG[limma_DEG_ANK1$lab == "UP"]
GroupA1 <- DESeq2_DEG_ANK1$ENSG[DESeq2_DEG_ANK1$lab == "DOWN"]
GroupB1 <- EdgeR_DEG_ANK1$ENSG[EdgeR_DEG_ANK1$lab == "DOWN"]
GroupC1 <- limma_DEG_ANK1$ENSG[limma_DEG_ANK1$lab == "DOWN"]
Up_genes <- intersect(GroupA,GroupB)
Down_genes <- intersect(GroupA1,GroupB1)

#id转换
Up_genes_id <- IDconvert(Up_genes,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #富集基因集合
Down_genes_id <- IDconvert(Down_genes,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #富集基因集合
DESeq2_DEG_ANK1_id <- IDconvert(DESeq2_DEG_ANK1$ENSG,from="ENSEMBL",to=c("SYMBOL","ENTREZID"))  #背景基因集合

#GO KEGG富集分析
ANK1_up_enrichGoKegg_result <- enrichGoAndKEGG(na.omit(Up_genes_id$ENTREZID),na.omit(DESeq2_DEG_ANK1_id$ENTREZID))
ANK1_dowm_enrichGoKegg_result <- enrichGoAndKEGG(na.omit(Down_genes_id$ENTREZID),na.omit(DESeq2_DEG_ANK1_id$ENTREZID))
#names(CASP8_up_enrichGoKegg_result) #[1] "go_enrich_results"   "kegg_enrich_results"
#存储Go KEGG富集分析结果
enrichPlot(ANK1_up_enrichGoKegg_result,"F:/ANK1.Up.Go.csv","F:/ANK1.Up.Go.pdf","F:/ANK1.Up.KEGG.csv","F:/ANK1.Up.KEGG.pdf")
enrichPlot(ANK1_dowm_enrichGoKegg_result,"F:/ANK1.Down.Go.csv","F:/ANK1.Down.Go.pdf","F:/ANK1.Down.KEGG.csv","F:/ANK1.Down.KEGG.pdf")


#-------------------------------------------------
#   对上调和下调基因分别进行HallMark富集分析
#-------------------------------------------------
#加载Hallmark基因集合
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/msigdb_genesets.RData")
#加载差异表达基因
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP8_DEG.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ANK1_DEG.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/SMAD2_DEG.RData")
DESeq2_DEG$lab <- ifelse(abs(DESeq2_DEG$log2FoldChange) > 1 & DESeq2_DEG$padj < 0.05,ifelse(DESeq2_DEG$log2FoldChange > 1,"UP","DOWN"),"NO")
EdgeR_DEG$lab <- ifelse(abs(EdgeR_DEG$logFC) > 1 & EdgeR_DEG$FDR < 0.05,ifelse(EdgeR_DEG$logFC > 1,"UP","DOWN"),"NO")
limma_DEG$lab <- ifelse(abs(limma_DEG$logFC) > 1 & limma_DEG$adj.P.Val < 0.05,ifelse(limma_DEG$logFC > 1,"UP","DOWN"),"NO")
DESeq2_DEG_SMAD2$lab <- ifelse(abs(DESeq2_DEG_SMAD2$log2FoldChange) > 1 & DESeq2_DEG_SMAD2$padj < 0.05,ifelse(DESeq2_DEG_SMAD2$log2FoldChange > 1,"UP","DOWN"),"NO")
EdgeR_DEG_SMAD2$lab <- ifelse(abs(EdgeR_DEG_SMAD2$logFC) > 1 & EdgeR_DEG_SMAD2$FDR < 0.05,ifelse(EdgeR_DEG_SMAD2$logFC > 1,"UP","DOWN"),"NO")
limma_DEG_SMAD2$lab <- ifelse(abs(limma_DEG_SMAD2$logFC) > 1 & limma_DEG_SMAD2$adj.P.Val < 0.05,ifelse(limma_DEG_SMAD2$logFC > 1,"UP","DOWN"),"NO")
DESeq2_DEG_ANK1$lab <- ifelse(abs(DESeq2_DEG_ANK1$log2FoldChange) > 1 & DESeq2_DEG_ANK1$padj < 0.05,ifelse(DESeq2_DEG_ANK1$log2FoldChange > 1,"UP","DOWN"),"NO")
EdgeR_DEG_ANK1$lab <- ifelse(abs(EdgeR_DEG_ANK1$logFC) > 1 & EdgeR_DEG_ANK1$FDR < 0.05,ifelse(EdgeR_DEG_ANK1$logFC > 1,"UP","DOWN"),"NO")
limma_DEG_ANK1$lab <- ifelse(abs(limma_DEG_ANK1$logFC) > 1 & limma_DEG_ANK1$adj.P.Val < 0.05,ifelse(limma_DEG_ANK1$logFC > 1,"UP","DOWN"),"NO")


temp <- msigdb_genesets$h.all
tempList <- list()
for(i in 1:length(temp)){
  tempList[[i]] <- cbind.data.frame(rep(names(temp[i]),length(temp[[i]])),temp[[i]])
}
TERM2GENE <- do.call(rbind,tempList)
colnames(TERM2GENE) <- c("term","geneID")
unique(TERM2GENE$term)
TERM2GENE$term <- as.character(TERM2GENE$term)
TERM2GENE$geneID <- as.character(TERM2GENE$geneID)
#geneID 转 geneSymbol
tt <- IDconvert(TERM2GENE$geneID,"ENTREZID","SYMBOL")
length(TERM2GENE$geneID)
dim(tt)
tt
TERM2GENE$geneID <- tt$SYMBOL
universe <- DESeq2_DEG$Gene.name   #表达谱有表达的基因 gene_symbol
#富集癌症hallmark基因集合
{
  DEseq2_CASP8upgenes <- DESeq2_DEG$Gene.name[DESeq2_DEG$lab == "UP"] #DESeq2上调基因 gene_symbol
  edgeR_CASP8upgenes <- EdgeR_DEG$Gene.name[EdgeR_DEG$lab == "UP"] #edgeR上调基因 gene_symbol
  upgenes <- intersect(DEseq2_CASP8upgenes,edgeR_CASP8upgenes)
  CASP8up_Hallmark <- enrichHallMark(upgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/CASP8up_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/CASP8up_Hallmark.pdf")
  
  DEseq2_CASP8downgenes <- DESeq2_DEG$Gene.name[DESeq2_DEG$lab == "DOWN"] #DESeq2下调基因 gene_symbol
  edgeR_CASP8downgenes <- EdgeR_DEG$Gene.name[EdgeR_DEG$lab == "DOWN"] #edgeR下调基因 gene_symbol
  downgenes <- intersect(DEseq2_CASP8downgenes,edgeR_CASP8downgenes)
  CASP8down_Hallmark <- enrichHallMark(downgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/CASP8down_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/CASP8down_Hallmark.pdf")
  
  DESeq2_SMAD2upgenes <- DESeq2_DEG_SMAD2$Gene.name[DESeq2_DEG_SMAD2$lab == "UP"] #DESeq2上调基因 gene_symbol
  edgeR_SMAD2upgenes <- EdgeR_DEG_SMAD2$Gene.name[EdgeR_DEG_SMAD2$lab == "UP"] #edgeR上调基因 gene_symbol
  SMAD2_upgenes <- intersect(DESeq2_SMAD2upgenes,edgeR_SMAD2upgenes)
  SMAD2up_Hallmark <- enrichHallMark(SMAD2_upgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/SMAD2up_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/SMAD2up_Hallmark.pdf")
  
  DESeq2_SMAD2downgenes <- DESeq2_DEG_SMAD2$Gene.name[DESeq2_DEG_SMAD2$lab == "DOWN"] #DESeq2下调基因 gene_symbol
  edgeR_SMAD2downgenes <- EdgeR_DEG_SMAD2$Gene.name[EdgeR_DEG_SMAD2$lab == "DOWN"] #edgeR下调基因 gene_symbol
  SMAD2_downgenes <- intersect(DESeq2_SMAD2downgenes,edgeR_SMAD2downgenes)
  SMAD2down_Hallmark <- enrichHallMark(SMAD2_downgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/SMAD2down_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/SMAD2down_Hallmark.pdf")
  
  DESeq2_ANK1upgenes <- DESeq2_DEG_ANK1$Gene.name[DESeq2_DEG_ANK1$lab == "UP"] #DESeq2上调基因 gene_symbol
  edgeR_ANK1upgenes <- EdgeR_DEG_ANK1$Gene.name[EdgeR_DEG_ANK1$lab == "UP"] #edgeR上调基因 gene_symbol
  ANK1_upgenes <- intersect(DESeq2_ANK1upgenes,edgeR_ANK1upgenes)
  ANK1up_Hallmark <- enrichHallMark(ANK1_upgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/ANK1up_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/ANK1up_Hallmark.pdf")
  
  DESeq2_ANK1downgenes <- DESeq2_DEG_ANK1$Gene.name[DESeq2_DEG_ANK1$lab == "DOWN"] #DESeq2上调基因 gene_symbol
  edgeR_ANK1downgenes <- EdgeR_DEG_ANK1$Gene.name[EdgeR_DEG_ANK1$lab == "DOWN"] #edgeR上调基因 gene_symbol
  ANK1_downgenes <- intersect(DESeq2_ANK1downgenes,edgeR_ANK1downgenes)
  ANK1down_Hallmark <- enrichHallMark(ANK1_downgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/ANK1down_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/ANK1down_Hallmark.pdf")
  
}

{
  DEseq2_CASP8upgenes <- DESeq2_DEG$Gene.name[DESeq2_DEG$lab == "UP"] #DESeq2上调基因 gene_symbol
  upgenes <- DEseq2_CASP8upgenes
  CASP8up_Hallmark <- enrichHallMark(upgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/CASP8up_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/CASP8up_Hallmark.pdf")
  
  DEseq2_CASP8downgenes <- DESeq2_DEG$Gene.name[DESeq2_DEG$lab == "DOWN"] #DESeq2下调基因 gene_symbol
  downgenes <- DEseq2_CASP8downgenes
  CASP8down_Hallmark <- enrichHallMark(downgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/CASP8down_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/CASP8down_Hallmark.pdf")
  
  DESeq2_SMAD2upgenes <- DESeq2_DEG_SMAD2$Gene.name[DESeq2_DEG_SMAD2$lab == "UP"] #DESeq2上调基因 gene_symbol
  SMAD2_upgenes <- DESeq2_SMAD2upgenes
  SMAD2up_Hallmark <- enrichHallMark(SMAD2_upgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/SMAD2up_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/SMAD2up_Hallmark.pdf")

  DESeq2_SMAD2downgenes <- DESeq2_DEG_SMAD2$Gene.name[DESeq2_DEG_SMAD2$lab == "DOWN"] #DESeq2下调基因 gene_symbol
  SMAD2_downgenes <- DESeq2_SMAD2downgenes
  SMAD2down_Hallmark <- enrichHallMark(SMAD2_downgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/SMAD2down_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/SMAD2down_Hallmark.pdf")

  DESeq2_ANK1upgenes <- DESeq2_DEG_ANK1$Gene.name[DESeq2_DEG_ANK1$lab == "UP"] #DESeq2上调基因 gene_symbol
  ANK1_upgenes <- DESeq2_ANK1upgenes
  ANK1up_Hallmark <- enrichHallMark(ANK1_upgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/ANK1up_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/ANK1up_Hallmark.pdf")
  
  DESeq2_ANK1downgenes <- DESeq2_DEG_ANK1$Gene.name[DESeq2_DEG_ANK1$lab == "DOWN"] #DESeq2上调基因 gene_symbol
  ANK1_downgenes <- DESeq2_ANK1downgenes
  ANK1down_Hallmark <- enrichHallMark(ANK1_downgenes,universe,TERM2GENE,p.cut=0.05,min=5,max=500,q.cut=0.2,csv.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/ANK1down_Hallmark.csv",pdf.path="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/ENSG/Hallmark/DESeq2/ANK1down_Hallmark.pdf")
}



