#############################################################################
####################整合结直肠癌driver基因###################################
#############################################################################
#读取colorectal的driver genes(来自于Lawrence et al)
DriverGenes1 <- c("APC", "TP53", "FBXW7", "SMAD4", "NRAS", "SMAD2", "TCF7L2", "BRAF", "KRAS", "PIK3CA", "PCBP1", "ARID1A", "ACVR1B", "ERBB3", "CASP8", "ELF3", "TRIM23", "CDC27", "B2M", "NTN4", "AXIN2", "SIRT4", "GOT1","RBM10", "BCLAF1", "BCOR", "MAP2K4", "IDH2", "PTEN", "ERBB2", "ARID2", "CTNNB1", "TRAF3", "CNBD1", "CD70");
DriverGenes2 <- read.csv("/home/shijian2015/colorectal_driver_gene/COSMIC.CGC.20180703.tsv",fill=TRUE,stringsAsFactors = FALSE,sep="\t",header=TRUE);
extract.cancer.genes <- function(cancer.type,CGC.file)
{
  ##########提取癌症基因中somatic mutation
  CGC.cancer.genes <- CGC.file[grep(cancer.type,CGC.file[,"Tumour.Types.Somatic."]),]
  CGC.cancer.genes.symbols <- CGC.cancer.genes[,"Gene.Symbol"]
  CGC.cancer.genes.ids <- CGC.cancer.genes[,"Entrez.GeneId"]
  
  #########分致癌基因(Dom,oncogene)和抑癌基因(Rec,TSG),两者有交集
  CGC.cancer.oncogenes.symbols <- CGC.cancer.genes[grep("oncogene",CGC.cancer.genes[,"Role.in.Cancer"]),"Gene.Symbol"]
  CGC.cancer.oncogenes.ids <- CGC.cancer.genes[grep("oncogene",CGC.cancer.genes[,"Role.in.Cancer"]),"Entrez.GeneId"]
  
  CGC.cancer.TSG.symbols <- CGC.cancer.genes[grep("TSG",CGC.cancer.genes[,"Role.in.Cancer"]),"Gene.Symbol"]
  CGC.cancer.TSG.ids <- CGC.cancer.genes[grep("TSG",CGC.cancer.genes[,"Role.in.Cancer"]),"Entrez.GeneId"]
  
  ########将上述六个集合存入到RData中
  save(CGC.cancer.genes.symbols,CGC.cancer.genes.ids,CGC.cancer.oncogenes.symbols,CGC.cancer.oncogenes.ids,CGC.cancer.TSG.symbols,CGC.cancer.TSG.ids,file=paste0("D:/R/Project/预后分析/课题重复/",cancer.type,".CGC.genesymbols.geneids.RData"))
}
colon.driver.symbols <- get(load("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/colon.CGC.genesymbols.geneids.RData"));
colorectal.driver.symbols <- get(load("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/colorectal.CGC.genesymbols.geneids.RData"));
cgc.driver.symbols <- get(load("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/CRC.CGC.genesymbols.geneids.RData"));
#CGC结直肠癌相关的Driver genes
DriverGenes2 <- unique(c(colon.driver.symbols,colorectal.driver.symbols,cgc.driver.symbols));
#整合CGC和Lawrence et al的Driver Genes
ColorectalDriverGenes <- unique(union(DriverGenes1,DriverGenes2));
save(ColorectalDriverGenes,file="D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData");

#加载结直肠癌driver genes
load("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData")