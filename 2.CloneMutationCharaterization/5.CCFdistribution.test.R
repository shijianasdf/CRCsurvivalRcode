#'-------------------------------------------------
#' 
#'-------------------------------------------------
#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
#加载结直肠癌driver gene
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData")
#加载结直肠癌临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData");

intersect.patients <- intersect(colorectal.clinic.data$PATIENT_ID,colorectal.select.table$patient)
colorectal.select.table <- colorectal.select.table[colorectal.select.table$patient %in% intersect.patients,]
colorectal.select.table <- colorectal.select.table[colorectal.select.table$Hugo_Symbol %in% ColorectalDriverGenes,]
head(colorectal.select.table)
dim(colorectal.select.table)
ColorectalDriverGenes[!ColorectalDriverGenes %in% unique(colorectal.select.table$Hugo_Symbol)]

KS.tests <- function(driver.gene){
  colorectal.select.table <- colorectal.select.table[colorectal.select.table$Hugo_Symbol==driver.gene,]
  if(nrow(colorectal.select.table) != 0){
    p.sig <- ks.test(colorectal.select.table$absolute.ccf,"punif")$p.value #KS检验每个基因的CCF值分布是否符合均匀分布
  }else{
    p.sig <- NULL
  }
  return(p.sig)
}
KS.tests.result <- sapply(ColorectalDriverGenes,function(x){KS.tests(x)})
KS.tests.result <- unlist(KS.tests.result)
q.value <- p.adjust(KS.tests.result,method = "BH")
q.fdr.value <- p.adjust(KS.tests.result,method = "fdr")
length(q.value)
which(q.value > 0.05) #没有一个基因突变的CCF值符合均匀分布
which(q.fdr.value > 0.05) #没有一个基因突变的CCF值符合均匀分布
which(q.value > 0.1) #没有一个基因突变的CCF值符合均匀分布
which(q.fdr.value > 0.1) #没有一个基因突变的CCF值符合均匀分布
