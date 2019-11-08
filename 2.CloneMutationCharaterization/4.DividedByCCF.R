#-------------------------------------------------
#' @author : shijian
#' @description :  依据每个driver基因的CCF值
#'                 给病人打标签打上相应的CCF值
#-------------------------------------------------
#加载结直肠癌克隆突变数据
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData"); #colorectal.select.table
#加载结直肠癌临床数据
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData"); #colorectal.clinic.data
#加载结直肠癌driver gene
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData"); #ColorectalDriverGenes
##提取需要的临床字段c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29)
{  
  colnames(colorectal.clinic.data)
  head(colorectal.clinic.data,1)
  head(colorectal.clinic.data[,c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29,30)])
  colorectal.clinic.data <- colorectal.clinic.data[,c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29,30)]
  head(colorectal.clinic.data)
}
## 对突变数据进行清洗
colnames(colorectal.select.table)[3] <- "Gene_Symbol"
colnames(colorectal.select.table)[1] <- "Patient_ID"
## 对临床数据进行清洗
colnames(colorectal.clinic.data)[1] <- "Patient_ID"
OS.significant <- c("ANK1","ARID1A","CASP8","GRIN2A","SMAD2")
## 给每个病人基于基因CCF值打标签
grouping.byCCF.oneMarker <- function(marker.symbol, clonality.data, clinical.data){
  #STEP 1: 处理样本，保留在突变数据和临床数据中共同出现的样本
  common.sample <- intersect(unique(clonality.data$Patient_ID), clinical.data$Patient_ID);
  clinical.data <- clinical.data[clinical.data$Patient_ID %in% common.sample, ];
  clonality.data <- clonality.data[clonality.data$Patient_ID %in% common.sample, ];
  
  #STEP 2:识别marker基因是否出现在克隆性数据中
  markers <- intersect(marker.symbol,unique(clonality.data$Gene_Symbol));
  if(length(markers) == 0)
  {
    #sample.label <- rep(paste(marker.symbol,"WT"),length(common.sample));
    sample.label <- rep(0,length(common.sample));
    return(cbind.data.frame(clinical.data,sample.label));
  }  
  #STEP 3: 得到克隆、亚克隆样本(不需要)
  # clonality.cellmarker.data <- unique(clonality.data[clonality.data$Gene_Symbol %in% marker.symbol,])   
  # clone.sample <- unique(clonality.cellmarker.data[which(clonality.cellmarker.data$"CI95.timing" == "Clonal"),1])
  # subclone.sample <- unique(clonality.cellmarker.data[which(clonality.cellmarker.data$"CI95.timing" == "Subclonal"),1])
  #step3: 得到对应基因的克隆突变数据,对于每个基因，CCF取最小值
  Gene.clonality.data <- clonality.data[clonality.data$Gene_Symbol %in% marker.symbol,]
  tempList <- split(Gene.clonality.data,Gene.clonality.data[,"Patient_ID"]) 
  tempMinList <- lapply(tempList,function(x){
    pos <- which.min(x$absolute.ccf)
    x[pos,]
  })
  clonality.min.data <- dplyr::bind_rows(tempMinList) #list to data.frame
  clinical.datas <- dplyr::full_join(clinical.data,clonality.min.data,by="Patient_ID")
  clinical.data <- clinical.datas[,!colnames(clinical.datas) %in% c("mutation_id","Gene_Symbol","Variant_Classification","CI95.timing","prob.clonal.timing")]
  clinical.data$absolute.ccf[which(is.na(clinical.data$absolute.ccf))] <- 0
  return(clinical.data)
}
# temp <- grouping.byCCF.oneMarker(ColorectalDriverGenes[1],colorectal.select.table,colorectal.clinic.data)
CRC.clinical.CCF.data.list <- lapply(ColorectalDriverGenes,grouping.byCCF.oneMarker,colorectal.select.table,colorectal.clinic.data)
length(CRC.clinical.CCF.data.list)
names(CRC.clinical.CCF.data.list) <- ColorectalDriverGenes 
for(i in names(CRC.clinical.CCF.data.list)){
  colnames(CRC.clinical.CCF.data.list[[i]])[18] <- i; #改名字
}
save(CRC.clinical.CCF.data.list,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.clinical.CCF.data.list.RData")



















