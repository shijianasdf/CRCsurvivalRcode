#'---------------------------------------------------------------------------------
#'                       
#'                       处理结直肠癌临床数据
#'                       
#'---------------------------------------------------------------------------------
#"D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/data_bcr_clinical_data_patient - 副本.txt" 最全的字段
#"D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/data_bcr_clinical_data_patient.txt"
colorectal.clinic.data.whole <- read.table("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/cbioportal/data_bcr_clinical_data_patient_whole.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F);
options(stringsAsFactors = F)
colorectal.clinic.data <- read.table("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/cbioportal/data_bcr_clinical_data_patient.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F);
table(colorectal.clinic.data$RESIDUAL_TUMOR,useNA = "always")
#读入firehorse数据，提取微卫星不稳定性信息并合并两个表格
{
  colorectal.clinic.data.firehorse <- read.table("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/gdac.broadinstitute.org_COADREAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0/All_CDEs.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F);
  head(colorectal.clinic.data.firehorse[1:5,1:5])
  mss_pos <- which(colorectal.clinic.data.firehorse[,1]=="mononucleotide_and_dinucleotide_marker_panel_analysis_status")
  head(t(colorectal.clinic.data.firehorse[c(1,mss_pos),]))
  colorectal.clinic.data.firehorse <- t(colorectal.clinic.data.firehorse[c(1,mss_pos),])[-1,]
  head(colorectal.clinic.data.firehorse)
  class(colorectal.clinic.data.firehorse)
  colorectal.clinic.data.firehorse <- as.data.frame(colorectal.clinic.data.firehorse)
  colorectal.clinic.data.firehorse[,1] <- toupper(colorectal.clinic.data.firehorse[,1])
  head(colorectal.clinic.data.firehorse)
  temp <- merge(colorectal.clinic.data,colorectal.clinic.data.firehorse,by.x="PATIENT_ID",by.y="1",all.x=TRUE)
  head(temp)
  colnames(temp)[26] <- "MSS_state"
  table(temp$MSS_state)
  temp$MSS_state[which(temp$MSS_state=="indeterminate")] <- NA
  temp$MSS_state <- toupper(temp$MSS_state)
  colorectal.clinic.data <- temp 
}

#规范化处理需要的临床信息
{
  table(colorectal.clinic.data[,6],useNA="always")
  colorectal.clinic.data [grep("T1", colorectal.clinic.data[,6] ),6] <- "T1";
  colorectal.clinic.data [grep("T2", colorectal.clinic.data[,6] ),6] <- "T2";
  colorectal.clinic.data [grep("T3", colorectal.clinic.data[,6] ),6] <- "T3";
  colorectal.clinic.data [grep("T4", colorectal.clinic.data[,6] ),6] <- "T4";
  colorectal.clinic.data [grep("TX", colorectal.clinic.data[,6] ),6] <- NA;
  colorectal.clinic.data [grep("Tis", colorectal.clinic.data[,6] ),6] <- NA;
  table(colorectal.clinic.data[,7],useNA="always")
  colorectal.clinic.data [grep("N0", colorectal.clinic.data[,7] ),7] <- "N0";
  colorectal.clinic.data [grep("N1", colorectal.clinic.data[,7] ),7] <- "N1";
  colorectal.clinic.data [grep("N2", colorectal.clinic.data[,7] ),7] <- "N2";
  colorectal.clinic.data [grep("N3", colorectal.clinic.data[,7] ),7] <- "N3";
  colorectal.clinic.data [grep("NX", colorectal.clinic.data[,7] ),7] <- NA;
  table(colorectal.clinic.data[,8],useNA="always")
  colorectal.clinic.data [grep("M0", colorectal.clinic.data[,8] ),8] <- "M0";
  colorectal.clinic.data [grep("M1", colorectal.clinic.data[,8] ),8] <- "M1";
  colorectal.clinic.data [grep("MX", colorectal.clinic.data[,8] ),8] <- NA;
  table(colorectal.clinic.data[,9],useNA="always")
  colorectal.clinic.data [grep("Stage IV", colorectal.clinic.data[,9] ),9] <- "Stage 4";
  colorectal.clinic.data [grep("Stage III", colorectal.clinic.data[,9] ),9] <- "Stage 3";
  colorectal.clinic.data [grep("Stage II", colorectal.clinic.data[,9] ),9] <- "Stage 2";
  colorectal.clinic.data [grep("Stage I", colorectal.clinic.data[,9] ),9] <- "Stage 1";
  colorectal.clinic.data [grep("Stage X", colorectal.clinic.data[,9] ),9] <- NA;
  colorectal.clinic.data [-c(grep("Stage", colorectal.clinic.data[,9] )),9] <- NA;
  table(colorectal.clinic.data$RESIDUAL_TUMOR,useNA="always")
  colorectal.clinic.data$RESIDUAL_TUMOR[colorectal.clinic.data$RESIDUAL_TUMOR=="RX"] <- NA
  #生存数据的清洗，早期先进性excel表格清洗
  table(colorectal.clinic.data$OS_STATUS,useNA="always")
  colorectal.clinic.data$OS_STATUS[which(colorectal.clinic.data$OS_STATUS == "LIVING")] <- 0;
  colorectal.clinic.data$OS_STATUS[which(colorectal.clinic.data$OS_STATUS == "DECEASED")] <- 1;
  table(colorectal.clinic.data$DFS_STATUS,useNA="always")
  colorectal.clinic.data$DFS_STATUS[which(colorectal.clinic.data$DFS_STATUS == "DiseaseFree")] <- 0;
  colorectal.clinic.data$DFS_STATUS[which(colorectal.clinic.data$DFS_STATUS == "Recurred/Progressed")] <- 1;
}
#整合tumour_site信息,根据tumour_site设置
{
  #读入结直肠癌肿瘤位置注释信息
  tumour_site <- read.csv("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/tumour_site.csv")
  table(colorectal.clinic.data$PRIMARY_SITE,useNA = "always")
  select <- match(colorectal.clinic.data$PRIMARY_SITE,tumour_site$anatomic_organ_subdivision)
  colorectal.clinic.data$tumour_site <- tumour_site$tumour_site[select]
  table(colorectal.clinic.data$tumour_site,useNA = "always")
}
{
  tumour_site <- read.csv("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/tumour_site_2.csv")
  table(colorectal.clinic.data$PRIMARY_SITE,useNA = "always")
  select <- match(colorectal.clinic.data$PRIMARY_SITE,tumour_site$anatomic_organ_subdivision)
  colorectal.clinic.data$tumour_site_2 <- tumour_site$tumour_site[select]
  table(colorectal.clinic.data$tumour_site_2,useNA = "always")
}
#整合亚型信息(分别来自2篇文章29622466，22810696)
{
  subtypes <- read.csv("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/coadREAD.subtype.csv")
  keep_subtype <- match(colorectal.clinic.data$PATIENT_ID,subtypes$pan.samplesID)
  colorectal.clinic.data$subtype <- subtypes$Subtype_Selected[keep_subtype]
  library(stringr)
  colorectal.clinic.data$subtype <- str_split(colorectal.clinic.data$subtype,"\\.",simplify = TRUE)[,2]
  colorectal.clinic.data$subtype[which(colorectal.clinic.data$subtype == "")] <- NA
  table(colorectal.clinic.data$subtype,useNA = "always")
  
  library(TCGAbiolinks)
  dataSubt <- TCGAquery_subtype(tumor = "COAD")
  dataSubt1 <- TCGAquery_subtype(tumor = "READ")
  subtypes <- unique(rbind.data.frame(dataSubt,dataSubt1))
  keep_subtype <- match(colorectal.clinic.data$PATIENT_ID,subtypes$patient)
  colorectal.clinic.data$subtype1 <- as.character(subtypes$methylation_subtype)[keep_subtype]
  colorectal.clinic.data$subtype1[which(colorectal.clinic.data$subtype1== "NA")] <- NA
  table(colorectal.clinic.data$subtype1,useNA = "always")
}
##提取需要的临床字段c(1,2,26,6,7,8,9,10,19,21,22,23,24,25)
{
  colnames(colorectal.clinic.data)
  head(colorectal.clinic.data,1)
  head(colorectal.clinic.data[,c(1,2,26,6,7,8,9,10,19,21,22,23,24,25,27,28,29)])
  colorectal.clinic.data <- colorectal.clinic.data[,c(1,2,26,6,7,8,9,10,19,21,22,23,24,25,27,28,29)]
}
##存储处理好的临床数据
save(colorectal.clinic.data,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData")
write.csv(colorectal.clinic.data,row.names = F, file="D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.csv")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData")
head(colorectal.clinic.data)
# PATIENT_ID    SEX                      RACE HISTORY_NEOADJUVANT_TRTYN AJCC_STAGING_EDITION
# 1 TCGA-3L-AA1B Female BLACK OR AFRICAN AMERICAN                        No                  7th
# 2 TCGA-4N-A93T   Male BLACK OR AFRICAN AMERICAN                        No                  7th
# 3 TCGA-4T-AA8H Female BLACK OR AFRICAN AMERICAN                        No                  7th
# 4 TCGA-5M-AAT4   Male BLACK OR AFRICAN AMERICAN                        No                  6th
# 5 TCGA-5M-AAT6 Female BLACK OR AFRICAN AMERICAN                        No                  6th
# 6 TCGA-5M-AATE   Male BLACK OR AFRICAN AMERICAN                        No                  7th
# AJCC_TUMOR_PATHOLOGIC_PT AJCC_NODES_PATHOLOGIC_PN AJCC_METASTASIS_PATHOLOGIC_PM
# 1                       T2                       N0                            M0
# 2                       T4                       N1                            M0
# 3                       T3                       N0                          <NA>
#   4                       T3                       N0                            M1
# 5                       T4                       N2                            M1
# 6                       T3                       N0                            M0
# AJCC_PATHOLOGIC_TUMOR_STAGE RESIDUAL_TUMOR VITAL_STATUS DAYS_TO_LAST_FOLLOWUP DAYS_TO_DEATH
# 1                     Stage 1             R0        Alive                   154            NA
# 2                     Stage 3             R0        Alive                     8            NA
# 3                     Stage 2             R0        Alive                   160            NA
# 4                     Stage 4             R0         Dead                  <NA>            49
# 5                     Stage 4             R2         Dead                  <NA>           290
# 6                     Stage 2             R0        Alive                  1200            NA
# KRAS_GENE_ANALYSIS_INDICATOR KRAS_MUTATION BRAF_GENE_ANALYSIS_INDICATOR
# 1                           NO          <NA>                           NO
# 2                           NO          <NA>                           NO
# 3                           NO          <NA>                           NO
# 4                           NO          <NA>                           NO
# 5                           NO          <NA>                         <NA>
#   6                          YES           YES                           NO
# BRAF_GENE_ANALYSIS_RESULT RADIATION_TREATMENT_ADJUVANT AGE     PRIMARY_SITE
# 1                      <NA>                         <NA>  61            Cecum
# 2                      <NA>                           NO  67  Ascending Colon
# 3                      <NA>                           NO  42 Descending Colon
# 4                      <NA>                           NO  74  Ascending Colon
# 5                      <NA>                           NO  40 Transverse Colon
# 6                      <NA>                           NO  76  Ascending Colon
# TUMOR_TISSUE_SITE OS_STATUS OS_MONTHS DFS_STATUS DFS_MONTHS MSS_state      tumour_site
# 1             Colon      <NA>        NA       <NA>         NA       MSS      right colon
# 2             Colon         0      4.80          0       4.80       MSS      right colon
# 3             Colon         0     12.65          0      12.65      <NA>       left colon
# 4             Colon         1      1.61       <NA>         NA       MSS      right colon
# 5             Colon         1      9.53          1       7.19     MSI-H Transverse Colon
# 6             Colon         0     39.42          1      26.61       MSS      right colon
# subtype
# 1     CIN
# 2     CIN
# 3     CIN
# 4     CIN
# 5     MSI
# 6     CIN
dim(colorectal.clinic.data)
#[1] 629  30