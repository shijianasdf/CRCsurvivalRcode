#'---------------------------------------------------
#'  93个驱动基因突变的MutationMap
#'---------------------------------------------------
#COAD克隆突变数据
coad.com.data <- get(load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/COAD.mut.cff.RData")); 
#READ克隆突变数据
read.com.data <- get(load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/READ.mut.cff.RData"));

#观察coad和read的样本号是不是有一样的,结论没有
intersect(coad.com.data$patient,read.com.data$patient)
#合并coad和read的克隆突变数据
colorectal.com.data <- rbind(coad.com.data,read.com.data)
dim(colorectal.com.data)
head(colorectal.com.data)
#对TCGA样本编号进行处理,选择需要的列
#patient、Variant_Classification、Hugo_Symbol、TCGA.purity、absolute.ccf、comb.timing、CI95.timing、prob.clonal.timing
colorectal.select.table <- colorectal.com.data[, c(1,2,3,4,5,6,15,18,21)];
colorectal.select.table$patient <- substring(colorectal.select.table$patient,1,12);
#删除Silent变异
colorectal.select.table <- colorectal.select.table[which(!(colorectal.select.table$Variant_Classification == "Silent")),];
#colorectal克隆突变数据
save(colorectal.select.table,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/mutationmap.data.RData");


#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/mutationmap.data.RData")
#加载结直肠癌driver gene
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData")
#加载结直肠癌临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData")

intersect.patients <- intersect(colorectal.clinic.data$PATIENT_ID,colorectal.select.table$patient) #临床数据有629个病人，536个病人
length(intersect.patients) #536个交集病人
colorectal.select.table <- colorectal.select.table[colorectal.select.table$patient %in% intersect.patients,]
colorectal.select.table <- colorectal.select.table[colorectal.select.table$Hugo_Symbol %in% ColorectalDriverGenes,]
head(colorectal.select.table)
table(colorectal.select.table$Variant_Classification)
length(unique(colorectal.select.table$mutation_id)) #4167个突变
length(unique(colorectal.select.table$patient)) #513个样本在结直肠癌driver基因上有突变
table(colorectal.select.table$CI95.timing) 
#Clonal Subclonal 
#3304       863
table(colorectal.select.table$patient)
table(colorectal.select.table$Hugo_Symbol) 
unique(colorectal.select.table$Hugo_Symbol)  #93个驱动基因在513个样本中有变异

#cbioportal mutation map需要的数据结构
# Sample_ID	Cancer_Type	Chromosome	Start_Position	End_Position	Reference_Allele	Variant_Allele
# TCGA-49-4494-01	Lung Adenocarcinoma	7	55249071	55249071	C	T
# TCGA-L9-A50W-01	Lung Adenocarcinoma	7	55249071	55249071	C	T

# Hugo_Symbol	HUGO symbol for the gene	TP53
# Protein_Change	Amino acid change	V600E
# Sample_ID	Tumor sample ID	TCGA-B5-A11E
# Mutation_Type	Translational effect of variant allele	Missense_Mutation, Nonsense_Mutation, etc.
# Chromosome	Chromosome number	X, Y, M, 1, 2, etc.
# Start_Position	Lowest numeric position of the reported variant on the genomic reference sequence	666
# End_Position	Highest numeric position of the reported variant on the genomic reference sequence	667
# Reference_Allele	The plus strand reference allele at this position	A
# Variant_Allele	Tumor sequencing (discovery) allele	C
# Validation_Status	Second pass results from orthogonal technology	Valid
# Mutation_Status	Mutation status	Somatic, Germline, etc.
# Center	Center/Institute reporting the variant	mskcc.org

head(colorectal.select.table)
# patient                   mutation_id Reference_Base Alternate_Base Hugo_Symbol Variant_Classification absolute.ccf
# 25  TCGA-3L-AA1B TCGA-3L-AA1B-01:14:65544658:G              G              A         MAX      Missense_Mutation         1.00
# 62  TCGA-3L-AA1B TCGA-3L-AA1B-01:1:115258747:C              C              T        NRAS      Missense_Mutation         1.00
# 111 TCGA-3L-AA1B TCGA-3L-AA1B-01:3:178936091:G              G              A      PIK3CA      Missense_Mutation         1.00
# 126 TCGA-3L-AA1B TCGA-3L-AA1B-01:5:112173917:C              C              T         APC      Nonsense_Mutation         1.00
# 127 TCGA-3L-AA1B TCGA-3L-AA1B-01:5:112175639:C              C              T         APC      Nonsense_Mutation         0.56
# 147 TCGA-3L-AA1B TCGA-3L-AA1B-01:7:140453146:G              G              C        BRAF      Missense_Mutation         1.00
# CI95.timing prob.clonal.timing
# 25       Clonal             Clonal
# 62       Clonal             Clonal
# 111      Clonal             Clonal
# 126      Clonal             Clonal
# 127   Subclonal          Subclonal
# 147      Clonal             Clonal
clone.gene.list <- split(colorectal.select.table,colorectal.select.table$Hugo_Symbol)
head(clone.gene.list[[1]],40)
#突变类型注释信息 https://genome-preview.ucsc.edu/cgi-bin/hgVai
table(colorectal.select.table$Variant_Classification)
mapper <- read.csv("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/mapper.csv",sep=",",header=T,fill=T,stringsAsFactors = F)
colorectal.select.table$Variant_Classification[grep("UTR",colorectal.select.table$Variant_Classification)] <- "other"
colorectal.select.table$Variant_Classification[grep("Frame_Shift",colorectal.select.table$Variant_Classification)] <- "frameshift"


tempList <- lapply(clone.gene.list,function(x){
  # Hugo_Symbol Sample_ID Mutation_Type Chromosome Start_Position End_Position Reference_Allele Variant_Allele
  # cancertype <- rep("colorectal",nrow(x))
  library(stringr)
  temp <- str_split(x$mutation_id,":",simplify = T)
  Chromosome <- temp[,2]
  Start_Position <- temp[,3]
  End_Position <- temp[,3]
  tt <- cbind.data.frame(x[,c(5,1,6)],Chromosome,Start_Position,End_Position,x[,c(3,4)])
  colnames(tt) <- c("Hugo_Symbol","Sample_ID","Mutation_Type","Chromosome","Start_Position","End_Position","Reference_Allele","Variant_Allele")
  tt
})
write.table(tempList[[1]],file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/test.txt",sep="\t",quote=F,row.names=F,col.names=T)




