load()

# Chromosome Start_Position End_Position Reference_Allele Variant_Allele
library(stringr)
by(colorectal.select.table,colorectal.select.table$Hugo_Symbol,function(x){
  str_split(x$mutation_id,":")    
})
x <- colorectal.select.table[which(colorectal.select.table$Hugo_Symbol == "APC"),]
str_split(x$mutation_id,":") 

#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
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

#对变异进行归类并剔除一些变异
table(colorectal.select.table$Variant_Classification)

# Missense Indel Nonsense Splicing
colorectal.select.table$Variant_Classification[grep("In_Frame|Frame_Shift",colorectal.select.table$Variant_Classification)] <- "Indel"
colorectal.select.table$Variant_Classification[grep("Splice",colorectal.select.table$Variant_Classification)] <- "Splicing"
colorectal.select.table$Variant_Classification[grep("Missense_Mutation",colorectal.select.table$Variant_Classification)] <- "Missense"
colorectal.select.table$Variant_Classification[grep("Nonsense_Mutation",colorectal.select.table$Variant_Classification)] <- "Nonsense"
# 其他变异类型变成others展示
colorectal.select.table$Variant_Classification[grep("UTR",colorectal.select.table$Variant_Classification)] <- "others"
colorectal.select.table$Variant_Classification[grep("Flank",colorectal.select.table$Variant_Classification)] <- "others"
colorectal.select.table$Variant_Classification[grep("Intron",colorectal.select.table$Variant_Classification)] <- "others"
colorectal.select.table$Variant_Classification[grep("Nonstop_Mutation",colorectal.select.table$Variant_Classification)] <- "others"
colorectal.select.table$Variant_Classification[grep("Translation_Start_Site",colorectal.select.table$Variant_Classification)] <- "others"
table(colorectal.select.table$Variant_Classification)

# cbioportal MutationMaper
head()




