#############################################################################
####################整合结肠癌和结肠癌克隆突变数据###########################
#############################################################################
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
colorectal.select.table <- colorectal.com.data[, c(1,2,5,6,15,18,21)]
colorectal.select.table$patient <- substring(colorectal.select.table$patient,1,12)
#删除Silent变异(其实还应该清楚intron，5'UTR等变异)
colorectal.select.table <- colorectal.select.table[which(!(colorectal.select.table$Variant_Classification == "Silent")),]
#colorectal克隆突变数据
save(colorectal.select.table,file="D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")

#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
head(colorectal.select.table)
#        patient                    mutation_id Hugo_Symbol Variant_Classification absolute.ccf
# 2 TCGA-3L-AA1B TCGA-3L-AA1B-01:10:128193330:C    C10orf90      Missense_Mutation         1.00
# 3 TCGA-3L-AA1B  TCGA-3L-AA1B-01:10:49997997:C       WDFY4      Missense_Mutation         0.99
# 4 TCGA-3L-AA1B   TCGA-3L-AA1B-01:10:6533674:C       PRKCQ      Missense_Mutation         1.00
# 5 TCGA-3L-AA1B TCGA-3L-AA1B-01:11:103153768:C     DYNC2H1      Missense_Mutation         1.00
# 6 TCGA-3L-AA1B TCGA-3L-AA1B-01:11:132289968:C       OPCML                  3'UTR         0.59
# 7 TCGA-3L-AA1B  TCGA-3L-AA1B-01:11:13984681:G       SPON1                    RNA         0.60
#   CI95.timing prob.clonal.timing
# 2      Clonal             Clonal
# 3      Clonal          Subclonal
# 4      Clonal             Clonal
# 5      Clonal             Clonal
# 6      Clonal          Subclonal
# 7   Subclonal          Subclonal

#absolute.ccf.0.95用于判断clone状态和亚克隆状态
#absolute.ccf

