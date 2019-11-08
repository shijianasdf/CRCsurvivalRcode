#'--------------------------------------------------------------------------------------------------
#'比较每个驱动基因克隆突变和亚克隆突变的变异类型差异，检验基因的克隆性和变异类型的关系
#'---------------------------------------------------------------------------------------------------
#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
#加载结直肠癌driver gene
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData")
#加载结直肠癌临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData")
##筛选具有临床数据的结直肠癌克隆突变数据，对结直肠癌克隆突变数据依据driver gene进行过滤
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
# 3'Flank                  3'UTR                5'Flank                  5'UTR 
# 2                    119                      1                     39 
# Frame_Shift_Del        Frame_Shift_Ins           In_Frame_Del           In_Frame_Ins 
# 441                    164                     31                      5 
# Intron      Missense_Mutation      Nonsense_Mutation       Nonstop_Mutation 
# 141                   2460                    670                      1 
# Splice_Site Translation_Start_Site 
# 88                      5
#Missense Indel Nonsense Splicing
colorectal.select.table$Variant_Classification[grep("In_Frame|Frame_Shift",colorectal.select.table$Variant_Classification)] <- "Indel"
colorectal.select.table$Variant_Classification[grep("Splice",colorectal.select.table$Variant_Classification)] <- "Splicing"
colorectal.select.table$Variant_Classification[grep("Missense_Mutation",colorectal.select.table$Variant_Classification)] <- "Missense"
colorectal.select.table$Variant_Classification[grep("Nonsense_Mutation",colorectal.select.table$Variant_Classification)] <- "Nonsense"
#过滤掉其余变异3'Flank 3'UTR  5'Flank  5'UTR Intron Nonstop_Mutation Translation_Start_Site
colorectal.select.table <- colorectal.select.table[-grep("UTR",colorectal.select.table$Variant_Classification),]
colorectal.select.table <- colorectal.select.table[-grep("Flank",colorectal.select.table$Variant_Classification),]
colorectal.select.table <- colorectal.select.table[-grep("Intron",colorectal.select.table$Variant_Classification),]
colorectal.select.table <- colorectal.select.table[-grep("Nonstop_Mutation",colorectal.select.table$Variant_Classification),]
colorectal.select.table <- colorectal.select.table[-grep("Translation_Start_Site",colorectal.select.table$Variant_Classification),]
unique(colorectal.select.table$Variant_Classification)

#将其余变异转变为others进行分析

write.table(colorectal.select.table,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.table.txt",sep="\t",quote=F,row.names=F,col.names=F)

##########绘制四种突变类型的克隆和亚克隆突变频率柱图#########
#根据基因分组
head(colorectal.select.table)
# patient                   mutation_id Hugo_Symbol Variant_Classification absolute.ccf CI95.timing prob.clonal.timing
# 25  TCGA-3L-AA1B TCGA-3L-AA1B-01:14:65544658:G         MAX               Missense         1.00      Clonal             Clonal
# 62  TCGA-3L-AA1B TCGA-3L-AA1B-01:1:115258747:C        NRAS               Missense         1.00      Clonal             Clonal
# 111 TCGA-3L-AA1B TCGA-3L-AA1B-01:3:178936091:G      PIK3CA               Missense         1.00      Clonal             Clonal
# 126 TCGA-3L-AA1B TCGA-3L-AA1B-01:5:112173917:C         APC               Nonsense         1.00      Clonal             Clonal
# 127 TCGA-3L-AA1B TCGA-3L-AA1B-01:5:112175639:C         APC               Nonsense         0.56   Subclonal          Subclonal
# 147 TCGA-3L-AA1B TCGA-3L-AA1B-01:7:140453146:G        BRAF               Missense         1.00      Clonal             Clonal
clone.gene.list <- split(colorectal.select.table,colorectal.select.table$Hugo_Symbol)
#第一种柱图，感觉画错了
{
  for(i in 1:length(clone.gene.list)){
    title <- names(clone.gene.list[i])
    tt <- as.data.frame(table(factor(clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),factor(clone.gene.list[[i]]$Variant_Classification,levels=c("Missense","Nonsense","Indel","Splicing"))))
    temp <- table(clone.gene.list[[i]]$CI95.timing)
    tt[tt$Var1 == "Clonal",3] <- tt[tt$Var1 == "Clonal",3]/temp["Clonal"]
    tt[tt$Var1 == "Subclonal",3] <- tt[tt$Var1 == "Subclonal",3]/temp["Subclonal"]
    colnames(tt) <- c("clonality","muttype","Freq")
    library(RColorBrewer)
    col <- brewer.pal(9,"Blues")[c(2,4,6,8)] #表示使用Blues的2,4,6,8颜色,这样可以区分开
    library(ggplot2)
    pdf(paste0("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/",title,".clonalAndVariable.pdf"),width = 9,height=9) #width height default 7
    group.barplot <- ggplot(data=tt,mapping=aes(x=muttype,y=Freq,fill=factor(clonality,levels=c("Subclonal","Clonal"))))+ #levels=c("Subclonal","Clonal","Clonal-Subclonal" 调整颜色顺序
      geom_bar(stat="identity",width=0.5,position=position_dodge())+ #default width 0.9 position="stack" position=position_dodge(0.8)
      scale_fill_manual(values = col[2:3])+ #设置颜色fill=factor(group,levels=c("Subclonal","Clonal","Clonal-Subclonal")) Subclonal对应col[1] Clonal对应col[2] Clonal-Subclonal对应col[3]
      theme_bw()+ 
      labs(title = title,x="Mutation Type",y="Frequency %",fill="clone status")+ #更改y轴题目和图例题目
      #x坐标文本字体角度angle，vjust垂直移动多少，hjust水平移动多少，size设置字体大小正常的0.4倍
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size = rel(1)))+
      theme(legend.position = "right",legend.key.size = unit(0.5, "cm"),legend.title=element_text(size=9),legend.text=element_text(size=9))  #调整图例位置,大小,调整图例标题大小,调整图例内部字体大小
      print(group.barplot) 
    dev.off()
  }
}
#克隆状态和突变类型的关系的整体费舍尔精确检验
{
  fisher.result <- list()
  fisher.result.p <- list()
  for(i in 1:length(clone.gene.list)){
    fisher.result[[i]] <- fisher.test(table(factor(clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),factor(clone.gene.list[[i]]$Variant_Classification,levels=c("Missense","Nonsense","Indel","Splicing"))))
    fisher.result.p[[i]] <- fisher.result[[i]]$p.value
  }
  names(fisher.result) <- names(clone.gene.list)
  names(fisher.result.p) <- names(clone.gene.list)
  unlist(fisher.result.p)[unlist(fisher.result.p) <= 0.05]
  save(fisher.result,fisher.result.p,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/fisher.result.RData")
  load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/fisher.result.RData")
}
#第二种柱图
{ ##   supp dose  len
  ## 1   VC D0.5  6.8
  ## 2   VC   D1 15.0
  ## 3   VC   D2 33.0
  ## 4   OJ D0.5  4.2
  ## 5   OJ   D1 10.0
  ## 6   OJ   D2 29.5
  head(clone.gene.list[[1]])
  table(factor(clone.gene.list[[1]]$CI95.timing,levels=c("Clonal","Subclonal")),factor(clone.gene.list[[1]]$Variant_Classification,levels=c("Missense","Nonsense","Indel","Splicing")))
  tt <- as.data.frame(table(factor(clone.gene.list[[1]]$CI95.timing,levels=c("Clonal","Subclonal")),factor(clone.gene.list[[1]]$Variant_Classification,levels=c("Missense","Nonsense","Indel","Splicing"))))
  tt[tt$Var2 == "Missense",3] <- tt[tt$Var2 == "Missense",3]/sum(tt[tt$Var2 == "Missense",3])
  tt[tt$Var2 == "Nonsense",3] <- tt[tt$Var2 == "Nonsense",3]/sum(tt[tt$Var2 == "Nonsense",3])
  tt[tt$Var2 == "Indel",3] <- tt[tt$Var2 == "Indel",3]/sum(tt[tt$Var2 == "Indel",3])
  tt[tt$Var2 == "Splicing",3] <- tt[tt$Var2 == "Splicing",3]/sum(tt[tt$Var2 == "Splicing",3])
  tt$Freq[is.nan(tt$Freq)] <- 0
  colnames(tt) <- c("clonality","muttype","Freq")
  for(i in 1:length(clone.gene.list)){
    title <- names(clone.gene.list[i])
    tt <- as.data.frame(table(factor(clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),factor(clone.gene.list[[i]]$Variant_Classification,levels=c("Missense","Nonsense","Indel","Splicing"))))
    #temp <- table(clone.gene.list[[i]]$CI95.timing)
    tt[tt$Var2 == "Missense",3] <- tt[tt$Var2 == "Missense",3]/sum(tt[tt$Var2 == "Missense",3])
    tt[tt$Var2 == "Nonsense",3] <- tt[tt$Var2 == "Nonsense",3]/sum(tt[tt$Var2 == "Nonsense",3])
    tt[tt$Var2 == "Indel",3] <- tt[tt$Var2 == "Indel",3]/sum(tt[tt$Var2 == "Indel",3])
    tt[tt$Var2 == "Splicing",3] <- tt[tt$Var2 == "Splicing",3]/sum(tt[tt$Var2 == "Splicing",3])
    tt$Freq[is.nan(tt$Freq)] <- 0
    colnames(tt) <- c("clonality","muttype","Freq")
    library(RColorBrewer)
    col <- brewer.pal(9,"Blues")[c(2,4,6,8)] #表示使用Blues的2,4,6,8颜色,这样可以区分开
    library(ggplot2)
    pdf(paste0("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/第二种柱图/",title,".clonalAndVariable.pdf"),width = 9,height=9) #width height default 7
    group.barplot <- ggplot(data=tt,mapping=aes(x=muttype,y=Freq,fill=factor(clonality,levels=c("Subclonal","Clonal"))))+ #levels=c("Subclonal","Clonal","Clonal-Subclonal" 调整颜色顺序
      geom_bar(stat="identity")+ #default width 0.9 position="stack" position=position_dodge(0.8)
      scale_fill_manual(values = col[2:3])+ #设置颜色fill=factor(group,levels=c("Subclonal","Clonal","Clonal-Subclonal")) Subclonal对应col[1] Clonal对应col[2] Clonal-Subclonal对应col[3]
      theme_bw()+ 
      labs(title = title,x="Mutation Type",y="Frequency %",fill="clone status")+ #更改y轴题目和图例题目
      #x坐标文本字体角度angle，vjust垂直移动多少，hjust水平移动多少，size设置字体大小正常的0.4倍
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size = rel(1)))+
      theme(legend.position = "right",legend.key.size = unit(0.5, "cm"),legend.title=element_text(size=9),legend.text=element_text(size=9))  #调整图例位置,大小,调整图例标题大小,调整图例内部字体大小
    print(group.barplot) 
    dev.off()
  }
}
#Missense vs NonMissense
#将Indel Nonsense Splicing变为others
Missense.colorectal.select.table <- colorectal.select.table
Missense.colorectal.select.table$Variant_Classification[grep("Indel",Missense.colorectal.select.table$Variant_Classification)] <- "others"
Missense.colorectal.select.table$Variant_Classification[grep("Nonsense",Missense.colorectal.select.table$Variant_Classification)] <- "others"
Missense.colorectal.select.table$Variant_Classification[grep("Splicing",Missense.colorectal.select.table$Variant_Classification)] <- "others"

Missense.clone.gene.list <- split(Missense.colorectal.select.table,Missense.colorectal.select.table$Hugo_Symbol)
sink(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/MissensesVsOthers.fisher.result.txt")
MissensesVsOthers.fisher.result <- list()
MissensesVsOthers.fisher.result.p <- list()
for(i in 1:length(Missense.clone.gene.list)){
  print(names(Missense.clone.gene.list[i]))
  print(table(factor(Missense.clone.gene.list[[i]]$Variant_Classification,levels=c("Missense","others")),factor(Missense.clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),dnn=c("muttype","clonality")))
  MissensesVsOthers.fisher.result[[i]] <- fisher.test( table(factor(Missense.clone.gene.list[[i]]$Variant_Classification,levels=c("Missense","others")),factor(Missense.clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),dnn=c("muttype","clonality")) )
  MissensesVsOthers.fisher.result.p[[i]] <- MissensesVsOthers.fisher.result[[i]]$p.value
  print(MissensesVsOthers.fisher.result[[i]])
  print(MissensesVsOthers.fisher.result.p[[i]])
}
names(MissensesVsOthers.fisher.result) <- names(Missense.clone.gene.list)
names(MissensesVsOthers.fisher.result.p) <- names(Missense.clone.gene.list)
sink()
save(MissensesVsOthers.fisher.result,MissensesVsOthers.fisher.result.p,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/MissensesVsOthers.fisher.result.RData")
unlist(MissensesVsOthers.fisher.result.p)[unlist(MissensesVsOthers.fisher.result.p) <= 0.05]

#Nonsense vs NonMissense
#将Indel Missense Splicing变为others   Indel Missense Nonsense Splicing
Nonsense.colorectal.select.table <- colorectal.select.table
Nonsense.colorectal.select.table$Variant_Classification[grep("Indel",Nonsense.colorectal.select.table$Variant_Classification)] <- "others"
Nonsense.colorectal.select.table$Variant_Classification[grep("Missense",Nonsense.colorectal.select.table$Variant_Classification)] <- "others"
Nonsense.colorectal.select.table$Variant_Classification[grep("Splicing",Nonsense.colorectal.select.table$Variant_Classification)] <- "others"

Nonsense.clone.gene.list <- split(Nonsense.colorectal.select.table,Nonsense.colorectal.select.table$Hugo_Symbol)
sink("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/NonsenseVsOthers.fisher.result.txt")
NonsenseVsOthers.fisher.result <- list()
NonsenseVsOthers.fisher.result.p <- list()
for(i in 1:length(Nonsense.clone.gene.list)){
  print(names(Nonsense.clone.gene.list[i])) 
  print(table(factor(Nonsense.clone.gene.list[[i]]$Variant_Classification,levels=c("Nonsense","others")),factor(Nonsense.clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),dnn=c("muttype","clonality"))) 
  NonsenseVsOthers.fisher.result[[i]] <- fisher.test( table(factor(Nonsense.clone.gene.list[[i]]$Variant_Classification,levels=c("Nonsense","others")),factor(Nonsense.clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),dnn=c("muttype","clonality")) )
  NonsenseVsOthers.fisher.result.p[[i]] <- NonsenseVsOthers.fisher.result[[i]]$p.value
  print(NonsenseVsOthers.fisher.result[[i]])
  print(NonsenseVsOthers.fisher.result.p[[i]])
}
names(NonsenseVsOthers.fisher.result) <- names(Nonsense.clone.gene.list)
names(NonsenseVsOthers.fisher.result.p) <- names(Nonsense.clone.gene.list)
sink()
save(NonsenseVsOthers.fisher.result,NonsenseVsOthers.fisher.result.p,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/NonsenseVsOthers.fisher.result.RData")
unlist(NonsenseVsOthers.fisher.result.p)[unlist(NonsenseVsOthers.fisher.result.p) <= 0.05]


#Indel vs NonIndel
#将Nonsense Missense Splicing变为others   Indel Missense Nonsense Splicing
Indel.colorectal.select.table <- colorectal.select.table
Indel.colorectal.select.table$Variant_Classification[grep("Nonsense",Indel.colorectal.select.table$Variant_Classification)] <- "others"
Indel.colorectal.select.table$Variant_Classification[grep("Missense",Indel.colorectal.select.table$Variant_Classification)] <- "others"
Indel.colorectal.select.table$Variant_Classification[grep("Splicing",Indel.colorectal.select.table$Variant_Classification)] <- "others"
table(Indel.colorectal.select.table$Variant_Classification)

Indel.clone.gene.list <- split(Indel.colorectal.select.table,Indel.colorectal.select.table$Hugo_Symbol)
sink("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/IndelVsOthers.fisher.result.txt")
IndelVsOthers.fisher.result <- list()
IndelVsOthers.fisher.result.p <- list()
for(i in 1:length(Indel.clone.gene.list)){
  print(names(Indel.clone.gene.list[i]))
  print(table(factor(Indel.clone.gene.list[[i]]$Variant_Classification,levels=c("Indel","others")),factor(Indel.clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),dnn=c("muttype","clonality")))
  IndelVsOthers.fisher.result[[i]] <- fisher.test( table(factor(Indel.clone.gene.list[[i]]$Variant_Classification,levels=c("Indel","others")),factor(Indel.clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),dnn=c("muttype","clonality")) )
  IndelVsOthers.fisher.result.p[[i]] <- IndelVsOthers.fisher.result[[i]]$p.value
  print(IndelVsOthers.fisher.result[[i]])
  print(IndelVsOthers.fisher.result.p[[i]])
}
names(IndelVsOthers.fisher.result) <- names(Indel.clone.gene.list)
names(IndelVsOthers.fisher.result.p) <- names(Indel.clone.gene.list)
sink()
save(IndelVsOthers.fisher.result,IndelVsOthers.fisher.result.p,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/IndelVsOthers.fisher.result.RData")
unlist(IndelVsOthers.fisher.result.p)[unlist(IndelVsOthers.fisher.result.p) <= 0.05]


#Splicing vs NonSplicing
#将Nonsense Missense Splicing变为others   Indel Missense Nonsense Splicing
Splicing.colorectal.select.table <- colorectal.select.table
Splicing.colorectal.select.table$Variant_Classification[grep("Nonsense",Splicing.colorectal.select.table$Variant_Classification)] <- "others"
Splicing.colorectal.select.table$Variant_Classification[grep("Missense",Splicing.colorectal.select.table$Variant_Classification)] <- "others"
Splicing.colorectal.select.table$Variant_Classification[grep("Indel",Splicing.colorectal.select.table$Variant_Classification)] <- "others"
table(Splicing.colorectal.select.table$Variant_Classification)

Splicing.clone.gene.list <- split(Splicing.colorectal.select.table,Splicing.colorectal.select.table$Hugo_Symbol)
sink("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/SplicingVsOthers.fisher.result.txt")
SplicingVsOthers.fisher.result <- list()
SplicingVsOthers.fisher.result.p <- list()
for(i in 1:length(Splicing.clone.gene.list)){
  print(names(Splicing.clone.gene.list[i]))
  print(table(factor(Splicing.clone.gene.list[[i]]$Variant_Classification,levels=c("Splicing","others")),factor(Splicing.clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),dnn=c("muttype","clonality")))
  SplicingVsOthers.fisher.result[[i]] <- fisher.test( table(factor(Splicing.clone.gene.list[[i]]$Variant_Classification,levels=c("Splicing","others")),factor(Splicing.clone.gene.list[[i]]$CI95.timing,levels=c("Clonal","Subclonal")),dnn=c("muttype","clonality")) )
  SplicingVsOthers.fisher.result.p[[i]] <- SplicingVsOthers.fisher.result[[i]]$p.value
  print(SplicingVsOthers.fisher.result[[i]])
  print(SplicingVsOthers.fisher.result.p[[i]])
}
names(SplicingVsOthers.fisher.result) <- names(Splicing.clone.gene.list)
names(SplicingVsOthers.fisher.result.p) <- names(Splicing.clone.gene.list)
sink()
save(SplicingVsOthers.fisher.result,SplicingVsOthers.fisher.result.p,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/SplicingVsOthers.fisher.result.RData")
unlist(SplicingVsOthers.fisher.result.p)[unlist(SplicingVsOthers.fisher.result.p) <= 0.05]

load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/fisher.result.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/SplicingVsOthers.fisher.result.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/MissensesVsOthers.fisher.result.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/NonsenseVsOthers.fisher.result.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/cloneAndMutationType/显著的/IndelVsOthers.fisher.result.RData")
unlist(fisher.result.p)[unlist(fisher.result.p) <= 0.05]
unlist(MissensesVsOthers.fisher.result.p)[unlist(MissensesVsOthers.fisher.result.p) <= 0.05]
unlist(SplicingVsOthers.fisher.result.p)[unlist(SplicingVsOthers.fisher.result.p) <= 0.05]
unlist(IndelVsOthers.fisher.result.p)[unlist(IndelVsOthers.fisher.result.p) <= 0.05]
unlist(NonsenseVsOthers.fisher.result.p)[unlist(NonsenseVsOthers.fisher.result.p) <= 0.05]


fisher.test(matrix(c(11,15,0,9),nrow=2))
