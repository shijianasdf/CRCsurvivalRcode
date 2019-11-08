##############探索克隆突变数据#####################
#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
#加载结直肠癌driver gene
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData")
#加载结直肠癌临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData");

###############对结直肠癌克隆突变数据，所有的突变进行探索和统计
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
setdiff(ColorectalDriverGenes,unique(colorectal.select.table$Hugo_Symbol)) #"RBM10","BCOR","DCAF12L2","ZMYM3"四个基因在样本中无变异
##统计具有不同突变状态的病人百分比
State <- by(colorectal.select.table,colorectal.select.table$patient,function(x){ifelse(length(unique(x$CI95.timing))==2,"clonal-subclonal",ifelse(unique(x$CI95.timing)=="Clonal","Clonal","Subclonal"))})
State <- as.character(State)
table(State)
library(ggplot2)
library(scales)#引入percent函数
library(RColorBrewer)
col <- brewer.pal(9,"Blues")[c(2,4,6,8)] #表示使用Blues的2,4,6,8颜色,这样可以区分开
pie.df <- data.frame(group = c("unmutated", "clonal-subclonal", "clonal","subclonal"),
                     value = c(536-513, 283, 224,6))
pie<- ggplot(pie.df, aes(x="", y=value, fill=factor(group,levels = c("unmutated", "subclonal","clonal","clonal-subclonal"))))+ #fill设置因子调颜色顺序
            geom_bar(stat = "identity",colour="black",size=1)+ #position=position_dodge(0.8),
            coord_polar(theta ="y", start=0)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )      
pie.out <- pie + scale_fill_manual(values = col) + blank_theme +
                 theme(axis.text.x=element_blank())+
                 labs(title = "% cases based on mutations:",x = NULL, y = NULL, fill = NULL)+
                 geom_text(aes(y = value/4 + c(0, cumsum(value)[-length(value)]), 
                 label = percent(value/sum(value))),size=5) #调文本位置,暂时没办法，用AI调吧
pie.out
ggsave(pie.out,filename = "D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.DriverGeneClone/pie.pdf")

## 统计93个driver基因的突变类型(克隆突变和亚克隆突以及)所在的样本数,因为4个基因没有突变
# Create some data
clone.gene.list <- split(colorectal.select.table,colorectal.select.table$Hugo_Symbol)
#temp.frame <- data.frame(gene=character(0),group=character(0),cases=charater(0))
temp.list <- lapply(clone.gene.list,function(x){
  temp.frame1 <- data.frame(patient=character(0),label=character(0))
  for(i in 1:length(unique(x$patient))){
    pos <- which(x$patient == unique(x$patient)[i])
    temp <- x[pos,c("patient","CI95.timing")]
    label <- ifelse(length(unique(temp$"CI95.timing"))==2,"Clonal-Subclonal",ifelse(unique(temp$"CI95.timing")=="Clonal","Clonal","Subclonal"))
    #num <- length(pos)
    temp1 <- cbind.data.frame(unique(x$patient)[i],label)
    temp.frame1 <- rbind.data.frame(temp.frame1,temp1)
  }
  tt <- table(factor(temp.frame1$label,levels=c("Clonal","Subclonal","Clonal-Subclonal")))
  temp.frame <- cbind.data.frame(gene=rep(unique(x$Hugo_Symbol),3),cases=tt)
  colnames(temp.frame)[2:3] <- c("group","cases")
  temp.frame
})
#待可视化的数据
df <- do.call(rbind,temp.list)
save(df,file = "D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.DriverGeneClone/")
# Create the group barplot
#按照克隆突变数排序
head(df)
level <- as.character(df[which(df$group == "Clonal"),]$gene[order(df[which(df$group == "Clonal"),]$cases,decreasing = T)])
df$gene <- factor(df$gene,levels=level)  #基因字段设置成因子类型
#按照亚克隆突变数排序
level_subclonal <- as.character(df[which(df$group == "Subclonal"),]$gene[order(df[which(df$group == "Subclonal"),]$cases,decreasing = T)])
df$gene <- factor(df$gene,levels=level)  #基因字段设置成因子类型
#按照克隆亚克隆图变数排序
level_subclonal_clonal <- as.character(df[which(df$group == "Clonal-Subclonal"),]$gene[order(df[which(df$group == "Clonal-Subclonal"),]$cases,decreasing = T)])
df$gene <- factor(df$gene,levels=level)  #基因字段设置成因子类型
library(ggplot2)
pdf("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/group.barplot1.pdf",width = 9,height=9) #width height default 7
group.barplot <- ggplot(data=df, mapping=aes(x=gene, y=cases, fill=factor(group,levels=c("Subclonal","Clonal","Clonal-Subclonal")))) + #levels=c("Subclonal","Clonal","Clonal-Subclonal" 调整颜色顺序
                    geom_bar(stat="identity",width=0.9,position=position_dodge(0.9)) + #default width 0.9 position="stack" position=position_dodge(0.8)
                    scale_fill_manual(values = col[2:4]) + #设置颜色fill=factor(group,levels=c("Subclonal","Clonal","Clonal-Subclonal")) Subclonal对应col[1] Clonal对应col[2] Clonal-Subclonal对应col[3]
                    theme_bw() + labs(x="gene",y="cases",fill="clone status") + #更改y轴题目和图例题目
                    #x坐标文本字体角度angle，vjust垂直移动多少，hjust水平移动多少，size设置字体大小正常的0.4倍
                    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size = rel(0.6))) +
                    theme(legend.position = "top",legend.key.size = unit(0.5, "cm")) +  #调整图例位置，大小
                    theme(legend.title=element_text(size=9)) +  #调整图例标题大小
                    theme(legend.text=element_text(size=9)) + #调整图例内部字体大小
                    theme(aspect.ratio = 0.3) #面板的长宽比
                    #coord_flip()
print(group.barplot) 
dev.off()
#ggsave(group.barplot,filename = "D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/group.barplot.pdf")

## 统计93个driver基因的突变类型(克隆突变和亚克隆突以及)所在的样本比例,因为4个基因没有突变
# Create some data
clone.gene.list <- split(colorectal.select.table,colorectal.select.table$Hugo_Symbol)
#统计每个基因clonal突变频率所占所有突变的比例
clonal.frequency <- lapply(clone.gene.list,function(x){
                              mutation_num <- table(x$CI95.timing)
                              pp <- mutation_num["Clonal"]/sum(mutation_num)
                              return(pp)
                    })
fivenum(as.numeric(clonal.frequency))

#temp.frame <- data.frame(gene=character(0),group=character(0),cases=charater(0))
temp.list <- lapply(clone.gene.list,function(x){
                      temp.frame1 <- data.frame(patient=character(0),label=character(0))
                      for(i in 1:length(unique(x$patient))){
                         pos <- which(x$patient == unique(x$patient)[i])
                         temp <- x[pos,c("patient","CI95.timing")]
                         label <- ifelse(length(unique(temp$"CI95.timing"))==2,"Clonal-Subclonal",ifelse(unique(temp$"CI95.timing")=="Clonal","Clonal","Subclonal"))
                         #num <- length(pos)
                         temp1 <- cbind.data.frame(unique(x$patient)[i],label)
                         temp.frame1 <- rbind.data.frame(temp.frame1,temp1)
                      }
                      tt <- table(factor(temp.frame1$label,levels=c("Clonal","Subclonal","Clonal-Subclonal")))
                      temp.frame <- cbind.data.frame(gene=rep(unique(x$Hugo_Symbol),3),
                                                     cases=tt/sum(tt)*100)
                      colnames(temp.frame)[2:3] <- c("group","cases")
                      temp.frame
             })
#统计clonal突变所占的突变样本比例
fivenum(as.numeric(lapply(temp.list,function(x){sum(x$cases[grep("Clonal",as.character(x$group))])})))


#待可视化的数据
df <- do.call(rbind,temp.list)
head(df)
save(df,file = "D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.DriverGeneClone/")
# Create the stack barplot
colnames(df)
pdf("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/stack.barplot1.pdf",width = 9,height=9) #width height default 7
stack.barplot <- ggplot(data=df, mapping=aes(x=gene, y=cases, fill=factor(group,levels=c("Subclonal","Clonal","Clonal-Subclonal")))) + #levels=c("Subclonal","Clonal","Clonal-Subclonal" 调整颜色顺序
                      geom_bar(stat="identity",width=0.7,position="stack") + #default width 0.9 position="stack" position=position_dodge(0.8)
                      scale_fill_manual(values = col[2:4]) +
                      theme_bw() + labs(fill="clone status") +
                      #x坐标文本字体角度angle，vjust垂直移动多少，hjust水平移动多少，size设置字体大小正常的0.6倍
                      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size = rel(0.6)),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
                      theme(legend.position = "top",legend.key.size = unit(0.5, "cm")) +  #调整图例位置，大小
                      theme(legend.title=element_text(size=9)) +  #调整图例标题大小
                      theme(legend.text=element_text(size=9)) + #调整图例内部字体大小
                      labs(y="% cases") + #y轴题目
                      theme(aspect.ratio = 0.3)
print(stack.barplot)
dev.off()

##统计97个driver基因每个变异Distribution of the CCF,4个基因没有突变
##creat dotplot data
head(colorectal.select.table)
col
pdf("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/dotplot.pdf",width = 9,height=9) #width height default 7
dot.plot <- ggplot(data=colorectal.select.table, mapping=aes(x=Hugo_Symbol, y=absolute.ccf,fill=factor(CI95.timing,levels = c("Subclonal","Clonal")))) + 
                geom_jitter(shape= 21,alpha = 1,size=1,width=0.2,color="black",stroke=0.4)+ #You'll have to use shapes from 21 to 25. These are the ones that have colour and fill properties
                labs(x="gene",y="CCF",fill = "clone status") + #设置legend题目为"Legend Title"
                theme_bw() + scale_fill_manual(values=col[c(1,4)],labels = c("Subclonal","Clonal")) +  #要先设置默认的主题风格
                theme(aspect.ratio = 0.3,legend.position = "top",legend.key.size = unit(0.5, "cm"),legend.title = element_text(size = 9),legend.text=element_text(size=9),axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size = rel(0.6)),panel.grid.major = element_line(size = 0.1),panel.grid.minor = element_line(size = 0.1))
print(dot.plot)
dev.off()
#ggsave(dot.plot,width=unit(15,"in"),height=unit(4,"in"),filename = "D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/dotplot.pdf")

##########探索突变互斥,只是一个可视化############
#Each column represents an colorectal case carrying at least 1 mutation in any of the studied genes.
{
  #所有变异数量不为0的样本
  head(colorectal.select.table)
  #要生成一个矩阵，行是基因(93个)，列是样本(513个)，不同颜色对应不同的变异类型，来进行可视化
  #1:WT 2:subclonal 3:clonal 4:subclonal-clonal
  heatmap.matrix <- matrix(0,nrow=93,ncol=513)
  rownames(heatmap.matrix) <- unique(colorectal.select.table$Hugo_Symbol)
  colnames(heatmap.matrix) <- unique(colorectal.select.table$patient)
  #step1 循环基因，判断该基因在每个样本的克隆状态，从而确定对应数值
  Hugo.list <- split(colorectal.select.table,colorectal.select.table$Hugo_Symbol)
  lapply(Hugo.list,function(x){
    
  })
  #step2 循环基因，判断该基因在每个样本的突变状态，从而确定对应数值
  #1:WT 2:MUT
  Hugo.list <- split(colorectal.select.table,colorectal.select.table$Hugo_Symbol)
  Hugo.patients.list <- lapply(Hugo.list,function(x){unique(x$patient)})
  Heatmap.row.list <- lapply(Hugo.patients.list,function(x){sign(colnames(heatmap.matrix) %in% x)+1})  #判断该基因在哪个病人中发生突变
  rownames(heatmap.matrix) <- names(Heatmap.row.list)
  for(i in 1:dim(heatmap.matrix)[1]){
    heatmap.matrix[i,] <- Heatmap.row.list[[i]]
  }
  #统计每个基因的突变频率
  Gene_frep <- apply(heatmap.matrix,1,function(x){length(which(x==2))/length(x)})
  cut_off <- 10/513
  Gene_frep <- Gene_frep[Gene_frep >= cut_off] #筛选至少在10个样本中发生突变的基因
  length(Gene_frep) #得到83个基因
  heatmap.matrix <- heatmap.matrix[rownames(heatmap.matrix) %in% names(Gene_frep), ]
  
  table(heatmap.matrix[1,],heatmap.matrix[2,])
  chisq.test(table(heatmap.matrix[1,],heatmap.matrix[2,]))
  fisher.test(table(heatmap.matrix[1,],heatmap.matrix[2,]))
}
# 用complexheatmap绘制突变热图矩阵
# patient                   mutation_id Hugo_Symbol Variant_Classification absolute.ccf CI95.timing
# 25  TCGA-3L-AA1B TCGA-3L-AA1B-01:14:65544658:G         MAX      Missense_Mutation         1.00      Clonal
# 62  TCGA-3L-AA1B TCGA-3L-AA1B-01:1:115258747:C        NRAS      Missense_Mutation         1.00      Clonal
# 111 TCGA-3L-AA1B TCGA-3L-AA1B-01:3:178936091:G      PIK3CA      Missense_Mutation         1.00      Clonal
# 126 TCGA-3L-AA1B TCGA-3L-AA1B-01:5:112173917:C         APC      Nonsense_Mutation         1.00      Clonal
# 127 TCGA-3L-AA1B TCGA-3L-AA1B-01:5:112175639:C         APC      Nonsense_Mutation         0.56   Subclonal
# 147 TCGA-3L-AA1B TCGA-3L-AA1B-01:7:140453146:G        BRAF      Missense_Mutation         1.00      Clonal
# prob.clonal.timing
# 25              Clonal
# 62              Clonal
# 111             Clonal
# 126             Clonal
# 127          Subclonal
# 147             Clonal
heatmap.matrix <- matrix(0,nrow=93,ncol=513)
rownames(heatmap.matrix) <- unique(colorectal.select.table$Hugo_Symbol)
colnames(heatmap.matrix) <- unique(colorectal.select.table$patient)
G_list <- split(colorectal.select.table,colorectal.select.table$Hugo_Symbol)
G_p_list <- lapply(G_list,function(x){
            temp <- c()
            for(i in unique(x$patient)){
               a <- ifelse(length(unique(x$CI95.timing[which(x$patient==i)]))==2,"ClonalSubclonal;",ifelse(unique(x$CI95.timing[which(x$patient==i)])=="Clonal","Clonal;","Subclonal;"))
               temp[i] <- a
            }
            temp
})
for(i in 1:length(G_p_list)){
  pos <- na.omit(match(names(G_p_list[[i]]),colnames(heatmap.matrix)))
  heatmap.matrix[names(G_p_list[i]),pos] <- unname(G_p_list[[i]])
  heatmap.matrix[names(G_p_list[i]),-pos] <- ""
}

# TCGA-05-4384-01 TCGA-05-4390-01 TCGA-05-4425-01
# KRAS "  "            "Clonal;"          "  "           
# HRAS "Subclonal;"            "  "            "  "           
# BRAF "  "            "  "            "Clonal-Subclonal;" 

# 突变可视化矩阵
mat <- heatmap.matrix
# 过滤表达频率在5%以上的基因
class(mat)
dim(mat)
rt <- apply(mat,1,function(x){
          length(x[x != ""])/513
      })
class(rt)
pos <- which(rt >= 0.05)
mat <- mat[pos,]
{
  library(ComplexHeatmap)
  library(RColorBrewer)
  cols <- brewer.pal(9,"Blues")[c(2,4,6,8)][2:4] #表示使用Blues的2,4,6,8颜色,这样可以区分开
  # "#9ECAE1" "#4292C6" "#08519C"
  #col = c("Clonal" = "red", "Subclonal" = "#008000","ClonalSubclonal"="blue")
  col = c("Clonal" = "#4292C6", "Subclonal" = "#9ECAE1", "ClonalSubclonal" = "#08519C")
  alter_fun = list(
    background = function(x, y, w, h){
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    Subclonal = function(x, y, w, h){
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#9ECAE1", col = NA))
    },
    Clonal = function(x, y, w, h){
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#4292C6", col = NA))
    },
    ClonalSubclonal = function(x, y, w, h){
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#08519C", col = NA))
    }
  )
  #In this case, we need to define a function to extract different alteration types and pass the function to get_type argument. 
  #The function should return a vector of alteration types.
  pdf("D:/Rsources/oncoprint.pdf",width = 12,height=8)
  oncoPrint( mat, get_type = function(x) strsplit(x, ";")[[1]], remove_empty_columns = TRUE,row_barplot_width = unit(10, "cm"),
             alter_fun = alter_fun, col = col, row_names_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 3),
             column_title = "OncoPrint for colorectal cancer clonality mutations",
             heatmap_legend_param = list(title = "Alternations", at = c("Clonal", "Subclonal","ClonalSubclonal"), 
                                         labels = c("Clonal", "Subclonal","ClonalSubclonal")) )
  dev.off()
}










