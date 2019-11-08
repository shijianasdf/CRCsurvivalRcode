#----------------------------------------------------------
#permutation test检测driver基因富集到克隆突变还是亚克隆突变
#----------------------------------------------------------
## 导入结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
## 导入driver基因信息
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData") #ColorectalDriverGenes
## 加载结直肠癌临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData");

permutation.test <- function(driver.gene,colorectal.select.table,times){
  intersect.patients <- intersect(colorectal.clinic.data$PATIENT_ID,colorectal.select.table$patient)
  colorectal.select.table <- colorectal.select.table[colorectal.select.table$patient %in% intersect.patients,]
  colorectal.select.table <- colorectal.select.table[colorectal.select.table$Hugo_Symbol %in% ColorectalDriverGenes,]
  observe.table <- colorectal.select.table[colorectal.select.table$Hugo_Symbol == driver.gene,] #基因参数
  mutation.type <- factor(observe.table$CI95.timing,levels = c("Clonal","Subclonal"))
  temp <- table(mutation.type)
  observation.clonal.ratio <- temp[1]/temp[2]
  observation.subclonal.ratio <- temp[2]/temp[1]
  num <- length(observe.table$mutation_id) #基因真实突变数
  pat <- unique(observe.table$patient) #该基因突变所在的病人集合
  background.clone.table <- colorectal.select.table[colorectal.select.table$patient %in% pat,]
  set.seed(123) #为了保证你随机生成的随机数前后一致
  random.clonal.ratios <- replicate(times,{
    pos <- sample(dim(background.clone.table)[1],num,replace = F) #随机抽一次
    random.mutation.type <- table(factor(background.clone.table[pos,"CI95.timing"],levels = c("Clonal","Subclonal")))
    random.clonal.ratio <- random.mutation.type[1]/random.mutation.type[2]
    #random.subclonal.ratio <- random.mutation.type[2]/random.mutation.type[1]
    return(random.clonal.ratio)
  })
  set.seed(123) #为了保证你随机生成的随机数前后一致
  random.subclonal.ratios <- replicate(times,{
    pos <- sample(dim(background.clone.table)[1],num,replace = F) #随机抽一次
    random.mutation.type <- table(factor(background.clone.table[pos,"CI95.timing"],levels = c("Clonal","Subclonal")))
    #random.clonal.ratio <- random.mutation.type[1]/random.mutation.type[2]
    random.subclonal.ratio <- random.mutation.type[2]/random.mutation.type[1]
    return(random.subclonal.ratio)
  })
  p.clonal.sig <- sum(random.clonal.ratios >= observation.clonal.ratio)/times
  p.subclonal.sig <- sum(random.subclonal.ratios >= observation.subclonal.ratio)/times
  result <- data.frame(obv=c(observation.clonal.ratio,observation.subclonal.ratio),p=c(p.clonal.sig,p.subclonal.sig))
  rownames(result) <- c("Clonal","Subclonal")
  return(result)
}
permutation.test <- function(driver.gene,colorectal.select.table,times){
  #' @param driver.gene: 驱动基因
  #' @param colorectal.select.table: 克隆突变数据
  #' @param times: 随机次数
  intersect.patients <- intersect(colorectal.clinic.data$PATIENT_ID,colorectal.select.table$patient)
  colorectal.select.table <- colorectal.select.table[colorectal.select.table$patient %in% intersect.patients,]
  colorectal.select.table <- colorectal.select.table[colorectal.select.table$Hugo_Symbol %in% ColorectalDriverGenes,]
  observe.table <- colorectal.select.table[colorectal.select.table$Hugo_Symbol == driver.gene,] #基因参数
  mutation.type <- factor(observe.table$CI95.timing,levels = c("Clonal","Subclonal"))
  temp <- table(mutation.type)
  observation.clonal.ratio <- temp[1]/temp[2]
  observation.subclonal.ratio <- temp[2]/temp[1]
  num <- length(observe.table$mutation_id) #基因真实突变数
  set.seed(123)   #为了保证你随机生成的随机数前后一致
  random.clonal.ratios <- replicate(times,{
    pos <- sample(dim(colorectal.select.table)[1],num,replace = F) #随机抽一次
    random.mutation.type <- table(factor(colorectal.select.table[pos,"CI95.timing"],levels = c("Clonal","Subclonal")))
    random.clonal.ratio <- random.mutation.type[1]/random.mutation.type[2]
    return(random.clonal.ratio)
  })
  set.seed(123)  #为了保证你随机生成的随机数前后一致
  random.subclonal.ratios <- replicate(times,{
    pos <- sample(dim(colorectal.select.table)[1],num,replace = F) #随机抽一次
    random.mutation.type <- table(factor(colorectal.select.table[pos,"CI95.timing"],levels = c("Clonal","Subclonal")))
    random.subclonal.ratio <- random.mutation.type[2]/random.mutation.type[1]
    return(random.subclonal.ratio)
  })
  p.clonal.sig <- sum(random.clonal.ratios >= observation.clonal.ratio)/times
  p.subclonal.sig <- sum(random.subclonal.ratios >= observation.subclonal.ratio)/times
  result <- data.frame(obv=c(observation.clonal.ratio,observation.subclonal.ratio),p=c(p.clonal.sig,p.subclonal.sig))
  rownames(result) <- c("Clonal","Subclonal")
  return(result)
}
## 检测每个driver基因突变倾向于克隆还是亚克隆
permutation.result <- lapply(ColorectalDriverGenes,function(x){
  permutation.test(x,colorectal.select.table,1000)
})
names(permutation.result) <- ColorectalDriverGenes
library(stringr)
temp <- cbind(rep(names(permutation.result),each=2),str_split(rownames(do.call(rbind.data.frame,permutation.result)),"\\.",simplify = T)[,2],do.call(rbind.data.frame,permutation.result))
colnames(temp)[1:2] <- c("gene","clonality")
clonal.pos <- seq(1,nrow(temp),by=2)
subclonal.pos <- seq(2,nrow(temp),by=2)
fdr <- vector(mode="numeric",length=nrow(temp))
clonal.fdr <- p.adjust(temp[clonal.pos,"p"],method = 'fdr') #使用fdr参数
subclonal.fdr <- p.adjust(temp[subclonal.pos,"p"],method = 'fdr') #使用fdr参数
fdr[clonal.pos] <- clonal.fdr
fdr[subclonal.pos] <- subclonal.fdr
temp <- cbind.data.frame(temp,fdr)
head(temp)
temp[temp$fdr<0.25,]
save(permutation.result,temp,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/permutation.test.result.RData")
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/permutation.test.result.RData")
temp <- na.omit(temp)
temp <- temp[!(temp$gene %in% as.character(temp$gene[is.infinite(temp$obv)])),]
library(ggplot2)
pdf("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/permutation.test.result5.pdf",width=11,height=7)
enrichment.plot <- ggplot(data=temp, aes(x=gene, y=obv, fill=factor(clonality,levels=c("Subclonal","Clonal")))) +
                        ylim(0,15)+ #y轴坐标0-15
                        geom_bar(stat="identity", position=position_dodge(),width = 0.6) +
                        geom_text(aes(label=ifelse(fdr <= 0.01,"***",ifelse(fdr <= 0.15,"**",ifelse(fdr <= 0.25,"*","")))), vjust=0.5, hjust=0.5,color="red", #fdr < 0.05 或者 < 0.1 添加文本
                                  position = position_dodge(1), size=4)+
                        scale_fill_brewer(palette="Paired")+
                        labs(fill="clonality status",y="clonality ratio",x="gene")+
                        theme_bw() + 
                        theme( aspect.ratio = 0.3, legend.position = "top", 
                               legend.key.size = unit(0.5, "cm"), 
                               legend.title = element_text(size = 9), 
                               legend.text=element_text(size=9), 
                               axis.text.x = element_text(angle = 90,size = rel(0.6)), 
                               panel.grid.major = element_line(size = 0.1), 
                               panel.grid.minor = element_line(size = 0.1) )
print(enrichment.plot)
dev.off()
#ggsave(enrichment.plot,filename ="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/permutation.test.result.pdf")
temp[which(temp$fdr <= 0.1),]
temp[which(temp$fdr <= 0.2),]
temp[which(temp$fdr <= 0.25),]


library(ggplot2)
df <- data.frame(random.ratios=random.ratios)
density.plot <- ggplot(df, aes(x=random.ratios)) + 
                          geom_histogram(aes(y=..density..), colour="black", fill="white") + #画频率分布直方图
                          geom_density(alpha=.2, fill="#FF6666") +    #画密度曲线图
                          geom_vline(xintercept = observation.ratio,  #加垂直线
                                       color = "red", size=1) +
                          theme_bw()
ggsave(density.plot,filename ="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/permutation.test.result.pdf")



