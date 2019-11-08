#'---------------------------------------------------
#'       93个驱动基因突变的互斥性分析
#'---------------------------------------------------
#导入OS生存数据                                                                           
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.SingleMarker.RData")
length(OS.byDriverGene.SingleMarker)
names(OS.byDriverGene.SingleMarker)
head(OS.byDriverGene.SingleMarker[[1]])
dim(OS.byDriverGene.SingleMarker[[1]])
#互斥性检验
library(dplyr)
temp <- full_join(OS.byDriverGene.SingleMarker[[1]],OS.byDriverGene.SingleMarker[[2]],
          by=c("Patient_ID","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","Residual","age","event","time","Subtype","Tumour_site"))
for(i in 3:length(OS.byDriverGene.SingleMarker)){
  #if(length(table(as.character(OS.byDriverGene.SingleMarker[[i]]$sample.label))) != 1)
    temp <- full_join(temp,OS.byDriverGene.SingleMarker[[i]],
                      by=c("Patient_ID","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","Residual","age","event","time","Subtype","Tumour_site"))
}
colnames(temp)[-(1:13)] <- names(OS.byDriverGene.SingleMarker)
temp <- temp[,!colnames(temp) %in% c("ZMYM3","DCAF12L2","BCOR","RBM10")]
temp <- temp[-(1:13)]
mint <- apply(temp,2,function(x){
  min(table(x))
})
pos <- mint >= 10
temp <- temp[,pos]
dim(temp) #536  83
head(temp)
library(forcats)
fisher.result.p <- matrix(0,nrow = 83,ncol=83,dimnames = list(colnames(temp),colnames(temp)))
fisher.result.OR <- matrix(0,nrow = 83,ncol=83,dimnames = list(colnames(temp),colnames(temp)))
for(i in 1:ncol(temp)){
  for(j in 1:ncol(temp)){
    tt <- fisher.test(table(fct_rev(temp[,i]),fct_rev(temp[,j])))
    fisher.result.p[i,j] <-  tt$p.value 
    fisher.result.OR[i,j] <-  tt$estimate
  }
}
#OR > 1 p < 0.05 - 共发生
#OR < 1 p < 0.05 + 互斥
fisher.result.q <- matrix(0,nrow = 83,ncol=83,dimnames = list(colnames(temp),colnames(temp)))
fisher.result.q <- p.adjust(as.numeric(fisher.result.p))  #默认什么方法进行多重比较
dim(fisher.result.q) <- c(83,83)
colnames(fisher.result.q) <- colnames(temp)
rownames(fisher.result.q) <- colnames(temp)
mat <- ifelse(fisher.result.q <= 0.05,ifelse(fisher.result.OR < 1,fisher.result.q,-fisher.result.q),1)
{
  library(pheatmap)
  library(RColorBrewer)
  col <- brewer.pal(11,"RdYlGn")
  #mat[upper.tri(mat)] <- NA
  pdf("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.刻画克隆性变异/mutualexclusive1.pdf")
  pheatmap(mat,cluster_rows=F,cluster_cols=F,
           border_color = "black",fontsize_row = 6,fontsize_col = 6,
           color = c("#1A9850","#006837","#D73027","#F46D43","white"),
           breaks = c(-0.05,-0.01,0,0.01,0.05,1)) #(-0.05,-0.01] (-0.01,0],(0,0.01],(0.01,0.05],(0.05,1]
  dev.off()
}


