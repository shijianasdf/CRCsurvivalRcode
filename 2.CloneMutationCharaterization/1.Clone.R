###################################################################################################
#生成clonal VS subclonal VS WT的临床生存数据，并且统计每个driver基因的clone突变数和subclone突变数
###################################################################################################
#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData"); #colorectal.select.table
#加载结直肠癌临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData"); #colorectal.clinic.data
#加载结直肠癌driver gene
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData"); #ColorectalDriverGenes

##提取需要的临床字段c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29)
{
  colnames(colorectal.clinic.data)
  head(colorectal.clinic.data,1)
  head(colorectal.clinic.data[,c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29,30)])
  colorectal.clinic.data <- colorectal.clinic.data[,c(1,2,26,6,7,8,9,10,19,22,23,24,25,27,28,29,30)]
  head(colorectal.clinic.data)
}
#对突变数据进行清洗
colnames(colorectal.select.table)[3] <- "Gene_Symbol";
colnames(colorectal.select.table)[1] <- "Patient_ID";
#对临床数据进行清洗
colnames(colorectal.clinic.data)[1] <- "Patient_ID";
#基于驱动基因克隆状态对样本分成3类(Clone,Subclone,WT)
grouping.byclonality.oneMarker <- function(marker.symbol, clonality.data, clinical.data)
{    
    #STEP 1: 处理样本，保留在突变数据和临床数据中共同出现的样本
    common.sample <- intersect(unique(clonality.data$Patient_ID), clinical.data$Patient_ID);
    clinical.data <- clinical.data[clinical.data$Patient_ID %in% common.sample, ];
    clonality.data <- clonality.data[clonality.data$Patient_ID %in% common.sample, ];
    
	  #STEP 2:识别marker基因是否出现在克隆性数据中
    markers <- intersect(marker.symbol,unique(clonality.data$Gene_Symbol));
    if(length(markers) == 0)
    {
  		sample.label <- rep(paste(marker.symbol,"WT"),length(common.sample));
  		return(cbind.data.frame(clinical.data,sample.label));
    }  
	
    #STEP 3: 得到克隆、亚克隆样本
    clonality.cellmarker.data <- unique(clonality.data[clonality.data$Gene_Symbol %in% marker.symbol,])   
    clone.sample <- unique(clonality.cellmarker.data[which(clonality.cellmarker.data$"CI95.timing" == "Clonal"),1])
    subclone.sample <- unique(clonality.cellmarker.data[which(clonality.cellmarker.data$"CI95.timing" == "Subclonal"),1])
    #####找到对于某个基因即是克隆又是亚克隆的样本,在我们分析中将这样的克隆亚克隆样本去掉,并且给样本打上标签
    temp.sample <- intersect(clone.sample,subclone.sample)
    if(length(temp.sample)==0){
  	 label1 <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "));
  	 label2 <- c(paste(markers, "WT",sep = " "),paste(markers, "Subclonal",sep = " "));
  	 level <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "),paste(markers, "Subclonal",sep = " ")); #设置因子顺序，对应画图颜色
  	 sample.label1 <- label1[as.numeric(clinical.data$Patient_ID %in% clone.sample)+1];
  	 sample.label2 <- label2[as.numeric(clinical.data$Patient_ID %in% subclone.sample)+1];
  	 pos2 <- grep("Subclonal",sample.label2);
  	 sample.label1[pos2] <- paste(markers,"Subclonal",sep=" ");
  	 clinical.data$'sample.label' <- factor(sample.label1,levels=level); #数据框新建一列	 
    }else{
     label1 <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "));
	   label2 <- c(paste(markers, "WT",sep = " "),paste(markers, "Subclonal",sep = " "));
     level <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "),paste(markers, "Subclonal",sep = " ")); #设置因子顺序，对应画图颜色	 
     sample.label1 <- label1[as.numeric(clinical.data$Patient_ID %in% clone.sample)+1];
	   sample.label2 <- label2[as.numeric(clinical.data$Patient_ID %in% subclone.sample)+1];
	   pos1 <- grep("Clonal",sample.label1);
	   pos2 <- grep("Subclonal",sample.label2);
	   Clonal.Subclonal.pos <- intersect(pos1,pos2);
	   sample.label1[pos2] <- paste(markers,"Subclonal",sep=" ");
	   sample.label1 <- sample.label1[-Clonal.Subclonal.pos] #删除对于该基因即有克隆突变又有亚克隆突变的标签
	   clinical.data <- clinical.data[!(clinical.data$Patient_ID %in% temp.sample),]; #删除对于该基因即有克隆突变又有亚克隆突变的样本
	   #sample.label1[Clonal.Subclonal.pos] <- paste(markers,"WT",sep = " "); 
	   clinical.data$'sample.label' <- factor(sample.label1,levels=level); #数据框新建一列	 
	}              
   return(clinical.data)
}
CRC.survival.data.list <- lapply(ColorectalDriverGenes,function(x){
	grouping.byclonality.oneMarker(x,colorectal.select.table,colorectal.clinic.data)
});
names(CRC.survival.data.list) <- ColorectalDriverGenes
#保存最初依据克隆状态分类的临床数据
save(CRC.survival.data.list,file="D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.survival.data.list.RData");

######################画出97个driver基因的clone突变样本数以及亚克隆突变样本数######################
CRC.clonal.subclonal.data.list <- list();
for(i in names(CRC.survival.data.list)){
  CRC.clonal.subclonal.data.list[[i]] <- CRC.survival.data.list[[i]]; #对list新加一个元素
  colnames(CRC.clonal.subclonal.data.list[[i]])[18] <- i; #改名字
}
names(CRC.clonal.subclonal.data.list) <- ColorectalDriverGenes;
#保存每个基因克隆状态对样本的分类数据，用于后续克隆性和临床变量研究的关系
save(CRC.clonal.subclonal.data.list,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.clonal.subclonal.data.list.RData")

#首先得到每个基因在536个样本中的克隆突变数和亚克隆突变数
CRC.clonal.subclonal.three.plot <- list(); #三种状态都有的基因
CRC.clonal.subclonal.two.plot <- list();  #只有一种变异的基因 AXIN1 POLE SMAD3 NTHL1都是clonal突变
CRC.clonal.subclonal.one.plot <- list(); #没有克隆变异的基因 RBM10 BCOR DCAF12L2 ZMYM3
for(i in names(CRC.clonal.subclonal.data.list)){
   #temp <- table(CRC.clonal.subclonal.data.list[[i]][,26]);
   temp <- table(as.character(CRC.clonal.subclonal.data.list[[i]][,18]));
   #print(temp)
   if(length(temp) == 3){
	  CRC.clonal.subclonal.three.plot[[i]] <- temp;
   }
   if(length(temp) == 2){
      #print(i);
	  CRC.clonal.subclonal.two.plot[[i]] <- temp;
   }
   if(length(temp) == 1){
      #print(i);
	  CRC.clonal.subclonal.one.plot[[i]] <- temp;
   }
}
#################################生成3种状态的基因可视化数据#########################################
picture.three.data <-  do.call(rbind,CRC.clonal.subclonal.three.plot);
picture.three.data <- as.data.frame(picture.three.data);
picture.three.data$genes <- rownames(picture.three.data);
picture.three.data <- picture.three.data[,c(4,1,2,3)];
colnames(picture.three.data) <- c("genes","Clonal","Subclonal","WT");
picture.three.data.plot.list <- list();
for(i in 1:nrow(picture.three.data)){
  temp <- data.frame(genes=rep(picture.three.data$genes[i],2),number=as.numeric(picture.three.data[i,])[2:3],status=as.character(colnames(picture.three.data[i,]))[2:3]);
	picture.three.data.plot.list[[picture.three.data$genes[i]]] <- temp;
}
picture.three.data <- do.call(rbind,picture.three.data.plot.list);
#######依据每个基因的clonal个数大小排列
level <- as.character(picture.three.data[which(picture.three.data$status == "Clonal"),]$genes[order(picture.three.data[which(picture.three.data$status == "Clonal"),]$number,decreasing = T)]);
picture.three.data$genes <- factor(picture.three.data$genes,levels=level);  #基因按照因子顺序排列
library(ggplot2);
picture <- ggplot(data=picture.three.data,mapping=aes(genes,number,fill=status))+ #数据框 
        		   #xlim(0, 80)+
            	 geom_bar(stat="identity", width = 0.5, position = position_dodge(0.9)) + #柱的宽度
        		   scale_fill_manual(values = c("RoyalBlue", "red"), labels=c("Clonal", "Subclonal"))+
        		   geom_text(mapping = aes(label= number))+ #加文本
        		   theme_set(theme_bw())+ #背景白色
        		   theme(panel.grid.major=element_line(colour=NA))+ #画板panel的格线去掉
        		   theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90))+
        		   coord_flip() #坐标翻转 
print(picture)
################################################################################################
##########################################生成2种状态的基因可视化数据######################################
picture.two.data <-  do.call(rbind,CRC.clonal.subclonal.two.plot);
picture.two.data <- as.data.frame(picture.two.data);
picture.two.data$genes <- rownames(picture.two.data);
picture.two.data <- picture.two.data[,c(3,1,2)];
colnames(picture.two.data) <- c("genes","Clonal","WT");
picture.two.data$'Subclonal' <- rep(0,nrow(picture.two.data));
picture.two.data <- picture.two.data[,c(1,2,4,3)];
picture.two.data.plot.list <- list();
for(i in 1:nrow(picture.two.data)){
  temp <- data.frame(genes=rep(picture.two.data$genes[i],2),number=as.numeric(picture.two.data[i,])[2:3],status=as.character(colnames(picture.two.data[i,]))[2:3]);
	picture.two.data.plot.list[[picture.two.data$genes[i]]] <- temp;
}
picture.two.data <- do.call(rbind,picture.two.data.plot.list);
##############################可视化两种状态的基因##################################
library(ggplot2);
#图层执行有顺序
picture <- ggplot(data=picture.two.data,aes(x=genes,y=number,fill=status))+ #数据框 
      		   ylim(0, 30)+
      		   #xlim(0, 50) +
          	   geom_bar(stat="identity", width = 0.5, position = position_dodge(0.9)) + #画柱
      		   scale_fill_manual(values = c("RoyalBlue", "red"), labels=c("Clonal", "Subclonal"))+
      		   geom_text(mapping = aes(label= number))+ #加文本信息number
      		   theme_set(theme_bw())+ #背景白色
      		   theme(panel.grid.major=element_line(colour=NA))+
      		   theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90))+
      		   coord_flip(); #坐标翻转
print(picture)		   
###################################################################################################




 
