#-----------------------------
#volcano plot
#-----------------------------
# ******************画图函数封装路径不要用中文，ggplot2注意print*************************
volcano.Deseq2.plot <- function(dat.plot,
                                gene,
                                filepath){
  suppressMessages(library(ggplot2))
  pdf(filepath,width=10) #"F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP8.volcano.pdf"
  this_tile <- paste0(gene,                    
                      'The number of up gene is ', nrow(dat.plot[dat.plot$lab =='UP',]),                      
                      'The number of down gene is ', nrow(dat.plot[dat.plot$lab =='DOWN',]))
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(dat.plot, aes(x = log2FoldChange, y = -log10(padj), color = lab))+  
              ggtitle(label = this_tile, subtitle = "Colored by fold-change direction") +
              geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
              #geom_text()
              scale_color_manual(name = "Directionality",
                                 values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
              theme_bw(base_size = 14) + # change overall theme
              theme(legend.position = "right") + # change the legend
              xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
              ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
              geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
              scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  print(vol)
  dev.off()
  return(vol)
}

volcano.EdgeR.plot <- function(dat.plot,
                               gene,
                               filepath){
  suppressMessages(library(ggplot2))
  pdf(filepath,width=10) #"F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP8.volcano.pdf"
  this_tile <- paste0(gene,                    
                      'The number of up gene is ', nrow(dat.plot[dat.plot$lab =='UP',]),                      
                      'The number of down gene is ', nrow(dat.plot[dat.plot$lab =='DOWN',]))
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(dat.plot, aes(x = logFC, y = -log10(FDR), color = lab))+  
            ggtitle(label = this_tile, subtitle = "Colored by fold-change direction") +
            geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
            #geom_text()
            scale_color_manual(name = "Directionality",
                               values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
            theme_bw(base_size = 14) + # change overall theme
            theme(legend.position = "right") + # change the legend
            xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
            ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
            geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
            scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  print(vol)
  dev.off()
  return(vol)
}

volcano.limma.voom.plot <- function(dat.plot,
                                    gene,
                                    filepath){
  suppressMessages(library(ggplot2))
  pdf(filepath,width=10) #"F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP8.volcano.pdf"
  this_tile <- paste0(gene,                    
                      'The number of up gene is ', nrow(dat.plot[dat.plot$lab =='UP',]),                      
                      'The number of down gene is ', nrow(dat.plot[dat.plot$lab =='DOWN',]))
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(dat.plot, aes(x = logFC, y = -log10(adj.P.Val), color = lab))+  
				ggtitle(label = this_tile, subtitle = "Colored by fold-change direction") +
				geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
				#geom_text()
				scale_color_manual(name = "Directionality",
								   values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
				theme_bw(base_size = 14) + # change overall theme
				theme(legend.position = "right") + # change the legend
				xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
				ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
				geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
				scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  print(vol)
  dev.off()
  return(vol)
}

box.jitter.plot <- function(dat,
                            filepath){
  # @param dat:可视化数据，Var2和value两个字段
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggsignif))
  pdf(filepath)
  p <- ggplot(data=dat,mapping = aes(x=Var2,y=value,colour=Var2))+
          geom_boxplot()+
          geom_jitter(shape=16, size=1,position=position_jitter(0.2))+
          labs(x="clonality",y="score value",colour="clonality")+
          geom_signif(comparisons = list(unique(as.character(dat$Var2))), test = wilcox.test, step_increase = 0.2)+
          theme_bw()
  print(p)
  dev.off()
  return(p)
}

box.jitter.ggpubr.plot <- function(dat,
                                   filepath){
  # @param dat:可视化数据，Var2和value两个字段  两列数据比较
  suppressMessages(library(ggpubr))
  suppressMessages(library(ggsignif))
  pdf(filepath)
  p <- ggboxplot(dat, x="Var2", y="value", color="Var2", 
            add = "jitter") +
            rotate_x_text(angle = 45) +
            geom_hline(yintercept = mean(dat$value), linetype = 2) + # Add horizontal line at base mean
            #stat_compare_means(label = "p.format", method = "wilcox.test")
            geom_signif(comparisons = list(levels(dat$Var2)[1:2]), test = wilcox.test, step_increase = 0.2)
  print(p)
  dev.off()
  return(p)
}

box.jitter.ggpubr.plot <- function(dat,
                                   filepath){
  # @param dat:可视化数据，Var2和value两个字段
  suppressMessages(library(ggpubr))
  suppressMessages(library(ggsignif))
  pdf(filepath)
  p <- ggboxplot(dat, x="Var2", y="value", color="Var2", 
            add = "jitter") +
            rotate_x_text(angle = 45) +
            geom_hline(yintercept = mean(dat$value), linetype = 2) + # Add horizontal line at base mean
            #stat_compare_means(label = "p.format", method = "wilcox.test")+
            geom_signif(comparisons = list(levels(dat$Var2)[1:2],levels(dat$Var2)[c(1,3)],levels(dat$Var2)[2:3]), test = wilcox.test, step_increase = 0.2)
  print(p)
  dev.off()
  return(p)
}

Venn.plot <- function(input,
                      category=category,
                      fill=c("skyblue", "pink1", "mediumorchid"),
                      filename=filename){
  #' @param input: 一个列表，可视化的数据
  #' @param fill: 设置颜色 e.g. c("skyblue", "pink1", "mediumorchid")
  #' @param catergory: 设置名字  e.g.  c("Deseq2", "EdgeR", "Limma-Voom")
  #' @param filename: 设置存储路径 e.g. "F:/venn.tiff"
        library(VennDiagram)
        #grid.newpage()
        venn.diagram(x=input,
                     fill=fill,
                     alpha=c(0.5,0.5,0.5),cex=2,
                     cat.fontface=4,fontfamily=3,
                     category = category,
                     filename=filename)
}

HeatMap <- function(datMatrix){
  library(pheatmap)
  pheatmap(datMatrix)
}



