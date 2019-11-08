#'-----------------------------------------------------
#'           统计表达矩阵，突变矩阵，表型矩阵
#'                  病人的交集情况
#'-----------------------------------------------------

load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData")
load(file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/colorectal.expression.ENSG.data.RData")

library(VennDiagram)
library(grDevices)
Group_clin <- unique(colorectal.clinic.data$PATIENT_ID) # 629
Group_mutation <- unique(colorectal.select.table$patient) # 536
colnames(temp.COADREAD.ReadCount.matrix) <- substring(colnames(temp.COADREAD.ReadCount.matrix),1,12)
Group_Exp <- unique(colnames(temp.COADREAD.ReadCount.matrix)) # 532
input <- list(Group_clin,Group_mutation,Group_Exp)
venn.diagram(x=input,
                     fill=c("skyblue", "pink1","goldenrod1"),
                     alpha=c(0.5,0.5,0.5),cex = 2.5,
                     cat.cex = 2.5,
                     lwd = c(2,2,2),
                     cat.fontface=2,
                     category = c("clinic", "mutation" , "Expr"),filename="D:/clinic.mutation.exp.tiff")

