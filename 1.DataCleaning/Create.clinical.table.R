#--------------------------------------------------
#               生成临床数据表格                                              
#--------------------------------------------------
#加载结直肠癌克隆突变数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clone.data/colorectal.select.clone.data.RData")
#加载结直肠癌driver gene
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.driver.gene.data/ColorectalDriverGenes.RData")
#加载结直肠癌临床数据
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.clinical.data/colorectal_clinical_data_patient.RData")

#筛选临床数据
intersect.patients <- intersect(colorectal.clinic.data$PATIENT_ID,colorectal.select.table$patient)
colorectal.clinic.data <- colorectal.clinic.data[colorectal.clinic.data$PATIENT_ID %in% intersect.patients,]

#生成刻画临床数据的表格
head(colorectal.clinic.data)
colorectal.clinic.data <- colorectal.clinic.data[,c(1,2,26,6,7,8,9,19,22,23,24,25,27,28,29,30)]
head(colorectal.clinic.data)
colorectal.clinic.data <- colorectal.clinic.data[,c(-13,-14,-15)]
head(colorectal.clinic.data)
table(colorectal.clinic.data$tumour_site_2,useNA = "ifany")
colnames(colorectal.clinic.data) <- c("Patient_ID","sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","OS_status","OS_months","DFS_status","DFS_months","Tumour_site")
colorectal.clinic.data$age <- ifelse(colorectal.clinic.data$age < 67.5,"<67.5",">=67.5")

library(tableone)
library(ReporteRs)
library(rJava)
library(DT)
catVars <- c("sex","MSS_state","T_stage","N_stage","M_stage","TNM_stage","age","Tumour_site")
catTableOverall <- CreateCatTable(vars = catVars, data = colorectal.clinic.data,includeNA = T)
tab_out <- print(catTableOverall)
## Save to a CSV file
write.csv(tab_out, file = "D:/clinical_table.csv")

file <- "D:/clinical_tables.html"
y <- DT::datatable(tab_out,escape = F,rownames=F)
DT::saveWidget(y,file)









