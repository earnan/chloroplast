# Chloroplot需要以下几种包才能正常运行
# 可通过以下命令安装
# install.packages("BiocManager")
# BiocManager::install("genbankr")
# BiocManager::install("coRdon")
# install.packages("circlize")


rm(list = ls())
args <- commandArgs(T) # args = commandArgs(trailingOnly = TRUE)
# 命令行参数依次为 输出文件夹路径/输入gbk文件路径/样本名/gbk类型(cp or mt)
wddir <- format(args[1])
gbkfile <- format(args[2])
filename <- paste(format(args[3]), "circular", sep = ".") # 类似python字符串的join()函数
type_flag <- format(args[4]) # cp/mt
print(paste("gbk:", gbkfile))
if (type_flag == "mt"){
  print("genome type: Mitochondrion")
}
if (type_flag == "cp"){
  print("genome type: Chloroplot")
}
print(paste("outdir:", wddir))
print(paste("outfile:", filename))

setwd(wddir)
library(dplyr)
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/color_complement.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/converse_ssc.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/detect_ir.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/GC_count.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/gene_color.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/gene_info.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/parse_gb_file.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/plot_genome.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/read_gb_file.R")
source("D:/OneDrive - cancer/Bioinfo_analysis/Chloroplot/R/test_parameters.R")
tmp <- PlotTab(gbfile = gbkfile, local.file = T)
if (type_flag == "mt"){
  PlotMitGenome(plot.tables = tmp, file.name = filename, cu.bias = FALSE)
}
if (type_flag == "cp"){
  PlotPlastidGenome(plot.tables = tmp, file.name = filename, cu.bias = FALSE)
}
# 去掉密码子偏好性 cu.bias = FALSE
# Rscript E:\OneDrive\jshy信息部\Script\chloroplast\R\Win_Chloroplot.R F:\ F:\Ustilago_esculenta_MT10.gbk Ustilago_esculenta_MT10 mt
# Rscript E:\OneDrive\jshy信息部\Script\chloroplast\R\Win_Chloroplot.R F:\ E:\Server\3171\complete_dir\04Annotation\Lithocarpus_litseifolius_Hance_Chun\Lithocarpus_litseifolius_Hance_Chun.gbk Lithocarpus_litseifolius_Hance_Chun cp