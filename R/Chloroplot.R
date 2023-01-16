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
source("/share/nas1/yuj/software/Chloroplot/R/color_complement.R")
source("/share/nas1/yuj/software/Chloroplot/R/converse_ssc.R")
source("/share/nas1/yuj/software/Chloroplot/R/detect_ir.R")
source("/share/nas1/yuj/software/Chloroplot/R/GC_count.R")
source("/share/nas1/yuj/software/Chloroplot/R/gene_color.R")
source("/share/nas1/yuj/software/Chloroplot/R/gene_info.R")
source("/share/nas1/yuj/software/Chloroplot/R/parse_gb_file.R")
source("/share/nas1/yuj/software/Chloroplot/R/plot_genome.R")
source("/share/nas1/yuj/software/Chloroplot/R/read_gb_file.R")
source("/share/nas1/yuj/software/Chloroplot/R/test_parameters.R")
tmp <- PlotTab(gbfile = gbkfile, local.file = T)
if (type_flag == "mt"){
  PlotMitGenome(plot.tables = tmp, file.name = filename, cu.bias = FALSE)
}
if (type_flag == "cp"){
  PlotPlastidGenome(plot.tables = tmp, file.name = filename, cu.bias = FALSE)
}
# 去掉密码子偏好性 cu.bias = FALSE
# Rscript E:/OneDrive/jshy信息部/分发项目/script/chloroplast/R/Chloroplot.R F:/Hibiscus_sabdariffa/ogdraw F:/Hibiscus_sabdariffa/gene_anno/Hibiscus_sabdariffa.gbk Hibiscus_sabdariffa
