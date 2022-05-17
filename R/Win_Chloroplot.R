# Chloroplot需要以下几种包才能正常运行
# 可通过以下命令安装
# install.packages("BiocManager")
# BiocManager::install("genbankr")
# BiocManager::install("coRdon")
# install.packages("circlize")


rm(list = ls())
args <- commandArgs(T) # args = commandArgs(trailingOnly = TRUE)
# 命令行参数依次为 输出文件夹路径/输入gbk文件路径/样本名
wddir <- format(args[1])
gbkfile <- format(args[2])
filename <- paste(format(args[3]), "circular", sep = ".") # 类似python字符串的join()函数
print(paste("outdir:", wddir))
print(paste("gbk:", gbkfile))
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
PlotPlastidGenome(plot.tables = tmp, file.name = filename, cu.bias = FALSE)
# 去掉密码子偏好性
# Rscript E:/OneDrive/jshy信息部/分发项目/script/chloroplast/R/Chloroplot.R F:/Hibiscus_sabdariffa/ogdraw F:/Hibiscus_sabdariffa/gene_anno/Hibiscus_sabdariffa.gbk Hibiscus_sabdariffa
