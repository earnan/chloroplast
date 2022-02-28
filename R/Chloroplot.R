# Chloroplot需要以下几种包才能正常运行
# 可通过以下命令安装
# install.packages("BiocManager")
# BiocManager::install("genbankr")
# BiocManager::install("coRdon")
# install.packages("circlize")



rm(list = ls())
args <- commandArgs(T) # args = commandArgs(trailingOnly = TRUE)
wddir = format(args[1])
gbkfile = format(args[2])
print(wddir)
print(gbkfile)
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
PlotPlastidGenome(plot.tables = tmp)
# Rscript E:/OneDrive/jshy信息部/分发项目/script/chloroplast/R/Chloroplot.R F:/test/222 F:/test/111/Hibiscus_sabdariffa.gbk
