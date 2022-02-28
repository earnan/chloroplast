library(optparse)
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="输入文件名[default= fault]", metavar="character"),
    make_option(c("--title"), type="character", default="bar", 
              help="图表的名称[default= fault]", metavar="character"),
    make_option(c("-b", "--bilv"), type="character", default="200", 
              help="比率，丰度在1/200以下的归为others[default= fault]", metavar="character"),
    make_option(c("--outtable"), type="character", default=NULL, 
              help="输出画图数据表格名[default= fault]", metavar="character"),
    make_option(c("--outpdf"), type="character", default=NULL, 
              help="输出PDF名[default= fault]", metavar="character")
  );
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$f)){
  print_help(opt_parser)
  stop("请输入表格文件", call.=FALSE)
  }
if (is.null(opt$outtable)){opt$outtable=paste(opt$f,'.xls',sep='')}
if (is.null(opt$outpdf)){opt$outpdf=paste(opt$f,'.pdf',sep='')}