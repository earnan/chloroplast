#/share/nas6/pub/pipline/genome-assembly-seq/chloroplast-genome-seq/v1.2/annotation/bin/chloroplast_annotaion.pl.v2.1.pl

	print STDERR "\n\n\n------------------------------------------------------------------------------\n"
	                  ."---------------------$sample -->Step 07 : OGDRAW ---------------------------------------\n"
	                  ."------------------------------------------------------------------------------\n\n\n";
	##画圈图	多线程下圈图容易错，最后再运行下
	open (CMD,">","$commands_dir/$sample/ogdraw.commands") or die $!;
	{
		$cmd = "Rscript /share/nas1/yuj/software/Chloroplot/Chloroplot.R $outdir/$sample/ogdraw $outdir/$sample/gene_anno/$sample.gbk $sample && convert -density 300 $outdir/$sample/ogdraw/$sample.circular.pdf  $outdir/$sample/ogdraw/$sample.dpi300.circular.png";
		print CMD $cmd,"\n";
	}
	close (CMD) ;
#	$cmd = " sh $commands_dir/$sample/ogdraw.commands ";
#	&PIPLINE($cmd,"ogdraw",$sample);




#/share/nas6/pub/pipline/genome-assembly-seq/chloroplast-genome-seq/v1.2/html_report/report2xml.yelvti.pl

#-----4.3.2 叶绿体基因组图谱
&EMPTY_TAG('h3','4.3.2 叶绿体基因组图谱','','type1');
&EMPTY_TAG('p','使用Chloroplot（https://irscope.shinyapps.io/Chloroplot/）制作叶绿体基因组图谱，如下图所示：','type1');

if(glob "$indir/04Annotation/*/ogdraw/*dpi300.circular.png"){
	&piclist("图$pid 叶绿体基因组图谱",'注：正向编码的基因位于圈内侧，反向编码的基因位于圈外侧。内部的灰色圈代表GC含量。',"$indir/04Annotation/*/ogdraw/*dpi300.circular.png");	
}else{
	&piclist("图$pid 叶绿体基因组图谱",'注：正向编码的基因位于圈内侧，反向编码的基因位于圈外侧。内部的灰色圈代表GC含量。',"$indir/04Annotation/*/ogdraw/*.circular.png");
}

