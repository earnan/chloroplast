	print STDERR "\n\n\n------------------------------------------------------------------------------\n"
	                  ."---------------------$sample -->Step 07 : OGDRAW ---------------------------------------\n"
	                  ."------------------------------------------------------------------------------\n\n\n";
	##画圈图	多线程下圈图容易错，最后再运行下
	open (CMD,">","$commands_dir/$sample/ogdraw.commands") or die $!;
	{
		$cmd = "Rscript /share/nas1/yuj/software/Chloroplot/Chloroplot.R $outdir/$sample/ogdraw $outdir/$sample/gene_anno/$sample.gbk $sample";
		
		print CMD $cmd,"\n";
	}
	close (CMD) ;
#	$cmd = " sh $commands_dir/$sample/ogdraw.commands ";
#	&PIPLINE($cmd,"ogdraw",$sample);
