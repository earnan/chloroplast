use File::Basename qw/basename/;


my %id_for_name;

for my $file(@ARGV){
	my $base = basename($file,".gbk");
	open IN,"$file" or die"$!";
	while(<IN>){
		if(/\/organism=\"(.*)\"/){
			$id_for_name{$base} = $1;
			last;
		}
	}
}

for my $id(sort keys %id_for_name){
	print "$id\t$id_for_name{$id}\n";
}

#用法 perl get_id_genbank.pl gbk/*.gbk > id.list