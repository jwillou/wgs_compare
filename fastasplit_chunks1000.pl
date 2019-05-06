#my $file    = "genome_format_unmask.fa";
#my $outfile = "genome_format_unmask_blocks.fa";
my $size    = 1000;

open (INFILE, "<", $ARGV[0]);
#open (OUTFILE, ">", $ARGV[1]);

my $header;
while(my $line = <INFILE>){
	chomp $line;
	
	#save header information
	if($line=~/^\>/){
		$header = $line;
		next;
	
	#if line is sequence data, split into chunks
	}else{
		my @chunks = $line =~ /[ATCGN]{$size}/g;
		my $i=0;
		foreach (@chunks) {
			$i++;
			
			#header printing can be modified to be easier to parse later if needed
  			print "$header";
  			print "_";
  			print "$i\n";
  			
  			#prints the sequence in whatever length specified by $size above
  			print "$_\n"; #
		}
	}
}
close INFILE;
#close OUTFILE;

exit 0;
