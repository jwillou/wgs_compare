my $file = "genome.fa";
my $outfile = "genome_format.fa";
open (INFILE, "<", $file);
open (OUTFILE, ">", $outfile);

my $i = 0; #counter for first sequence
my $seq;
while(my $line = <INFILE>){
	$i++;
	chomp $line;
	
	#print sequence and new header
	if($line=~/^\>/){
		if($i > 1){
			print OUTFILE "$seq\n";
		}
		print OUTFILE "$line\n";
		$seq = undef;
		next;
	
	#if line is sequence data concatonate all seq lines
	}else{
		if($seq){
			$seq = "$seq$line";
			next;
		}
		else{
			$seq = $line;
			next;
		}
	}
}
print OUTFILE "$seq\n";

close INFILE;
close OUTFILE;

exit 0;