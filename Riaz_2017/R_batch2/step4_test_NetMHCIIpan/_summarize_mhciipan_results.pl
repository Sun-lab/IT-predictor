$finput  = $ARGV[0];
$foutput = $ARGV[1];
$peptide_length = 9;

open(FIN, $finput) or die "Cannot find file $finput \n";
open(FOUT, ">", $foutput) or die "Cannot create file $foutput\n";
print FOUT "seq\tallele\tpeptide\tof\tcore\tcore_rel\tidentity\tscore_EL\tpercent_rank_EL\texp_bind\tbind_level\n";

# the header of each result table includes 12 columns

# Seq Allele Peptide Identity Pos Core Core_Rel 1-log50k(aff) Affinity(nM) %Rank 
# Exp_Bind BindingLevel

while(<FIN>){
	chomp;
	s/^\s+//;
	
	@items = split(/\s+/, $_);

	if($items[0] eq "Pos" && $items[1] eq "MHC"){
	  # print "@items\n";
		$line = <FIN>;
		chomp($line);
		$line = ~ s/^\s+//;
		
		if(! $line =~ /^---------------------------------------/){
			die "expected a line of '----...' aftr the headers\n";
		}
		
		$counter  = 1;
		
		while(<FIN>){
			chomp;
	    s/^\s+//;

		  if($_ =~ /^---------------------------------------/){ last; }
		  
			@items = split(/\s+/);
			
			if($items[0] != $counter){ print "@items\n"; die "unexpected value for column 'Pos'\n"; }
			
			if($counter == 1){
				$hla      = $items[1];
				$identity = $items[6];
			}else{
			  if($hla      != $items[1]){ die "HLA do not match\n"; }
			  if($identity != $items[6]){ 
			    print "identity is $identity, and items[3] is $items[3]\n";
			    die "identity does not match\n"; 
			  }
			}
			
			print FOUT "$items[0]\t$items[1]\t$items[2]\t$items[3]\t$items[4]\t$items[5]\t$items[6]\t$items[8]\t$items[9]\t$items[10]\n";
			$counter += 1;
		}
		
	}

}
close(FIN);
close(FOUT);
