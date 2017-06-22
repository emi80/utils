#! /usr/bin/perl
use strict; use warnings; 
use Getopt::Long;

#CCK & AB <15/may/2015>

# default values
my $asta;	# option variable with default value
my $ssj;	
my $ssc;	
my $event;
my $threshold = 10;
my $out="astalavista.psi.out";
my $verbose;
my $help;

if (!@ARGV) {
    print "$0: Argument required.\n";
    usage();
}


GetOptions(
    'asta|a=s' => \$asta,
    'ssj|b=s' => \$ssj,
    'ssc|c:s' => \$ssc,
    'event|e=s' => \$event,
    'Threshold|t:i' => \$threshold,
    'out|o:s' => \$out,
    'help|h|?' => \$help
   ) or die("Error in command line arguments\n\n*** --help *** for usage\n\n");

if($help){
    usage();
}

sub usage {
    print STDERR<<USAGE;

  Description: compute PSI values of pairwise AStalavista AS events based on ipsa splice site 
               junctions and splice site counts (boundary).
               Output: the AStalavista file with the PSI as the last field and a TSV file id event and PSI
	               in one folder per event.

  Usage: perl $0
    --asta|a      : pairwise AS events from astalavista output gtf file (1-based coordinate)
    --ssj|b       : splice site junction quantifications from ipsa pipeline BED file (1-based coordinate)
    --ssc|c       : splice site boundary quantifications from ipsa pipeline BED file (1-based coordinate)
                    only required for intron retention
    --event|e     : list of events to be processed (comma-separated)
                    ESS = exon skipping single
		    ESM = exon skipping multiple
	            IR  = intron retention
		    ME  = mutually exclusive exons
		    AD  = alternative donor
		    AA  = alternative acceptor
    --threshold|t : minimum number of read counts for the sum of all splice site junctions and boundaries
                    PSI=NA if the sum is below the threshold  (default = 10) 
		    e.g. single exon skipping: PSI=a+b/a+b+2c, a+b+c > threshold
    --out|o       : output file name [optional] (without extension)
    --outdir      : folder to write output files
#    --verbose|v   : foreach event print type_of_event structure PSI
    --help
    
  Example:
    perl $0 --asta astalavista.gtf --ssj sample.ssj.bed --ssc sample.ssc.bed --event ESS,IR --outdir

  Information on PSI calculation for the different AS events:
    display ~cklein/utils/usage/psi.png

USAGE
    exit();
}
###########################################################################################
print STDERR "#============= PSI LOCAL: ASTALAVISTA + IPSA ==================\n\n";

#=====================================================================
# warn "Splice chain from second transcript has two equal coordinates"
#=====================================================================
my $warn=0;

#========
# events
#========
my @events= split(/,/,$event);
print STDERR "SELECTED EVENTS:\n";
foreach my $e (@events){
    print STDERR "$e\n";
    `mkdir -p $e`;
}

#=============================
# output files for each event
#=============================
my $essg;
my $esst;
my $irg;
my $irt;
my $meg;
my $met;
my $adg;
my $adt;
my $aag;
my $aat;
my $esmg;
my $esmt;

if("ESS" ~~ @events){
    open($essg,">ESS/$out.gtf"); #open output file
    open($esst,">ESS/$out.tsv"); #open output file
}
if("IR" ~~ @events){
    open($irg,">IR/$out.gtf"); #open output file
    open($irt,">IR/$out.tsv"); #open output file
}
if("ME" ~~ @events){
    open($meg,">ME/$out.gtf"); #open output file
    open($met,">ME/$out.tsv"); #open output file
}
if("AD" ~~ @events){
    open($adg,">AD/$out.gtf"); #open output file
    open($adt,">AD/$out.tsv"); #open output file
}
if("AA" ~~ @events){
    open($aag,">AA/$out.gtf"); #open output file
    open($aat,">AA/$out.tsv"); #open output file
}
if("ESM" ~~ @events){
    open($esmg,">ESM/$out.gtf"); #open output file
    open($esmt,">ESM/$out.tsv"); #open output file
}


#===================
# READING SSJ FILE
#===================
my %ssj;
my $cSSJ=0;
print STDERR "\nREADING SSJ INPUT FILE: $ssj\n";
open(J,"<$ssj" || die "can't open file $ssj: $!");
while(my $l = <J>){
    chomp($l);
    my @tmp = split(/\t/,$l);
	$ssj{$tmp[0]}=$tmp[1];
	$cSSJ=$cSSJ+1;
}
close(J);
print STDERR "$cSSJ input lines\n"; 

#======================================
# READING SSC FILE IF NEEDED (only IR)
#======================================
 my %ssc;
if("IR" ~~ @events){
    if(defined $ssc){
	my $cSSC=0;
	print STDERR "\nREADING SSC INPUT FILE: $ssc\n";
	open(C,"<$ssc" || die "can't open file $ssc: $!");
	while(my $l = <C>){
	    chomp($l);
	    my @tmp = split(/\t/,$l);
	    $ssc{$tmp[0]}=$tmp[1];
	    $cSSC=$cSSC+1;
	}
	close(C);
	print STDERR "$cSSC input lines\n"; 
    }else{ 
	die("\nEXITING: SSC file is required for IR.\n\n\n");
    }
}
#=================================================
# READING AStalavista GTF and compute PSI values
#=================================================
print STDERR "\nREADING AStalavista INPUT FILE: $asta\n";

open(A,"<$asta" || die "can't open file $asta: $!");
while(my $line = <A>){
    chomp($line);
    my @tmp = split(/\t|\"/,$line);
    my $chr=$tmp[0];
    my $l = @tmp;
    my $id=$tmp[$l-2];

    #check if flanks are genomic sorted
    my $flank; 
    if(($tmp[3] =~ /^[0-9]+$/ ) and ($tmp[4] =~ /^[0-9]+$/ )){
	if($tmp[3] < $tmp[4]){
	    $flank = "$tmp[3]_$tmp[4]";
	}else{
	    $flank = "$tmp[4]_$tmp[3]";
	    print STDERR "Inverted flanks $tmp[4]_$tmp[3]\n";
	}
    }else{
	$flank = "$tmp[3]_$tmp[4]";
	print STDERR "Flanks $tmp[3]_$tmp[4] are not numeric\n";
    }

    my $strand = $tmp[6];
    my $structure = $tmp[15]; 
    my $splice_chain = $tmp[17]; 


    #-------------------------
    # EXON SKIPPING SINGLE
    #-------------------------
    if (($structure =~ /^0\,1\-2\^$/) && ("ESS" ~~ @events)) {
	# psi=(a+b)/(a+b+2c)
	# only need ssj
	# complete:psi=1
	# skipped:psi=0
	# a=f1,s1;b=s2,f2; (sorted: always from sc2)
	# c=f1,f2;

	my ($f1,$f2) = split(/\_/,$flank); #flanking 1 and 2
	
	# confirm that the splice chain has the expected format
	if($splice_chain =~ /^,\d+\-\d+\^$/){
	
	    my ($sc1,$sc2) = split(/\,/,$splice_chain); #splice_chain from 1st/2nd transcripts
	    
	    # first transcript(s) should be always empty, i.e. it is the exon skipped transcript
	    # my @sc1 = split(/\^|\-|\[|\]/,$sc1); #s1 is $sc1[0] and s2 is $sc1[1] (from the 1st transcripts)
	    
	    my @tmpSC2 = split(/\^|\-|\[|\]/,$sc2); #s1 is $sc2[0] and s2 is $sc2[1] (from the 2nd transcripts)
	    my @sc2;
	    # second transcript(s)
	    my $flag_compute_psi=0;
	    if($tmpSC2[0] == $tmpSC2[1]){$warn=$warn+1;}
	    elsif(($strand eq "+") and ($tmpSC2[0] < $tmpSC2[1])){
		$sc2[0]=$tmpSC2[0];
		$sc2[1]=$tmpSC2[1];
		$flag_compute_psi=1;
	    }elsif(($strand eq "-") and ($tmpSC2[1] < $tmpSC2[0])){
		$sc2[0]=$tmpSC2[1];
		$sc2[1]=$tmpSC2[0];
		$flag_compute_psi=1;
	    }else{
		print STDERR "\n**********\nStrand $strand not available or it does not match the expected genomic position\n$line\n****\n";
	    }
	    
	    # compute PSI
	    my $psi;
	    if($flag_compute_psi){
		my $a = 0;
		my $b = 0;
		my $c = 0;
		
		if(defined $ssj{$chr."_".$f1."_".$sc2[0]."_".$strand}){
		    $a = $ssj{$chr."_".$f1."_".$sc2[0]."_".$strand};
		}
		if(defined $ssj{$chr."_".$sc2[1]."_".$f2."_".$strand}){
		    $b = $ssj{$chr."_".$sc2[1]."_".$f2."_".$strand};
		}
		if(defined $ssj{$chr."_".$flank."_".$strand}){
		    $c = $ssj{$chr."_".$flank."_".$strand}; 
		}
		
		
	    
		if( ($a+$b+$c) > $threshold ){
		    $psi = (($a + $b)/($a + $b + 2*$c));	
		}
		else{
		    $psi="NA";
		}
	    }else{
		$psi="NA";
	    }

	    # print gtf output
	    print $essg $line."psi \"";
	    if($psi ne "NA"){
		printf $essg ("%.4f", $psi);
	    }else{
		print $essg $psi;
	    }
	    print $essg "\";\n";
	    
            # print tsv output
	    print $esst "$id\t";
	    if($psi ne "NA"){
		printf $esst ("%.4f", $psi);
	    }else{
		print $esst $psi;
	    }
	    print $esst "\n";
	    
	}else{
	    print STDERR "Splice chain $splice_chain does not have the expected format $structure\n";
	}
    }
    

    #------------------
    # INTRON RETENTION
    #------------------
    elsif($structure =~/^0\,1\^2\-$/ && ("IR" ~~ @events)){
	# psi=(2a)/(b+c+2a)
	# no intron:psi=1
	# retained:psi=0
	# requires ssj (a) and ssc (b and c)
	# a=s1,s2; (sorted: always from sc2)
	# b=s1; c=s2;
	
	my ($f1,$f2) = split(/\_/,$flank); #flanking 1 and 2
	my ($sc1,$sc2) = split(/\,/,$splice_chain); #splice_chain from 1st/2nd transcripts

	# confirm that the splice chain has the expected format
	if($splice_chain =~ /^,\d+\^\d+\-$/){
	    
	    # first transcript(s) should be always empty, i.e. it is the intron retained transcript
	    # my @sc1 = split(/\^|\-|\[|\]/,$sc1); #s1 is $sc1[0] and s2 is $sc1[1] (from the 1st transcripts)
	    
	    # second transcript(s)
	    my @tmpSC2 = split(/\^|\-|\[|\]/,$sc2); #s1 is $sc2[0] and s2 is $sc2[1] (from the 2nd transcripts)
	    my @sc2;
	    if(($strand eq "+") and ($tmpSC2[0] < $tmpSC2[1])){
		$sc2[0]=$tmpSC2[0];
		$sc2[1]=$tmpSC2[1];
	    }elsif(($strand eq "-") and ($tmpSC2[1] < $tmpSC2[0])){
		$sc2[0]=$tmpSC2[1];
		$sc2[1]=$tmpSC2[0];
	    }else{
		print STDERR "\n**********\nStrand $strand not available or it does not match the expected genomic position\n$line\n****\n";
	    }
	    
	    my $a = 0;
	    my $b = 0;
	    my $c = 0;
	    if(defined $ssj{$chr."_".$sc2[0]."_".$sc2[1]."_".$strand}){
		$a = $ssj{$chr."_".$sc2[0]."_".$sc2[1]."_".$strand};
	    }
	    if(defined $ssc{$chr."_".$sc2[0]."_".$strand}){
		$b = $ssc{$chr."_".$sc2[0]."_".$strand};
	    }
	    if(defined $ssc{$chr."_".$sc2[1]."_".$strand}){
		$c = $ssc{$chr."_".$sc2[1]."_".$strand};
	    }
	    

	    # compute PSI
	    my $psi;
	    
	    if(($a+$b+$c) > $threshold ){
		$psi = ((2*$a)/($b + $c + 2*$a));		
	    }else{
		$psi="NA";
	    }
	    
	    # print gtf output
	    print $irg $line."psi \"";
	    if($psi ne "NA"){
		printf $irg ("%.4f", $psi);
	    }else{
		print $irg $psi;
	    }
	    print $irg "\";\n";
    
	    # print tsv output
	    print $irt "$id\t";
	    if($psi ne "NA"){
		printf $irt ("%.4f", $psi);
	    }else{
		print $irt $psi;
	    }
	    print $irt "\n";
	    
	}else{
	    print STDERR "Splice chain $splice_chain does not have the expected format $structure\n";
	}
    } 
    

    #-------------------
    # MUTUALLY EXCLUSIVE
    #-------------------
    elsif($structure=~/^1\-2\^\,3\-4\^$/ && ("ME" ~~ @events) ){
	# psi=(c+d)/(c+d+a+b)
	# requires only ssj
	# c=f1,s1;d=s2,f2; (sorted: always from sc1)
	# a=f1,s3; b=s4,f2;(sorted: always from sc2)

	
	# confirm that the splice chain has the expected format
	if($splice_chain =~ /^\d+\-\d+\^,\d+\-\d+\^$/){
	    my ($f1,$f2) = split(/\_/,$flank); #flanking 1 and 2
	    my ($sc1,$sc2) = split(/\,/,$splice_chain); #splice_chain from 1st/2nd transcripts
	    
	    # first transcript(s)
	    my @tmpSC1 = split(/\^|\-|\[|\]/,$sc1); #s1 is $sc1[0] and s2 is $sc1[1] (from the 1st transcripts)
	    # second transcript(s)
	    my @tmpSC2 = split(/\^|\-|\[|\]/,$sc2); #s1 is $sc2[0] and s2 is $sc2[1] (from the 2nd transcripts)
	    
   	    my @sc1;
	    my @sc2;
	    my $flag=0; # true if conditions satisfied
	    
	    if(($strand eq "+") and ($tmpSC1[0] < $tmpSC1[1]) and ($tmpSC2[0] < $tmpSC2[1])){
		$sc1[0]=$tmpSC1[0];
		$sc1[1]=$tmpSC1[1];
		
		$sc2[0]=$tmpSC2[0];
		$sc2[1]=$tmpSC2[1];
		
		$flag = 1;
		
	    }
	    # invert sc1 and sc2
	    elsif(($strand eq "-") and ($tmpSC1[1] < $tmpSC1[0]) and ($tmpSC2[1] < $tmpSC2[0])){
		
		if($tmpSC2[0] < $tmpSC1[0]){
		    $sc1[0]=$tmpSC2[1];
		    $sc1[1]=$tmpSC2[0];
		    
		    $sc2[0]=$tmpSC1[1];
		    $sc2[1]=$tmpSC1[0];
		    
		    $flag =1;
		    
		}else{
		    print STDERR "Sorting Strand $strand to genomic position: transcripts are not inverted\n";
		}
		
	    }else{
		print STDERR "\n**********\nStrand $strand not available or it does not match the expected genomic position\n$line\n****\n";
	    }
	    
	    if($flag){
		my $c = 0;
		my $d = 0;
		my $a = 0;
		my $b = 0;
		if(defined $ssj{$chr."_".$f1."_".$sc1[0]."_".$strand}){
		    $c = $ssj{$chr."_".$f1."_".$sc1[0]."_".$strand};
		}
		if(defined $ssj{$chr."_".$sc1[1]."_".$f2."_".$strand}){
		    $d = $ssj{$chr."_".$sc1[1]."_".$f2."_".$strand};
		}
		if(defined $ssj{$chr."_".$f1."_".$sc2[0]."_".$strand}){
		    $a = $ssj{$chr."_".$f1."_".$sc2[0]."_".$strand};
		}
		if(defined $ssj{$chr."_".$sc2[1]."_".$f2."_".$strand}){
		    $b = $ssj{$chr."_".$sc2[1]."_".$f2."_".$strand};
		}
		
		# compute PSI
		my $psi;
		
		if(($a+$b+$c+$d) > $threshold ){
		    $psi = (($c + $d)/($c + $d + $a+ $b));		
		}else{
		    $psi="NA";
		}

		# print gtf output
		print $meg $line."psi \"";
		if($psi ne "NA"){
		    printf $meg ("%.4f", $psi);
		}else{
		    print $meg $psi;
		}
		print $meg "\";\n";
		
		# print tsv output
		print $met "$id\t";
		if($psi ne "NA"){
		    printf $met ("%.4f", $psi);
		}else{
		    print $met $psi;
		}
		print $met "\n";
		
	    }
	}else{
	    print STDERR "Splice chain $splice_chain does not have the expected format $structure\n";
	}
    }

    #-------------------
    # ALTERNATIVE DONOR
    #-------------------
    elsif($structure=~/^1\^\,2\^$/  && ("AD" ~~ @events)){
	# psi=(a)/(a+b)
	# requires only ssj
	# a=s1,f2; (sorted: always from sc1)
	# b=s2,f2; (sorted: always from sc2)
	
	#check if the splice chain has the expected format
	if($splice_chain=~/^(\d+)\^\,(\d+)\^$/){
	    my $tmp1 = $1; #splice_chain from 1st transcript
	    my $tmp2 = $2; #splice_chain from 2nd transcript
	    my $sc1;
	    my $sc2;
	    
	    if(($strand eq "+") and ($tmp1 < $tmp2)){
		$sc1=$tmp1;
		$sc2=$tmp2;
	    }elsif(($strand eq "-") and ($tmp1 > $tmp2)){
		$sc1=$tmp2;
		$sc2=$tmp1;
	    }else{
		print STDERR "\n**********\nStrand $strand not available or it does not match the expected genomic position\n$line\n****\n";
	    }
	    
	    my ($f1,$f2) = split(/\_/,$flank); #flanking 1 and 2
	    
	    my $a = 0;
	    my $b = 0;
	    if(defined $ssj{$chr."_".$sc1."_".$f2."_".$strand}){
		$a = $ssj{$chr."_".$sc1."_".$f2."_".$strand};
	    }
	    if(defined $ssj{$chr."_".$sc2."_".$f2."_".$strand}){
		$b = $ssj{$chr."_".$sc2."_".$f2."_".$strand};
	    }
	    
	    # compute PSI
	    my $psi;
	    
	    if( ($a+$b) > $threshold ){
		$psi = (($a)/($a+ $b));		
	    }else{
		$psi="NA";
	    }
	    # print gtf output
	    print $adg $line."psi \"";
	    if($psi ne "NA"){
		printf $adg ("%.4f", $psi);
	    }else{
		print $adg $psi;
	    }
	    print $adg "\";\n";
	    
	    # print tsv output
	    print $adt "$id\t";
	    if($psi ne "NA"){
		printf $adt ("%.4f", $psi);
	    }else{
		print $adt $psi;
	    }
	    print $adt "\n";
	    	    
	}
       
    }

    #---------------------------
    # ALTERNATIVE ACCEPTOR
    #---------------------------
    elsif($structure=~/^1\-\,2\-$/  && ("AA" ~~ @events)){
	# psi=(a)/(a+b)
	# requires only ssj
	# a=f1,s1; (sorted: always from sc1)
	# b=f1,s2; (sorted: always from sc2)	
	
	#check if the splice chain has the expected format
	if($splice_chain=~/^(\d+)\-\,(\d+)\-$/){
	    my $tmp1 = $1; #splice_chain from 1st transcript
	    my $tmp2 = $2; #splice_chain from 2nd transcript
	    my $sc1;
	    my $sc2;

	    if(($strand eq "+") and ($tmp1 < $tmp2)){
		$sc1=$tmp1;
		$sc2=$tmp2;
	    }elsif(($strand eq "-") and ($tmp1 > $tmp2)){
		$sc1=$tmp2;
		$sc2=$tmp1;
	    }else{
		print STDERR "\n**********\nStrand $strand not available or it does not match the expected genomic position\n$line\n****\n";
	    }   
	    
	    my ($f1,$f2) = split(/\_/,$flank); #flanking 1 and 2
	    
	    my $a = 0;
	    my $b = 0;
	    if(defined $ssj{$chr."_".$f1."_".$sc1."_".$strand}){
		$a = $ssj{$chr."_".$f1."_".$sc1."_".$strand};
	    }
	    if(defined $ssj{$chr."_".$f1."_".$sc2."_".$strand}){
		$b = $ssj{$chr."_".$f1."_".$sc2."_".$strand};
	    }
	    
	    # compute PSI
	    my $psi;
	    
	    if( ($a+$b) > $threshold ){
		$psi = (($a)/($a + $b));		
	    }else{
		$psi="NA";
	    }
	    # print gtf output
	    print $aag $line."psi \"";
	    if($psi ne "NA"){
		printf $aag ("%.4f", $psi);
	    }else{
		print $aag $psi;
	    }
	    print $aag "\";\n";
	    
	    # compute tsv output
	    print $aat "$id\t";
	    if($psi ne "NA"){
		printf $aat ("%.4f", $psi);
	    }else{
		print $aat $psi;
	    }
	    print $aat "\n";
	    
	}
	
    }

    #-----------------------
    # EXON SKIPPING MULTIPLE
    #-----------------------
    elsif($structure=~/^0\,1\-2\^3\-4\^/ && ("ESM" ~~ @events)){
	# pattern matching at least two exons skipped
	# number of junctions is precisely the number of splice sites (n) 
	# (example with 2 skipped exons => 4 new splice sites s1..s4, 
	# 3 junctions (a,b,c) with the alternative aplice sites ie a=f1,s1;s2,s3;s4,f2
	# and one when exons are skipped e=f1,f2)
	# psi=(a+b+c+d+..)/(a+b+c+d+..+(n-1)*e)
	# requires only ssj
	# a=f1,s1;b=s2,s3;c=s4,f2; (sorted: always from sc2)
	# e=f1,f2;

    
	#check whether it is only MULTI exon skipping and not a more complex pattern
	my @matches = ($structure =~ /(\d+\-\d+\^)/g);
	my $l = @matches; #number of skipped exons
	my $n = 2*$l; # number of splice sites
	my $string = $structure;
	$string =~ s/^0\,//; #remove "0," from the beginning of the sequence
	$string =~ s/(\d+\-\d+\^)//g; #remove as many times as the pattern appears
	
	if($string eq ""){ #if nothing is left, it is only MULTI exon skipping

	    my ($f1,$f2) = split(/\_/,$flank); #flanking 1 and 2
	    my ($sc1,$sc2) = split(/\,/,$splice_chain); #splice_chain from 1st/2nd transcripts
	    # first transcript is empty, i.e. multi exon skipping
	    if($sc1 eq ""){
		
		my @tmpSC2 = split(/\^|\-|\[|\]/,$sc2); #s1 is $sc2[0] and s2 is $sc2[1] (from the 2nd transcripts)
		my @sc2;
		
		if($strand eq "+"){
		    @sc2 = @tmpSC2;
		}else{
		    @sc2 = reverse @tmpSC2;
		}
		
		my $psi;
		my $numeratorPSI = 0;
		my $denominatorPSI = 0; #will be the same as the numerator + the last value (n-1)*ssj{f1_f2}
		my $sum = 0; #test for the threshold: a+b+..+n i.e. numeratorPSI + the last term only for denominator withou "(n-1)*"
		
		#first junction
		if(defined $ssj{$chr."_".$f1."_".$sc2[0]."_".$strand}){
		    $numeratorPSI = $ssj{$chr."_".$f1."_".$sc2[0]."_".$strand}; #first junction
		}
		
		for(my $i=1; $i<($n-2); $i++){#interior juntions, numerator and denominator get the same values
		    if(defined  $ssj{$chr."_".$sc2[$i]."_".$sc2[$i+1]."_".$strand}){
			$numeratorPSI = $numeratorPSI + $ssj{$chr."_".$sc2[$i]."_".$sc2[$i+1]."_".$strand};
		    }
		}
		#last junction(s_n,f2)
		if(defined $ssj{$chr."_".$sc2[$n-1]."_".$f2."_".$strand}){
		    $numeratorPSI = $numeratorPSI + $ssj{$chr."_".$sc2[$n-1]."_".$f2."_".$strand};
		}
		#last term for the denominator only (n-1)*juntion(f1,f2)
		$denominatorPSI = $numeratorPSI;
		$sum = $numeratorPSI;
		if(defined $ssj{$chr."_".$f1."_".$f2."_".$strand}){
		    $denominatorPSI = $denominatorPSI + (($n-1)*($ssj{$chr."_".$f1."_".$f2."_".$strand}));
		    $sum = $sum + $ssj{$chr."_".$f1."_".$f2."_".$strand};
		}
		if($sum > $threshold){
		    $psi = $numeratorPSI / $denominatorPSI; 
		}else{
		    $psi="NA";
		}
		
		# print gtf output
		print $esmg $line."psi \"";
		if($psi ne "NA"){
		    printf $esmg ("%.4f", $psi);
		}else{
		    print $esmg $psi;
		}
		print $esmg "\";\n";
		
		# print tsv output
		print $esmt "$id\t";
		if($psi ne "NA"){
		    printf $esmt ("%.4f", $psi);
		}else{
		    print $esmt $psi;
		    }
		print $esmt "\n";
		
	    }
	}
    } # close the last case of structure, otherwise complex AS/DSP/VST event and PSI=NA
}

# close files
close(A);
if("ESS" ~~ @events){
    close($essg);
    close($esst);
}
if("IR" ~~ @events){
    close($irg);
    close($irt);
}
if("ME" ~~ @events){
    close($meg);
    close($met);
}
if("AD" ~~ @events){
    close($adg);
    close($adt);
}
if("AA" ~~ @events){
    close($aag);
    close($aat);
}
if("ESM" ~~ @events){
    close($esmg);
    close($esmt);
}

# warnings
print STDERR "\nWARNING: $warn Splice chains from second transcript have two equal coordinates\n\n\n";
