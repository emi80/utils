#! /usr/bin/perl
use strict; use warnings; 
use Getopt::Long;

#CCK & AB <15/may/2015>

# default values
my $asta;	# option variable with default value
my $ssj;	# option variable with default value
my $ssc;	# option variable with default value
my $threshold = 10;
my $out="astalavista.psi.out.gtf";
my $verbose;
my $help;

GetOptions(
    'asta|a=s' => \$asta,
    'ssj|b=s' => \$ssj,
    'ssc|c=s' => \$ssc,
    'Threshold|t:i' => \$threshold,
    'out|o:s' => \$out,
    'verbose|v' => \$verbose,
    'help|h|?' => \$help
   ) or die("Error in command line arguments\n\n*** --help *** for usage\n\n");

if($help){
    usage();
}

sub usage {
    print STDERR<<USAGE;

  Description: compute PSI values of pairwise AStalavista AS events based on ipsa splice site 
               junctions and splice site counts (boundary).
               Output: the AStalavista file with the PSI as the last field.

  Usage: perl $0
    --asta|a      : pairwise AS events from astalavista output gtf file (1-based coordinate)
    --ssj|b       : splice site junction quantifications from ipsa pipeline BED file (1-based coordinate)
    --ssc|c       : splice site boundary quantifications from ipsa pipeline BED file (1-based coordinate)
    --threshold|t : minimum number of read counts for the sum of all splice site junctions and boundaries
                    PSI=NA if the sum is below the threshold  (default = 10) 
		    e.g. single exon skipping: PSI=a+b/a+b+2c, a+b+c > threshold
    --out|o       : output file name [optional] (default "astalavista.psi.out.gtf")
    --verbose|v   : foreach event print type_of_event structure PSI
    --help
    
  Example:
    perl $0 --asta astalavista.gtf --ssj sample.ssj.bed --ssc sample.ssc.bed

  Information on PSI calculation for the different AS events:
    display ~cklein/utils/usage/psi.jpg

USAGE
    exit();
}
###########################################################################################


# ssj bed
my %ssj;
my $cSSJ=0;
print STDERR "\nReading ssj input file: $ssj\n";
open(J,"<$ssj" || die "can't open file $ssj: $!");
while(my $l = <J>){
    chomp($l);
    my @tmp = split(/\t/,$l);
    if($tmp[0] =~ /^[^\_]+\_(\d+)\_(\d+)\_[+|-]$/){
	my $k = $1."_".$2; #key=coordinate_coordinate; value=counts;
	$ssj{$k}=$tmp[1];
	$cSSJ=$cSSJ+1;
    }else{
	print STDERR "Problem reading ssj file $ssj\n";
    }
}
close(J);
print STDERR "$cSSJ input lines from $ssj\n"; 


# ssc bed
my %ssc;
my $cSSC=0;
print STDERR "\nReading ssc input file: $ssc\n";
open(C,"<$ssc" || die "can't open file $ssc: $!");
while(my $l = <C>){
    chomp($l);
    my @tmp = split(/\t/,$l);
    if($tmp[0] =~ /^[^\_]+\_(\d+)\_[+|-]$/){
	#key=coordinate; value=counts;
	$ssc{$1}=$tmp[1];
	$cSSC=$cSSC+1;
    }
}
close(C);
print STDERR "$cSSC input lines from $ssc\n"; 


# astalavista GTF
print STDERR "\nReading asta input file: $asta\n";

open(O,">$out"); #open output file
open(A,"<$asta" || die "can't open file $asta: $!");
while(my $line = <A>){
    chomp($line);
    my @tmp = split(/\t|\"/,$line);
    
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
    
    if ($structure =~ /^0\,1\-2\^$/) {# exon skipping SINGLE 
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
	    if(defined $ssj{$f1."_".$sc2[0]}){
		$a = $ssj{$f1."_".$sc2[0]};
	    }
	    if(defined $ssj{$sc2[1]."_".$f2}){
		$b = $ssj{$sc2[1]."_".$f2};
	    }
	    if(defined $ssj{$flank}){
		$c = $ssj{$flank}; 
	    }
	    
	    my $psi;
	    
	    if( ($a+$b+$c) > $threshold ){
		$psi = (($a + $b)/($a + $b + 2*$c));	
	    }
	    else{
		$psi="NA";
	    }
	    #print output
	    print O $line."psi \"$psi\";\n";
	    if($verbose){
		print "exon_skipping_single\t$structure\t$psi\n";
	    }
	}else{
	    print STDERR "Splice chain $splice_chain does not have the expected format $structure\n";
	}
    }
    
    elsif($structure =~/^0\,1\^2\-$/){# intron retention
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
	    if(defined $ssj{$sc2[0]."_".$sc2[1]}){
		$a = $ssj{$sc2[0]."_".$sc2[1]};
	    }
	    if(defined $ssc{$sc2[0]}){
		$b = $ssc{$sc2[0]};
	    }
	    if(defined $ssc{$sc2[1]}){
		$c = $ssc{$sc2[1]};
	    }
	    
	    my $psi;
	    
	    if(($a+$b+$c) > $threshold ){
		$psi = ((2*$a)/($b + $c + 2*$a));		
	    }else{
		$psi="NA";
	    }
	    #print output
	    print O $line."psi \"$psi\";\n";
	    
	    if($verbose){
		print "intron_retention\t$structure\t$psi\n";
	    }	
	}else{
	    print STDERR "Splice chain $splice_chain does not have the expected format $structure\n";
	}
    } 
    
    elsif($structure=~/^1\-2\^\,3\-4\^$/){# mutually exclusive
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
		if(defined $ssj{$f1."_".$sc1[0]}){
		    $c = $ssj{$f1."_".$sc1[0]};
		}
		if(defined $ssj{$sc1[1]."_".$f2}){
		    $d = $ssj{$sc1[1]."_".$f2};
		}
		if(defined $ssj{$f1."_".$sc2[0]}){
		    $a = $ssj{$f1."_".$sc2[0]};
		}
		if(defined $ssj{$sc2[1]."_".$f2}){
		    $b = $ssj{$sc2[1]."_".$f2};
		}
		
		my $psi;
		
		if(($a+$b+$c+$d) > $threshold ){
		    $psi = (($c + $d)/($c + $d + $a+ $b));		
		}else{
		    $psi="NA";
		}
		#print output
		print O $line."psi \"$psi\";\n";
		
		if($verbose){
		    print "mutually_exclusive\t$structure\t$psi\n";
		}
	    }
	}else{
	    print STDERR "Splice chain $splice_chain does not have the expected format $structure\n";
	}
    }
    elsif($structure=~/^1\^\,2\^$/){# alternative donor
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
	    if(defined $ssj{$sc1."_".$f2}){
		$a = $ssj{$sc1."_".$f2};
	    }
	    if(defined $ssj{$sc2."_".$f2}){
		$b = $ssj{$sc2."_".$f2};
	    }
	    
	    my $psi;
	    
	    if( ($a+$b) > $threshold ){
		$psi = (($a)/($a+ $b));		
	    }else{
		$psi="NA";
	    }
	    #print output
	    print O $line."psi \"$psi\";\n";
	    
	    if($verbose){
		print "alternative_donor\t$structure\t$psi\n";
	    }	    
	}
       
    }
    elsif($structure=~/^1\-\,2\-$/){# alternative acceptor
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
	    if(defined $ssj{$f1."_".$sc1}){
		$a = $ssj{$f1."_".$sc1};
	    }
	    if(defined $ssj{$f1."_".$sc2}){
		$b = $ssj{$f1."_".$sc2};
	    }
	    
	    my $psi;
	    
	    if( ($a+$b) > $threshold ){
		$psi = (($a)/($a + $b));		
	    }else{
		$psi="NA";
	    }
	    #print output
	    print O $line."psi \"$psi\";\n";
	    
	    if($verbose){
		print "alternative_acceptor\t$structure\t$psi\n";
	    }
	}
	
    }
    elsif($structure=~/^0\,1\-2\^3\-4\^/){# exon skipping MULTI
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
		if(defined $ssj{$f1."_".$sc2[0]}){
		    $numeratorPSI = $ssj{$f1."_".$sc2[0]}; #first junction
		}
		
		for(my $i=1; $i<($n-2); $i++){#interior juntions, numerator and denominator get the same values
		    if(defined  $ssj{$sc2[$i]."_".$sc2[$i+1]}){
			$numeratorPSI = $numeratorPSI + $ssj{$sc2[$i]."_".$sc2[$i+1]};
		    }
		}
		#last junction(s_n,f2)
		if(defined $ssj{$sc2[$n-1]."_".$f2}){
		    $numeratorPSI = $numeratorPSI + $ssj{$sc2[$n-1]."_".$f2};
		}
		#last term for the denominator only (n-1)*juntion(f1,f2)
		$denominatorPSI = $numeratorPSI;
		$sum = $numeratorPSI;
		if(defined $ssj{$f1."_".$f2}){
		    $denominatorPSI = $denominatorPSI + (($n-1)*($ssj{$f1."_".$f2}));
		    $sum = $sum + $ssj{$f1."_".$f2};
		}
		if($sum > $threshold){
		    $psi = $numeratorPSI / $denominatorPSI; 
		}else{
		    $psi="NA";
		}
		
		#print output
		print O $line."psi \"$psi\";\n";
		
		if($verbose){
		    print "exon_skipping_multi\t$structure\t$psi\n";
		}
	    }
	}
    } # close the last case of structure, otherwise complex AS/DSP/VST event and PSI=NA
}
close(A);
close(O);
