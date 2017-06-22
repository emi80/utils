#! /usr/bin/perl
use strict; use warnings; 
use Getopt::Long;

#CCK & AB <15/may/2015>

# default values
my $asta;	# option variable with default value
my $ssj;	# option variable with default value
my $ssc;	# option variable with default value
my $out="astalavista.psi.out.gtf";
my $verbose;
my $help;

GetOptions(
    'asta|a=s' => \$asta,
    'ssj|b=s' => \$ssj,
    'ssc|c=s' => \$ssc,
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
    --asta|a    : pairwise AS events from astalavista output gtf file (1-based coordinates)
    --ssj|b     : splice site junction quantifications from ipsa pipeline BED file (1-based offset)
    --ssc|c     : splice site boundary quantifications from ipsa pipeline BED file (1-based offset)
    --out|o     : output file name [optional] (default "astalavista.psi.out.gtf")
    --verbose|v : foreach event print type_of_event structure PSI
    --help
    
  Exemple:
    perl $0 --asta astalavista.gtf --ssj sample.ssj.bed --ssc sample.ssc.bed

  Information on PSI calculation for the different AS events:
    display ~cklein/utils/usage/psi.jpg

USAGE
    exit();
}
###########################################################################################


# ssj bed
my %ssj;
open(J,"<$ssj" || die "can't open file $ssj: $!");
while(my $l = <J>){
    chomp($l);
    my @tmp = split(/\t/,$l);
    my $k = "$tmp[1]_$tmp[2]"; #key=coordinate_coordinate; value=counts;
    $ssj{$k}=$tmp[6];
}
close(J);


# ssc bed
my %ssc;
open(C,"<$ssc" || die "can't open file $ssc: $!");
while(my $l = <C>){
    chomp($l);
    my @tmp = split(/\t/,$l);
    #key=coordinate; value=counts;
    $ssc{$tmp[1]}=$tmp[6];
}
close(C);


# astalavista GTF
open(O,">$out"); #open output file
open(A,"<$asta" || die "can't open file $asta: $!");
while(my $line = <A>){
    chomp($line);
    my @tmp = split(/\t|\"/,$line);
    my $flank = "$tmp[3]_$tmp[4]"; 
    my $structure = $tmp[15]; 
    my $splice_chain = $tmp[17]; 

    if($tmp[2] eq "as_event"){
	my ($psi,$flag) = computePSI(\%ssj,\%ssc,$flank,$structure,$splice_chain);
	if($flag){#print only the simple AS events for which we are computing PSI (see usage)
	    print O $line."psi \"$psi\";\n";
	}
    }

}
close(A);
close(O);

###################

sub computePSI {
    my %j = %{shift()};
    my %c = %{shift()};
    my $f = shift();
    my $s = shift();
    my $sc=shift();

    my ($f1,$f2) = split(/\_/,$f); #flanking 1 and 2
    my ($sc1,$sc2) = split(/\,/,$sc); #splice_chain from 1st/2nd transcripts
    my @sc1 = split(/\^|\-|\[|\]/,$sc1); #s1 is $sc1[0] and s2 is $sc1[1] (from the 1st transcripts)
    my @sc2 = split(/\^|\-|\[|\]/,$sc2); #s1 is $sc2[0] and s2 is $sc2[1] (from the 2nd transcripts)

    #return variables from this subroutine
    my $p; #psi value
    my $flagSimpleEvent=0; #flag to output only simple events
    

    
    for ($s) {
	
	if (/^0\,1\-2\^$/) {# exon skipping SINGLE 
	    # psi=(a+b)/(a+b+2c)
	    # only need ssj
	    # complete:psi=1
	    # skipped:psi=0
	    # a=f1,s1;b=s2,f2; (sorted: always from sc2)
            # c=f1,f2;
	    
	    $flagSimpleEvent=1;

	    my $a = $j{$f1."_".$sc2[0]};
	    my $b = $j{$sc2[1]."_".$f2};
	    my $c = $j{$f}; 
	    
	    if(defined $a && defined $b && defined $c){
		$p = (($a + $b)/($a + $b + 2*$c));	
	    }
	    else{
		$p="NA";
	    }
	    
	    if($verbose){
		print "exon_skipping_single\t$s\t$p\n";
	    }
	}

	elsif(/^0\,1\^2\-$/){# intron retention
	    # psi=(2a)/(b+c+2a)
	    # no intron:psi=1
	    # retained:psi=0
	    # requires ssj (a) and ssc (b and c)
	    # a=s1,s2; (sorted: always from sc2)
            # b=s1; c=s2;

	    $flagSimpleEvent=1;

	    my $a = $j{$sc2[0]."_".$sc2[1]};
	    my $b = $c{$sc2[0]};
	    my $c = $c{$sc2[1]};
	    
	    if(defined $a && defined $b && defined $c){
		$p = ((2*$a)/($b + $c + 2*$a));		
	    }else{
		$p="NA";
	    }

	    if($verbose){
		print "intron_retention\t$s\t$p\n";
	    }

	} 
	    
	elsif(/^1\-2\^\,3\-4\^$/){# mutually exclusive
	    # psi=(c+d)/(c+d+a+b)
	    # requires only ssj
	    # c=f1,s1;d=s2,f2; (sorted: always from sc1)
            # a=f1,s3; b=s4,f2;(sorted: always from sc2)
	    
	    $flagSimpleEvent=1;

	    my $c = $j{$f1."_".$sc1[0]};
	    my $d = $j{$sc1[1]."_".$f2};
	    my $a = $j{$f1."_".$sc2[0]};
	    my $b = $j{$sc2[1]."_".$f2};

	    if(defined $a && defined $b && defined $c && defined $d){
		$p = (($c + $d)/($c + $d + $a+ $b));		
	    }else{
		$p="NA";
	    }
	    
	    if($verbose){
		print "mutually_exclusive\t$s\t$p\n";
	    }

	}
	elsif(/^1\^\,2\^$/){# alternative donor
	    # psi=(a)/(a+b)
	    # requires only ssj
	    # a=s1,f2; (sorted: always from sc1)
            # b=s2,f2; (sorted: always from sc2)

	    $flagSimpleEvent=1;

	    my $a = $j{$sc1[0]."_".$f2};
	    my $b = $j{$sc2[0]."_".$f2};
	    
	    if(defined $a && defined $b){
		$p = (($a)/($a+ $b));		
	    }else{
		$p="NA";
	    }
	    
	    if($verbose){
		print "alternative_donor\t$s\t$p\n";
	    }	    

	}
    	elsif(/^1\-\,2\-$/){# alternative acceptor
	    # psi=(a)/(a+b)
	    # requires only ssj
	    # a=f1,s1; (sorted: always from sc1)
            # b=f1,s2; (sorted: always from sc2)
	    
	    $flagSimpleEvent=1;

	    my $a = $j{$f1."_".$sc1[0]};
	    my $b = $j{$f1."_".$sc2[0]};
	    
	    if(defined $a && defined $b){
		$p = (($a)/($a + $b));		
	    }else{
		$p="NA";
	    }

	    if($verbose){
		print "alternative_acceptor\t$s\t$p\n";
	    }

	}
	elsif(/^0\,1\-2\^3\-4\^/){# exon skipping MULTI
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
	    my @matches = ($s =~ /(\d+\-\d+\^)/g);
	    my $l = @matches; #number of skipped exons
	    my $n = 2*$l; # number of splice sites
	    my $string = $s;
	    $string =~ s/^0\,//; #remove "0," from the beginning of the sequence
	    $string =~ s/(\d+\-\d+\^)//g; #remove as many times as the pattern appears
	   
	    if($string eq ""){ #if nothing is left, it is only MULTI exon skipping

		$flagSimpleEvent=1;

		if(defined $j{$f1."_".$sc2[0]}){ #if first junction is defined in the ssj file

		    my $numeratorPSI = $j{$f1."_".$sc2[0]}; #first junction
		    my $denominatorPSI = $j{$f1."_".$sc2[0]}; #first junction
		    my $flag=1; #flag for defined internal junctions

		    for(my $i=1; $i<($n-2); $i++){#interior juntions, numerator and denominator get the same values
			if(defined  $j{$sc2[$i]."_".$sc2[$i+1]}){
			    $numeratorPSI = $numeratorPSI + $j{$sc2[$i]."_".$sc2[$i+1]};
			    $denominatorPSI = $denominatorPSI + $j{$sc2[$i]."_".$sc2[$i+1]};
			}else{
			    $flag=0; # any internal juntion not defined
			} 
		    }
		    
		    #last junction and check if internal junctions are defined 
		    if(defined $j{$sc2[$n-1]."_".$f2} && defined $j{$f1."_".$f2} && $flag eq 1){

			#for numerator and denominator junction(s_n,f2)
			$numeratorPSI = $numeratorPSI + $j{$sc2[$n-1]."_".$f2};
			$denominatorPSI = $denominatorPSI + $j{$sc2[$n-1]."_".$f2};

			#for denominator only (n-1)*juntion(f1,f2)
			$denominatorPSI = $denominatorPSI + (($n-1)*($j{$f1."_".$f2})); 
			$p = $numeratorPSI / $denominatorPSI; 
		    }else{
			$p="NA";
		    }
		
		}
		else{#first juntion not defined
		    $p="NA";
		}

		if($verbose){
		    print "exon_skipping_multi\t$s\t$p\n";
		}
		
	    }
	    else{#a more complex pattern than multi exon skipping
		$p = "NA";
	    }
	    
	} # close the last case of structure, otherwise complex AS/DSP/VST event and PSI=NA
	
	else { 
	    $p="NA";
	}

    }
  
    return ($p,$flagSimpleEvent);
}
