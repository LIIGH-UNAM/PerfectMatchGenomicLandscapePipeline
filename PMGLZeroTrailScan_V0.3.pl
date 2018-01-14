#!/usr/bin/perl -w

use strict;
use Getopt::Long;

#print "arg\n";
my $arg = $ARGV[0];
if (!$arg){
 $arg = '';
}

if ($arg eq "-h" || $arg eq "-help") {
    open HELP, "| more";
    print HELP <<End_of_Help;
###
### AUTHOR : 
###     Jair Santiago Garcia Sotelo, send comments to jsgarcia\@liigh.unam.mx
###     Kim Palacios Flores, send comments to kimpalaciosflores\@gmail.com
###
### NAME : PMGLZeroTrailScan.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Locates signatures of variation along a Perfect Match Genomic Landscape (PMGL) using the zero-trail scan.
###	
### OUTPUT : File containing the PMGL rows corresponding to the Downstream Recovery String of each signature of variation.
###
### USAGE : perl PMGLZeroTrailScan.pl -inputFile PMGL#_RG#_all_PMnCR_SV.tab -minPMn # -maxPMn_1 # -CRn_1 # -lowComplexity -zeroTrail
###      
### OPTIONS :
###	-inputfile : The complete PMGL file.
###	-minPMn : minimum number of normalized perfect matches at position n.
###	-maxPMn_1 : maximum number of normalized perfect matches at position n-1.
###	-CRn_1 : Count Reference at position n-1.
###	-lowComplexity (optional) : Filter to eliminate signatures of variation with a low complexity Downstream Recovery String.
### -zeroTrail (optional) : Filter to eliminate signatures of variation with less than 16 near-zero values in the zero-trail.
### -outputDir : Path of output directory.
###
###	-h or -help
###	
### DATE : 01/05/2016
###
### Requirements : 
### 1)
### 	Perl's  Libraries
###     Getopt::Long;
###
###
End_of_Help
    close HELP;
    exit(0);
}

my %opts;

### Parameters 
GetOptions (\%opts,
	'inputFile=s', 
	'minPMn=i', 
	'maxPMn_1=i', 
	'CRn_1=i', 
	'lowComplexity', 
	'zeroTrail', 
	'outputDir=s');

&ReadArguments();
&main();

	my @structureInputFile;

sub main {

	my $inputFile =  $opts{inputFile};
	my $minPMn = $opts{minPMn};
	my $maxPMn_1 = $opts{maxPMn_1};
	my $CRn_1 = $opts{CRn_1};
	my $lowComplexity = $opts{lowComplexity};
	my $zeroTrail = $opts{zeroTrail};
	my $outputDir = $opts{outputDir};


    print "\nLoading:  ${inputFile}\n";

	open(IN, $inputFile) or die("Can't read ${inputFile} file");
 		@structureInputFile = <IN>;
	close(IN);

	$inputFile =~ /(.*\/)(.*)\.tab$/;
  	my $path = $1;
  	my $nameFile = $2;

	open(OUTPUT_FILE_SIGNATURES, ">".$outputDir.$nameFile."_ZeroTrailScan.tab") or die("Can't create ".$inputFile."_ZeroTrailScan.tab file");
	
		foreach my $position (2 .. $#structureInputFile) {
			my $currentLine = $structureInputFile[$position];
			my $previousLine = $structureInputFile[$position-1];
			my $printID= "true";

			chop $currentLine;
			chop $previousLine;


			if ($currentLine !~ /^#/) {
	        	my @currentLineArray = split("\t", $currentLine);
	        	my @previousLineArray = split("\t", $previousLine);
	        	my $totalZero=0;

	        	if (defined $currentLineArray[0] && defined $currentLineArray[1] && defined $currentLineArray[2] && defined $currentLineArray[3] && defined $currentLineArray[4] && defined $currentLineArray[5] && defined $previousLineArray[5]){
					if ( $currentLineArray[4] >= $minPMn && $previousLineArray[4] <= $maxPMn_1 && $previousLineArray[1] == $CRn_1 && $currentLineArray[5] ne "NA"){
						if ($lowComplexity){
							$printID = &LowComplexity($currentLineArray[2]);
						}
						if ($zeroTrail  && $printID eq "true"){
							$totalZero = &ZeroTrailScan($position);
							if ($totalZero >= 16){
								print OUTPUT_FILE_SIGNATURES $currentLineArray[0]."\t".$currentLineArray[1]."\t".$currentLineArray[2]."\t".$currentLineArray[3]."\t".$currentLineArray[4]."\t".$currentLineArray[5]."\t".$totalZero."\n";
							}
		        		}elsif($printID eq "true"){
							print OUTPUT_FILE_SIGNATURES $currentLineArray[0]."\t".$currentLineArray[1]."\t".$currentLineArray[2]."\t".$currentLineArray[3]."\t".$currentLineArray[4]."\t".$currentLineArray[5]."\n";
		        		}
					}
				}
			}
    	}
	close(OUTPUT_FILE_SIGNATURES);
}


### Zero-trail filter
sub ZeroTrailScan {
	my $position= shift;
	my $totalZero=0;

	for (my $positionFile = ($position-1); $positionFile >= 0; $positionFile--) {
		my $currentLine = $structureInputFile[$positionFile];
		chop $currentLine;
		my @currentLineArray = split("\t", $currentLine);
		if ($currentLineArray[4] <= 2){
			$totalZero++;
		}
       last if ($currentLineArray[4] > 2);
    }
	return $totalZero;
}

### Read arguments
sub ReadArguments {

	my $mandatoryParameters = 'true';

	if (!$opts{inputFile}){
		$opts{inputFile} = '';
	}
	if (!$opts{minPMn}){
		$opts{minPMn} = '';
	}
	if (!$opts{maxPMn_1}){
		$opts{maxPMn_1} = '';
	}
	if (!$opts{CRn_1}){
		$opts{CRn_1} = '';
	}
	if (!$opts{outputDir}){
		$opts{outputDir} = '';
	}

	### Mandatory parameters
	if ($opts{inputFile} eq ''){
		print ("Needs the -inputFile parameter \n");
		$mandatoryParameters = 'false';
	}elsif (! -f $opts{inputFile}){
		print ("\t-inputFile File does not exist: $opts{inputFile}\n");
		$mandatoryParameters = 'false';
	}elsif ($opts{inputFile} !~ /PMGL\d+_RG\d+_all_PMnCR_SV\.tab$/){
		print ("\t-inputFile Required the following nomenclature: PMGL#_RG#_all_PMnCR_SV.tab \n");
		$mandatoryParameters = 'false';
	}

	if ($opts{minPMn} eq ''){
		print ("Needs the -minPMn parameter \n");
		$mandatoryParameters = 'false';
	}

	if ($opts{maxPMn_1} eq ''){
		print ("Needs the -maxPMn_1 parameter \n");
		$mandatoryParameters = 'false';
	}

	if ($opts{CRn_1} eq ''){
		print ("Needs the -CRn_1 parameter \n");
		$mandatoryParameters = 'false';
	}

	if ($opts{outputDir} eq ''){
		print ("Needs the -outputDir parameter \n");
		$mandatoryParameters = 'false';
	}elsif (! -d $opts{outputDir}){
		print ("\t-outputDir Directory does not exist: $opts{outputDir}\n");
		$mandatoryParameters = 'false';
	}

    if ($mandatoryParameters eq 'false'){
  		exit(0);
	}
}

### Low complexity filter
sub LowComplexity {
    my $sequence = shift;
    my $lowComplexity= "false";
	my $length = length($sequence);
	my $a=0;
	my $t=0;
	my $g=0;
	my $c=0;

	for (my $countLetters = 0; $countLetters < $length; $countLetters++){
		my $character = substr($sequence,$countLetters,1);
		if (uc $character eq "A"){
			$a++;
		}elsif (uc $character eq "T"){
			$t++;
		}elsif (uc $character eq "G"){
			$g++;
		}elsif (uc $character eq "C"){
			$c++;
		}
	}

	my $rule = 0;

	if ($a >= 1){
		$rule ++;
	}
	if ($t >= 1){
		$rule ++;
	}
	if ($g >= 1){
		$rule ++;
	}
	if ($c >= 1){
		$rule ++;
	}
	if ($rule > 2 && $a < 15 && $t < 15  && $g < 15  && $c < 15 ){
			$lowComplexity = "true";
	}

    return $lowComplexity;
}

