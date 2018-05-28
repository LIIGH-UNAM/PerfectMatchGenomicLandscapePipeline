#!/usr/bin/perl -w

use strict;
use LWP::Simple;
use File::Copy;
use Getopt::Long;
use Math::Round;
use Math::BigFloat;
use warnings;

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
### NAME : signatureValue_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Generates the SV column of the PMGL to be generated (each Reference String's signature value). 
###   
### OUTPUT : File containing the RGSL, PM, PMnCR, and SV columns of the PMGL to be generated. 
###     
###
### USAGE : perl signatureValue_v1.0.pl -inputFile PMGL#_RG#_all_PMnCR.tab
###      
### OPTIONS :
###   -inputFile : File containing the RGSL, PM, and PMnCR columnns of the PMGL to be generated. 
###   -h or -help
###   
### DATE : 01/10/2017
###
### Requirements : 
### 1)
###   Perl's  Libraries
###     Getopt::Long;
###
###
End_of_Help
    close HELP;
    exit(0);
}

my %opts;

### Parameters 
GetOptions (\%opts,'inputFile=s');

&readArguments();
&main();

sub main {  

	my $inputFile =  $opts{inputFile};

    print "\n\nLoading:  ${inputFile}\n";

	my @structureInputFile;

	open(INPUTFILE, $inputFile) or die("Can't read ${inputFile} file");
    	@structureInputFile = <INPUTFILE>;
    close(INPUTFILE);

	$inputFile =~ s/.tab//;

	open(OUTPUT_FILE, ">".$inputFile."_SV.tab") or die("Can't create ${inputFile}_SV.tab file");

		print OUTPUT_FILE "ID\tCR\tSEQ\tPM\tPMnCR\tSV\n";

		my $currentLine = $structureInputFile[0];
		chop $currentLine;

		my @currentLineArray = split("\t", $currentLine);

		print OUTPUT_FILE $currentLineArray[0]."\t".$currentLineArray[1]."\t".$currentLineArray[2]."\t".$currentLineArray[3]."\t".$currentLineArray[4]."\t"."NA"."\n";

        foreach my $structureInputFileIndex (1 .. $#structureInputFile) {

        	$currentLine = $structureInputFile[$structureInputFileIndex];
			my $previousLine = $structureInputFile[$structureInputFileIndex-1];
			
			chop $currentLine;
			chop $previousLine;

			@currentLineArray = split("\t", $currentLine);
			my @previousLineArray = split("\t", $previousLine);

			if (defined $currentLineArray[0] && defined $currentLineArray[1] && defined $currentLineArray[2] && defined $currentLineArray[3] && defined $currentLineArray[4]){
				if (defined $previousLineArray[4]){
					my $currentID = $currentLineArray[0];
					my $previousID = $previousLineArray[0];

					#$currentID =~ m/.*_(\d\d)_(\d+)/;
					$currentID =~ m/.*_(\d+)_(\d+)/;
					my $currentChromosome = $1;
					my $currentPosition = $2;

					#$previousID =~ m/.*_(\d\d)_(\d+)/;
					$previousID =~ m/.*_(\d+)_(\d+)/;
					my $previousChromosome = $1;
					my $previousPosition = $2;
					
					if ($previousChromosome == $currentChromosome && ($previousPosition+1) == $currentPosition){
						my $sv = ($currentLineArray[4]+1)/($previousLineArray[4]+1);
		               	$sv = sprintf("%.2f",$sv);
			           
						print OUTPUT_FILE $currentLineArray[0]."\t".$currentLineArray[1]."\t".$currentLineArray[2]."\t".$currentLineArray[3]."\t".$currentLineArray[4]."\t".$sv."\n";
					}else{
						print OUTPUT_FILE $currentLineArray[0]."\t".$currentLineArray[1]."\t".$currentLineArray[2]."\t".$currentLineArray[3]."\t".$currentLineArray[4]."\t"."NA"."\n";	
					}
				}else{
					print OUTPUT_FILE $currentLineArray[0]."\t".$currentLineArray[1]."\t".$currentLineArray[2]."\t".$currentLineArray[3]."\t".$currentLineArray[4]."\t"."NA"."\n";	
				}

			}	
        }
        close(OUTPUT_FILE);
}

### Read arguments
sub readArguments {
	my $mandatoryParameters = 'true';

	if (!$opts{inputFile}){
		$opts{inputFile} = '';
	}

	### Mandatory parameters
	if ($opts{inputFile} eq ''){
		print ("Needs the -inputfile parameter \n");
		$mandatoryParameters = 'false';
	}

    if ($mandatoryParameters eq 'false'){
		exit(0);
	}
}
