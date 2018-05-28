#!/usr/bin/perl -w

#use strict;
use LWP::Simple;
use File::Copy;
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
### NAME : normalizedByCountReference_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Generates the PMnCR column of the PMGL to be generated (each Reference String's PM normalized by its CR). 
###   
### OUTPUT : File containing the RGSL, PM, and PMnCR columns of the PMGL to be generated. 
###     
###
### USAGE : perl normalizedByCountReference_v1.0.pl -inputFile PMGL#_RG#_all.tab
###      
### OPTIONS :
###   -inputFile : File containing the RGSL and PM columnn of the PMGL to be generated. 
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
	my @structureInputFile;

    print "\n\nLoading:  ${inputFile}\n";

	open(INPUT_FILE, $inputFile) or die("Can't read ${inputFile} file");
 		@structureInputFile = <INPUT_FILE>;
	close(INPUT_FILE);

	$inputFile =~ s/.tab//;
	
	open(OUTPUT_FILE, ">".$inputFile."_PMnCR.tab") or die("Can't create ${inputFile}_PMnCR.tab file");

	foreach my $structureInputFileIndex (0 .. $#structureInputFile) {    	
		my $currentLine = $structureInputFile[$structureInputFileIndex];
		chop $currentLine;
			if ($currentLine !~ /^#/) {
        		my @currentLineArray = split("\t", $currentLine);
				if (defined $currentLineArray[0] && defined $currentLineArray[1]&& defined $currentLineArray[2 ] && defined $currentLineArray[3]){
					my $pMnCR = ($currentLineArray[3])/$currentLineArray[1];
					$pMnCR = sprintf("%.2f", $pMnCR);
					print OUTPUT_FILE $currentLineArray[0]."\t".$currentLineArray[1]."\t".$currentLineArray[2]."\t".$currentLineArray[3]."\t".$pMnCR."\n";
				}else{
					print OUTPUT_FILE $currentLine."\n";
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