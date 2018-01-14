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
### NAME : orderFamily_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Orders the IDF column of the RGSL by chromosome and by position within chromosomes. 
###   
### OUTPUT : The RGSL with the ordered IDF column.
###     
###     perl orderFamily_v1.0.pl -inputFile /url/bowtie/RGSL##_##_strings_kbp.landscape_sort -rgslId RGSL# -kmerFile /url/Kmers/RG#_#_strings_kbp.fna -outputFile /url/
###      
### OPTIONS :
###   -inputFile : The RGSL file.
###   -outputFile : Path of output file.
###   -rgslId : Reference Genome Self Landscape unique identifier (RGSL#, the # indicates the numeric component of the RGSL unique identifier).
###   -kmerFile : Length of kmer.
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
GetOptions (\%opts, 
	'inputFile=s', 
	'rgslId=s', 
	'kmerFile=s', 
	'outputFile=s');

&readArguments();
&main();

sub main {

	my $inputFile =  $opts{inputFile};
	my $outputFile = $opts{outputFile};
	my $kmerFile = $opts{kmerFile};
	my $rgslId = $opts{rgslId};

	my @structureInputFile;
	my @structureKmer; 
	my %idFoundKmer = ();

    print "\n\nLoading:  ${inputFile}\n";

	open(ININPUTFILE, $inputFile) or die("Can't read ${inputFile} file");
 		@structureInputFile = <ININPUTFILE>;
	close(ININPUTFILE);

    print "\n\nLoading:  ${kmerFile}\n";

	open(KMER_FILE, $kmerFile) or die("Can't read ${kmerFile} file");
 		@structureKmer = <KMER_FILE>;
	close(KMER_FILE);

    foreach my $indexKmer (0 .. $#structureKmer) {
        my $currentLine = $structureKmer[$indexKmer];
        chop $currentLine;

        if ($currentLine =~ /^>/) {
			my $nextLine = $structureKmer[$indexKmer+1];
        	chop $nextLine;
           	$idFoundKmer{$currentLine} = $nextLine; 
        }     

    }

    open(OUTPUT_FILE, ">".${outputFile}) or die("Can't create file");	

	foreach my $indexInputFile (0 .. $#structureInputFile) {
        	
		my $currentLine = $structureInputFile[$indexInputFile];
		chop $currentLine;

    	my @currentLineArray = split("\t", $currentLine);
		my @ids = split(",", $currentLineArray[4]);
		my %idFound;

		foreach my $posIDs (0 .. $#ids) {
				
				my $id = $ids[$posIDs];
				my @gr = split("_", $id);
				
				$idFound{$gr[1]}{$gr[2]}   = $gr[2];
		}
		my $seq = "";
		if ($idFoundKmer{">".$currentLineArray[1]}){
			$seq= $idFoundKmer{">".$currentLineArray[1]};
		}
		my $numIDsOrdered = 0;
		my $currentId = $currentLineArray[1];
		my @currentIdSplit = split("_", $currentId);

		print OUTPUT_FILE $rgslId."_".$currentIdSplit[1]."_".$currentIdSplit[2]."\t".$currentLineArray[2]."\t".$seq."\t";

		foreach my $chromosome (sort keys %idFound) {
			for $position ( sort {$a<=>$b} keys %{$idFound{$chromosome}} ) {
        		print OUTPUT_FILE  $rgslId."_".$chromosome."_".$idFound{$chromosome}{$position};
        		if ( $numIDsOrdered  < $#ids ) { 
      				print OUTPUT_FILE  ",";
   				}
   				$numIDsOrdered++;
    		}
		}
		print OUTPUT_FILE  "\n";
	}	
	close(OUTPUT_FILE);
}

### Read arguments
sub readArguments {
	my $mandatoryParameters = 'true';

	if (!$opts{inputFile}){
		$opts{inputFile} = '';
	}
	if (!$opts{outputFile}){
		$opts{outputFile} = '';
	}
	if (!$opts{kmerFile}){
		$opts{kmerFile} = '';
	}
	if (!$opts{rgslId}){
        $opts{rgslId} = '';
    }

	### Mandatory parameters
	if ($opts{inputFile} eq ''){
		print ("Needs the -inputFile parameter \n");
		$mandatoryParameters = 'false';
	}
	if ($opts{outputFile} eq ''){
		print ("Needs the -outputFile parameter \n");
		$mandatoryParameters = 'false';
	}	
	if ($opts{kmerFile} eq ''){
		print ("Needs the -kmerFile parameter \n");
		$mandatoryParameters = 'false';
	}	
	if ($opts{rgslId} eq ''){
    	print ("Needs the -rgslId parameter \n");
    	$mandatoryParameters = 'false';
    }

    if ($mandatoryParameters eq 'false'){
  		exit(0);
	}
}
