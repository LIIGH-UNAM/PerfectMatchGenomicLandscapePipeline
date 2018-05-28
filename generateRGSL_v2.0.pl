#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Copy;
use FindBin qw( $RealBin );
require "${RealBin}/library.pm";


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
### NAME : generateRGSL_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Generates a Reference Genome Self Landscape (RGSL) from a Reference Genome sequence in fasta format. 
###   
### OUTPUT : RGSL#.tab file. The RGSL reports each Reference String's unique identifier (ID), the number of times its sequence is present 
###	in the entire Reference Genome (CR), its DNA sequence (SEQ), and the unique identifiers of all Reference Strings in the entire Reference 
### Genome that share the same sequence (IDF). 
###     
### USAGE : perl generateRGSL_v1.0.pl -binDir /url/ -fastaDir /url/RG#_#.fasta -bowtieDir /url/ -outputDir /url/ -kmerLength # -rgslId RGSL# -memory #
###      
### OPTIONS :
###   -binDir : Path of directory containing all PMGL pipeline scripts.
###   -fastaDir : Path of directory containing the Reference Genome sequence per chromosome in fasta format. 
###			The files contained in fastaDir require the following nomenclature: RG#_#.fasta (The first # indicates the numeric component 
###			of the RG unqiue identifier, the second # indicates the chromosome number).
###   -bowtieDir : Path of directory containing bin for bowtie execution.
###   -outputDir  : Path of output directory.
###   -kmerLength : Length of kmer.
###   -rgslId : Reference Genome Self Landscape unique identifier (RGSL#, the # indicates the numeric component of the RGSL unique identifier).
###   -memory : RAM memory to be used for bowtie execution.
###   -h or -help
###   
### DATE : 01/10/2017
###
### Requirements : 
### 1)
###		Perl's  Libraries
###     	Getopt::Long;
### 2)
###  	Bowtie (0.12.7)
###		
###	
###
End_of_Help
    close HELP;
    exit(0);
}

my %opts;

### Parameters 
GetOptions (\%opts,
	'binDir=s', 
	'fastaDir=s', 
	'bowtieDir=s', 
	'outputDir=s', 
	'kmerLength=i', 
	'rgslId=s', 
	'memory=i'
	);

&readArguments();
&main();

sub main {

my $binDir = $opts{binDir};
my $fastaDir = $opts{fastaDir};
my $bowtieDir = $opts{bowtieDir};
my $outputDir = $opts{outputDir};
my $kmerLength = $opts{kmerLength};
my $rgslId = $opts{rgslId};
my $memory = $opts{memory};

my @fastaFileArray;

#### Generating Directory

	my $outFasta = $outputDir."fasta/";
	my $outBowtie = $outputDir."bowtie/";
	my $outKMers = $outputDir."Kmers/";
	my $outLN = $outputDir."LN/";
	my $outRGSL = $outputDir."RGSL/";
	
	system ("rm -rf ${outFasta}");
	system ("rm -rf ${outBowtie}");
	system ("rm -rf ${outKMers}");
	system ("rm -rf ${outLN}");
	system ("rm -rf ${outRGSL}");

	mkdir $outFasta;
	mkdir $outBowtie;
	mkdir $outKMers;
	mkdir $outLN;
	mkdir $outRGSL;

#### Concatenating in a single fasta file
print "\nConcatenating in a single fasta file\n";

	opendir(my $fastaFileDirIndex, $fastaDir) || die "Error :( $! \n";  
    while(readdir $fastaFileDirIndex){  
        if (-f $fastaDir . "/" . $_) { 
        	if ($_ =~ /.*fasta/) { 
                     push (@fastaFileArray, $_); 
            }
        } 
    }
    closedir $fastaFileDirIndex;  

    $fastaFileArray[0] =~ /(RG\d+)_.*/;
	my $nameAllFasta = $1."_allChr.fasta";

	foreach my $fastaFileArrayIndex (0 .. $#fastaFileArray) {  
	    print "\tCopying ".$fastaFileArray[$fastaFileArrayIndex]." to ".$nameAllFasta."\n";
		`cat $fastaDir$fastaFileArray[$fastaFileArrayIndex] >> ${outFasta}${nameAllFasta}`
	}

#### Generating kmers
print "\nGenerating kmers\n";

	my $numberOfFastaFiles = @fastaFileArray;

	foreach my $index (0 .. $#fastaFileArray) {
	    print "\tGenerating kmer file ".$fastaFileArray[$index]."\n";

	    $fastaFileArray[$index] =~ /RG\d+_(\d+).fasta/;

	    my $chromosome = $1;

	    $fastaFileArray[$index] =~ /(RG\d+_\d+).fasta/;

	    #cutting sequence into kmers
		system ("perl ${binDir}makeRGKmers_v1.0.pl -inputFile ${fastaDir}$fastaFileArray[$index] -kmerLength $kmerLength -rgslId ${rgslId}_${chromosome} -outputFile ${outKMers}${1}_strings_kbp.fna");
		
		#Generating bash: Report all perfect-match occurrences of each k-mer in the bowtie database
		system ("echo ${bowtieDir}bowtie -v 0 -a ${nameAllFasta} -f ${outKMers}${1}_strings_kbp.fna ${outBowtie}${rgslId}_${chromosome}_strings_kbp.out >> ${outKMers}JOBS_bowtieALL_RGSL.bash");


		#Generating bash: instructions to generate landscapes from bowtie files
		system ("echo perl ${binDir}makeRGSLFile_v1.0.pl -idFamily -inputFile ${outBowtie}${rgslId}_${chromosome}_strings_kbp.out -outputFile ${outBowtie}${rgslId}_${chromosome}_strings_kbp.landscape >> ${outBowtie}JOBS_makeLandscape.bash");

		#Generating bash: ordering all landscape files by position
		system ("echo  sort -nk1 ${outBowtie}${rgslId}_${chromosome}_strings_kbp.landscape \\> ${outBowtie}${rgslId}_${chromosome}_strings_kbp.landscape_sort >> ${outBowtie}JOBS_sortLandscape.bash");

		#Generating bash: ordering IDF
		system ("echo perl ${binDir}orderFamily_v1.0.pl -rgslId ${rgslId} -kmerFile ${outKMers}${1}_strings_kbp.fna -inputFile ${outBowtie}${rgslId}_${chromosome}_strings_kbp.landscape_sort -outputFile ${outLN}${rgslId}_${chromosome}_landscape.tab >> ${outBowtie}JOBS_Order_family.bash");

		if ($?) {
		   exit(0);
		}
	}

####  Generating bowtie database
print "\nGenerating bowtie database\n";

	system ("echo ${bowtieDir}bowtie-build ${outFasta}${nameAllFasta} ${nameAllFasta} > ${outFasta}JOBS_bowtie-build.bash");
	system ("perl ${binDir}makeSGE_v1.0.pl -inputFile ${outFasta}JOBS_bowtie-build.bash -memory ${memory}");

	if ($?) {
	   	exit(0);
	}

	### Executing JOB in cluster
 	print "\tExecuting JOB in cluster \n";
 	
	my $qsubBD = `qsub -t 1-1:1 ${outFasta}JOBS_bowtie-build.sge`;

	$qsubBD =~ /Your job-array (\d+)\.\d+-\d+:\d+ \(\"JOBS_bowtie-build\"\) has been submitted/;
    $qsubBD = $1;

    &waitingForTheCluster($qsubBD);

	`cp ${outFasta}*.ebwt ${bowtieDir}indexes/.`;

### Executing bash: Reporting all perfect-match occurrences of each K-mer in the bowtie database
print "\nReporting all perfect-match occurrences of each K-mer in the bowtie database\n";
	
	system ("perl ${binDir}makeSGE_v1.0.pl -inputFile ${outKMers}JOBS_bowtieALL_RGSL.bash -memory ${memory}");

	print "\tExecuting JOB in cluster \n";
	my $qsubkmerDB = `qsub -t 1-${numberOfFastaFiles}:1 ${outKMers}JOBS_bowtieALL_RGSL.sge`;

	$qsubkmerDB =~ /Your job-array (\d+)\.\d+-\d+:\d+ \(\"JOBS_bowtieALL_RGSL\"\) has been submitted/;
    $qsubkmerDB = $1;

	&waitingForTheCluster($qsubkmerDB);

### Executing bash: instructions to generate landscapes from bowtie files
print "\nExecuting instructions to generate landscapes from bowtie files\n";

	system ("perl ${binDir}makeSGE_v1.0.pl -inputFile ${outBowtie}JOBS_makeLandscape.bash -memory ${memory} ");

	print "\tExecuting JOB in cluster \n";
	my $qsubMakeLandscape = `qsub -t 1-${numberOfFastaFiles}:1 ${outBowtie}JOBS_makeLandscape.sge`;

	$qsubMakeLandscape =~ /Your job-array (\d+)\.\d+-\d+:\d+ \(\"JOBS_makeLandscape\"\) has been submitted/;
    $qsubMakeLandscape = $1;

    &waitingForTheCluster($qsubMakeLandscape);

### Executing bash: ordering all landscape files by position
print "\nOrdering all landscape files by position\n";

	system ("tcsh ${outBowtie}JOBS_sortLandscape.bash");

### Executing bash: ordering IDF
print "\nOrdering IDF\n";

	system ("perl ${binDir}makeSGE_v1.0.pl -inputFile ${outBowtie}JOBS_Order_family.bash -memory ${memory} ");

	print "\tExecuting JOB in cluster \n";
	my $qsubOrderFamily = `qsub -t 1-${numberOfFastaFiles}:1 ${outBowtie}JOBS_Order_family.sge`;

	$qsubOrderFamily =~ /Your job-array (\d+)\.\d+-\d+:\d+ \(\"JOBS_Order_family\"\) has been submitted/;
    $qsubOrderFamily = $1;

    &waitingForTheCluster($qsubOrderFamily);

#Concatenating in a single fasta file
print "\nConcatenating in a single fasta file\n";

	#Header
	system ("echo -e \"ID\\tCR\\tSEQ\\tIDF\" >> ${outRGSL}${rgslId}.tab ");
	my %lnFound = ();
	opendir(my $lnFileDirIndex, $outLN) || die "Error :( $! \n";  
    while(readdir $lnFileDirIndex){  
        if (-f $outLN . "/" . $_) { 
        	if ($_ =~ /.*_landscape.tab/) { 
        		$_ =~ /${rgslId}_(\d+)_landscape\.tab$/;
  				my $chr = $1;
        		$lnFound{$chr} = $_;
            }
        } 
    }
    closedir $lnFileDirIndex;  

    foreach my $cromosoma (sort {$a<=>$b} keys %lnFound) {
    	system ("cat ${outLN}$lnFound{$cromosoma}  >> ${outRGSL}${rgslId}.tab ");
    }  

	print "\nRGSL creado\n\n";
}

### Read arguments
sub readArguments {

	my $mandatoryParameters = 'true';

	if (!$opts{fastaDir}){
		$opts{fastaDir} = '';
	}
	if (!$opts{bowtieDir}){
		$opts{bowtieDir} = '';
	}
	if (!$opts{kmerLength}){
		$opts{kmerLength} = '';
	}
	if (!$opts{outputDir}){
		$opts{outputDir} = '';
	}
	if (!$opts{binDir}){
		$opts{binDir} = '';
	}
	if (!$opts{rgslId}){
		$opts{rgslId} = '';
	}
	if (!$opts{memory}){
		$opts{memory} = '';
	}

	### Mandatory parameters
	if ($opts{fastaDir} eq ''){
		print ("\tNeeds the -fastaDir parameter \n");
		$mandatoryParameters = 'false';
	}elsif (! -d $opts{fastaDir}){
		print ("\t-fastaDir Directory does not exist: $opts{fastaDir}\n");
		$mandatoryParameters = 'false';
	}

	if ($opts{bowtieDir} eq ''){
		print ("\tNeeds the -bowtieDir parameter \n");
		$mandatoryParameters = 'false';
	}elsif (! -d $opts{bowtieDir}){
		print ("\t-bowtieDir Directory does not exist: $opts{bowtieDir}\n");
		$mandatoryParameters = 'false';
	}

	if ($opts{kmerLength} eq ''){
		print ("\tNeeds the -kmerLength parameter \n");
		$mandatoryParameters = 'false';
	}

	if ($opts{outputDir} eq ''){
		print ("\tNeeds the -outputDir parameter \n");
		$mandatoryParameters = 'false';
	}elsif (! -d $opts{outputDir}){
		print ("\t-outputDir Directory does not exist: $opts{outputDir}\n");
		$mandatoryParameters = 'false';
	}

	if ($opts{binDir} eq ''){
		print ("\tNeeds the -binDir parameter \n");
		$mandatoryParameters = 'false';
	}elsif (! -d $opts{binDir}){
		print ("\t-binDir Directory does not exist: $opts{binDir}\n");
		$mandatoryParameters = 'false';
	}

	if ($opts{rgslId} eq ''){
		print ("\tNeeds the -rgslId parameter \n");
		$mandatoryParameters = 'false';
	}elsif($opts{rgslId} !~ /^RGSL\d+$/){
		print ("\t-rgslId Required the following nomenclature: RGSL# \n");
		$mandatoryParameters = 'false';
	}

	if ($opts{memory} eq ''){
		print ("\tNeeds the -memory parameter \n");
		$mandatoryParameters = 'false';
	}

    if ($mandatoryParameters eq 'false'){
  		exit(0);
	}
}


