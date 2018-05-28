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
### NAME : generatePMGL_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Generates a Perfect Match Genomic Landscape (PMGL) from a Reference Genome Self Landscape (RGSL) and a Query Genome (sequence reads in fastq format).
###   
### OUTPUT : PMGL#_RG#_all_PMnCR_SV.tab file. The PMGL reports the RGSL (except the IDF column) plus the number of perfect match occurrences of each Reference String in 
### the Read Strings Dataset (PM), each Reference String's PM normalized by its CR (PMnCR), each Reference String's signature value (SV).
###     
###
### USAGE : perl generatePMGL_v1.0.pl -binDir /url/ -fastqFile /url/QG#.fastq -fastaRgDir /url/RG#_#.fasta -lnRgDir /url/ -outputDir /url/ -kmerLength # -pmglId PMGL# -memory #
###      
### OPTIONS :
###   -binDir : Path of directory containing all PMGL pipeline scripts.
###   -fastqFile : Path of directory containing the Query Genome sequence reads in a single fastq file. 
###			The Query Genome sequence reads must be contained in a single file in fastq format and require the following 
###			nomenclature: QG#.fastq (the # indicates the numeric component of the QG unique identifier).
###   -fastaRgDir : Path of directory containing the Reference Genome sequence per chromosome in fasta format. 
###   -lnRgDir : Path of directory containing the RGSL per chromosome. 
###   -outputDir : Path of output directory.
###   -kmerLength : Length of kmer.
###   -pmglId : Perfect Match Genomic Landscape unique identifier (PGML#, the # indicates the numeric component of the PMGL unique identifier).
###   -memory : RAM memory to be used for jellyfish execution.
###   -h or -help
###   
### DATE : 01/10/2017
###
### Requirements : 
### 1)
###   Perl's  Libraries
###     Getopt::Long;
### 2)
###  	Jellyfish (1.1.10)
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
	'fastqFile=s', 
	'fastaRgDir=s', 
	'lnRgDir=s', 
	'outputDir=s', 
	'kmerLength=i', 
	'pmglId=s', 
	'memory=i');

&readArguments();
&main();

sub main {

my $binDir =  $opts{binDir};
my $fastqFile =  $opts{fastqFile};
my $fastaRgDir =  $opts{fastaRgDir};
my $lnRgDir =  $opts{lnRgDir};
my $outputDir =  $opts{outputDir};
my $kmerLength =  $opts{kmerLength};

my $pmglId =  $opts{pmglId};
my $memory =  $opts{memory};

my @fastaFileArray;
my @lnFileArray;

#### Generating Directory

	my $outBD = $outputDir."BD/";
	my $outKmerCovPlot = $outputDir."kmerCovPlot/";
	my $outLN = $outputDir."landscape/";
	my $outPMGL = $outputDir."PMGL/";
	my $nameRGSL ="";

	system ("rm -rf ${outBD}");
	system ("rm -rf ${outKmerCovPlot}");
	system ("rm -rf ${outLN}");
	system ("rm -rf ${outPMGL}");

	mkdir $outBD;
	mkdir $outKmerCovPlot;
	mkdir $outLN;
	mkdir $outPMGL;

#### Generating jellyfish database
print "\nGenerating jellyfish database\n";

	system ("echo jellyfish count -o ${pmglId}_BD -m ${kmerLength} --both-strands -s ${memory}G -t 12 ${fastqFile} >> ${outBD}JOBS_makeDB_Jellifish.bash ");
	system ("perl ${binDir}makeSGE_v1.0.pl -i ${outBD}JOBS_makeDB_Jellifish.bash -memory ${memory} -emos");

	if ($?) {
		   exit(0);
	}

	### Executing JOB in cluster

 	print "\tExecuting JOB in cluster \n";
 	
	my $qsubBD = `qsub -t 1-1:1 ${outBD}JOBS_makeDB_Jellifish.sge`;

	$qsubBD =~ /Your job-array (\d+)\.\d+-\d+:\d+ \(\"JOBS_makeDB_Jellifish\"\) has been submitted/;
    $qsubBD = $1;

    &waitingForTheCluster($qsubBD);

#Comparing jellyfish database with Reference Genome sequence
print "Comparing jellyfish database with Reference Genome sequence\n";

	opendir(my $fastaRgDirIndex, $fastaRgDir) || die "Error :( $! \n"; 
    while(readdir $fastaRgDirIndex){  
        if (-f $fastaRgDir . "/" . $_) { 
        	if ($_ =~ /.*fasta/) { 
                     push (@fastaFileArray, $_); 
            }
        } 
    }
    closedir $fastaRgDirIndex; 

	system ("echo \\#\!/bin/bash >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash");
	system ("echo \\#\$ -cwd >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo \\#\$ -j y >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo \\#\$ -S /bin/bash >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo \\#\$ -e ${outKmerCovPlot}JOBS_kmer-cov-plot.error >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo \\#\$ -o ${outKmerCovPlot}JOBS_kmer-cov-plot.salida >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo \\#\$ -N JOBS_kmer-cov-plot >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo \\#\$ -l qname=all.q >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo . /etc/profile.d/modules.sh >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo module load amos/3.1.0 >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	system ("echo module load jellyfish/1.1.10 >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");

	$fastaFileArray[0] =~ /(.*)_\d+.fasta/;
	$nameRGSL = $1;

	foreach my $fastaFileArrayIndex (0 .. $#fastaFileArray) {  
		$fastaFileArray[$fastaFileArrayIndex] =~ /.*_(\d+).fasta/;
	    my $chromosome = $1;
		system ("echo kmer-cov-plot --jellyfish -s ${outBD}${pmglId}_BD_0 \\< ${fastaRgDir}$fastaFileArray[$fastaFileArrayIndex] \\> ${outKmerCovPlot}${pmglId}_${nameRGSL}_${chromosome}.kmer-cov-plot >> ${outKmerCovPlot}JOBS_kmer-cov-plot.bash ");
	}

	my $qsubKmerCovPlot = `qsub ${outKmerCovPlot}JOBS_kmer-cov-plot.bash`; 

 	print "\tExecuting JOB in cluster \n";
        
	$qsubKmerCovPlot =~ /Your job (\d+) \(\"JOBS_kmer-cov-plot\"\) has been submitted/;
    $qsubKmerCovPlot = $1;

    &waitingForTheCluster($qsubKmerCovPlot);

print "Concatenating RGSL with jellyfish output\n";

	opendir(my $lnRgDirIndex, $lnRgDir) || die "Error :( $! \n"; 
    while(readdir $lnRgDirIndex){ 
        if (-f $lnRgDir . "/" . $_) {
        	if ($_ =~ /.*tab/) { 
                     push (@lnFileArray, $_); 
            }
        } 
    }
    closedir $lnRgDirIndex; 

    my $numberLNReferenceFile = @lnFileArray;

	foreach my $lnFileArrayIndex (0 .. $#lnFileArray) {  
		$lnFileArray[$lnFileArrayIndex] =~ /.*_(\d+)_landscape.tab/;
	    my $chromosome = $1;

		system ("echo  perl ${binDir}jointRefKmer-cov_v1.0.pl -referenceFile ${lnRgDir}$lnFileArray[$lnFileArrayIndex] -kmerCovFile ${outKmerCovPlot}${pmglId}_${nameRGSL}_${chromosome}.kmer-cov-plot -outputFile ${outLN}${pmglId}_${nameRGSL}_${chromosome}_landscape.tab  >> ${outKmerCovPlot}JOBS_Final_LN.bash");
		
	}

	system ("perl ${binDir}makeSGE_v1.0.pl -i ${outKmerCovPlot}JOBS_Final_LN.bash -memory ${memory} ");

	print "\tExecuting JOB in cluster \n";

	my $qsubLN = `qsub -t 1-${numberLNReferenceFile}:1 ${outKmerCovPlot}JOBS_Final_LN.sge`;

	$qsubLN =~ /Your job-array (\d+)\.\d+-\d+:\d+ \(\"JOBS_Final_LN\"\) has been submitted/;
    $qsubLN = $1;

    &waitingForTheCluster($qsubLN);

#Concatenating all chromosomes in a single file
print "\nConcatenating all chromosomes in a single file\n";

	system ("cat ${outLN}${pmglId}*_landscape.tab  >> ${outPMGL}${pmglId}_${nameRGSL}_all.tab ");

#PMnCR normalization
print "\nPMnCR normalization\n";

	system ("echo  perl ${binDir}normalizedByCountReference_v1.0.pl -inputfile ${outPMGL}${pmglId}_${nameRGSL}_all.tab  > ${outPMGL}JOBS_normalized_string.bash");

	system ("perl ${binDir}makeSGE_v1.0.pl -i ${outPMGL}JOBS_normalized_string.bash -memory ${memory} ");

	my $qsubNCR = `qsub -t 1-1:1 ${outPMGL}JOBS_normalized_string.sge`;

	$qsubNCR =~ /Your job-array (\d+)\.\d+-\d+:\d+ \(\"JOBS_normalized_string\"\) has been submitted/;
    $qsubNCR = $1;

    &waitingForTheCluster($qsubNCR);

#signature value 
print "\nsignature value\n";

	system ("echo  perl ${binDir}signatureValue_v2.0.pl -inputfile ${outPMGL}${pmglId}_${nameRGSL}_all_PMnCR.tab  > ${outPMGL}JOBS_signature_value.bash");

	system ("perl ${binDir}makeSGE_v1.0.pl -i ${outPMGL}JOBS_signature_value.bash -memory ${memory} ");

	my $qsubSV = `qsub -t 1-1:1 ${outPMGL}JOBS_signature_value.sge`;

	$qsubSV =~ /Your job-array (\d+)\.\d+-\d+:\d+ \(\"JOBS_signature_value\"\) has been submitted/;
    $qsubSV = $1;

    &waitingForTheCluster($qsubSV);


print "\nPMGL created\n\n";

}

### Read arguments
sub readArguments {

	my $mandatoryParameters = 'true';

	if (!$opts{fastqFile}){
		$opts{fastqFile} = '';
	}
	if (!$opts{fastaRgDir}){
		$opts{fastaRgDir} = '';
	}
	if (!$opts{lnRgDir}){
		$opts{lnRgDir} = '';
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
	if (!$opts{pmglId}){
		$opts{pmglId} = '';
	}
	if (!$opts{memory}){
		$opts{memory} = '';
	}

	### Mandatory parameters
    if ($opts{fastqFile} eq ''){
            print ("Needs the -fastqFile parameter \n");
            $mandatoryParameters = 'false';
    }elsif (! -f $opts{fastqFile}){
        print ("\t-fastqFile File does not exist: $opts{fastqFile}\n");
        $mandatoryParameters = 'false';
    }elsif ($opts{fastqFile} !~ /QG\d+\.fastq$/){
        print ("\t-fastqFile Required the following nomenclature: QG#.fastq \n");
        $mandatoryParameters = 'false';
    }

	if ($opts{fastaRgDir} eq ''){
		print ("\tNeeds the -fastaRgDir parameter \n");
		$mandatoryParameters = 'false';
	}elsif (! -d $opts{fastaRgDir}){
		print ("\t-fastaRgDir Directory does not exist: $opts{fastaRgDir}\n");
		$mandatoryParameters = 'false';
	}

	if ($opts{lnRgDir} eq ''){
		print ("\tNeeds the -lnRgDir parameter \n");
		$mandatoryParameters = 'false';
	}elsif (! -d $opts{lnRgDir}){
		print ("\t-lnRgDir Directory does not exist: $opts{lnRgDir}\n");
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


	if ($opts{pmglId} eq ''){
		print ("\tNeeds the -pmglId parameter \n");
		$mandatoryParameters = 'false';
	}elsif($opts{pmglId} !~ /^PMGL\d+$/){
		print ("\t-pmglId Required the following nomenclature: PMGL# \n");
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
