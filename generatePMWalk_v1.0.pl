#!/usr/bin/perl -w

#use strict;
use LWP::Simple;
use File::Copy;
use Getopt::Long;
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
### NAME : generatePMWalk_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Extends existing alignments.
###   
### OUTPUT : Generates four types of output:
### 1) The Query Genome sequence reads containing a perfect match with the corresponding kmer. 
### 2) The different Query Genome sequences (Read Families) that contain a perfect match with the corresponding kmer, plus 25 nucleotides upstream of the Query Genome, and that have a minimum number of 
###     occurrences in the Query Genome (determined by the minCountFamily parameter).
### 3) The alignment files generated by the MUSCLE Multiple Sequence Alignment tool.
### 4) Alignment files in fasta format.
###     
###
### USAGE : perl generatePMWalk_v1.0.pl -binDir /url/ -readsDir /url/reads -familyDir /url/family/ -muscleDir /url/muscle/ -alignmentDir /url/alignment/ -familyFile RGSL#_#_#_F#.fasta  -rgslFile /url/RGSL#.tab -fastqFile /url/QG#.fastq -rawRgDir /url/RG#.txt -kmerLength # -rgId RG# -minCountFamily # -maxFamily #
###      
### OPTIONS :
###   -binDir : Path of directory containing all PMGL pipeline scripts.
###   -readsDir : Path of directory containing reads attracted by kmers from successive alignment extensions. 
###   -familyDir : Path of directory containing generated Read Families.
###   -muscleDir : Path of directory containing MUSCLE alignments. 
###   -alignmentDir : Path of directory containing alignments in fasta format. 
###   -familyFile : File containing growing alignment fasta file. 
###   -rgslFile : Path of file containing the RGSL.
###   -fastqFile : Path of file containing the Query Genome sequence reads in a single fastq file. 
###   -rawRgDir : Path of directory containing the raw RG sequence per chromosome in txt format. 
###   -kmerLength : Length of kmer.
###   -rgId : Reference Genome unique identifier.
###   -minCountFamily : Minimum number of identical cut reads per Family.
###   -maxFamily : Maximum number of different Read Families.
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
    'binDir=s', 
    'readsDir=s',
    'familyDir=s', 
    'muscleDir=s', 
    'alignmentDir=s', 
    'familyFile=s', 
    'rgslFile=s',
    'fastqFile=s', 
    'rawRgDir=s', 
    'kmerLength=i', 
    'rgId=s', 
    'minCountFamily=i', 
    'maxFamily=s',
    'h|help');

&readArguments();


my @structura;
&main();

sub main {

my $binDir =  $opts{binDir};
my $readsDir =  $opts{readsDir};
my $familyDir =  $opts{familyDir};
my $muscleDir =  $opts{muscleDir};
my $alignmentDir =  $opts{alignmentDir};
my $familyFile =  $opts{familyFile};
my $rgslFile =  $opts{rgslFile};
my $fastqFile =  $opts{fastqFile};
my $rawRgDir =  $opts{rawRgDir};
my $kmerLength =  $opts{kmerLength};
my $rgId =  $opts{rgId};
my $minCountFamily =  $opts{minCountFamily};
my $maxFamily =  $opts{maxFamily};




### Generating bash: search Query Genome sequence reads containing kmer
print "\nGenerating bash: search Query Genome sequence reads containing kmer\n";

    $familyFile =~ s/.fasta//;

    open(IN, ${alignmentDir}.${familyFile}."_Query.fasta") or die("Can't read file");
        @structuraFile = <IN>;
    close(IN);

     my $nameWalkFile = "";

    foreach my $position (0 .. $#structuraFile) {
        my $current_row = $structuraFile[$position];
        chop $current_row;

        $current_row =~ /\w{10}(\w{25})(\w+)/;
    
        my $sequenceCore = $1;
        my $seqQueryAdd = $2;
        $nameWalkFile = $familyFile;

        my $chromosome = $familyFile;
        $chromosome =~ /^.*_(\d+)_(\d+)_.*/;
        $chromosome = $1;
        my $numWalk = 0;

        if ($nameWalkFile =~ /F\d+_W\d+F\d+/) {
            $nameWalkFile =~ /.*W(\d+)F(\d+)$/;
            my $walk = $1;
            my $family = $2;
            $walk = $walk+1;
            $nameWalkFile = $nameWalkFile."_W${walk}";
            $numWalk = $walk+1;
        }else{
            $nameWalkFile = $nameWalkFile."_W1";
            $numWalk  = 2;
        }
		system ("bash ${binDir}findReads_v1.0.sh ${sequenceCore} ${fastqFile} ${nameWalkFile} ${rgslFile} ${rawRgDir}${rgId}_".$chromosome.".txt ${readsDir} ${familyDir} ${binDir} ${minCountFamily} ${alignmentDir} ${kmerLength} ${numWalk} ${maxFamily} ${seqQueryAdd} ");
    }

### Executing MUSCLE tools
print "\nExecuting MUSCLE tools\n";

    my @familyFileArray;

    opendir(my $recorrer_dir_family, $familyDir) || die "Error :( $! \n";  
    while(readdir $recorrer_dir_family){ 
        if (-f $familyDir . "/" . $_) { 
            if ($_ =~ /^${nameWalkFile}.*\.fasta/) {
                     push (@familyFileArray, $_); 
            }
        } 
    }
    closedir $recorrer_dir_family; 

    for (my $familyFileArrayIndex = 0; $familyFileArrayIndex < @familyFileArray;  $familyFileArrayIndex++) { 

        print "\n\nFile:  $familyFileArray[$familyFileArrayIndex]\n";        
        system ("perl ${binDir}muscle_lwp.pl ${familyDir}$familyFileArray[$familyFileArrayIndex] --email jsgarciasotelo\@gmail.com --outfile ${muscleDir}$familyFileArray[$familyFileArrayIndex]");

        if (! -f $muscleDir.$familyFileArray[$familyFileArrayIndex].".aln-clustalw.clw"){
            $familyFileArrayIndex --;
        }

    }


### Formatting alignment files
print "\nFormatting alignment files\n";

    my @muscleDirFile;

    opendir(my $muscleDirIndex, ${muscleDir}) || die "Error :( $! \n";  
    while(readdir $muscleDirIndex){  
        if (-f ${muscleDir} . "/" . $_) { 
            if ($_ =~ /^${nameWalkFile}.*\.aln-clustalw\.clw/) { 
                     push (@muscleDirFile, $_);
            }
        } 
    }
    closedir $muscleDirIndex; 

     &cutAlignment($muscleDir, $alignmentDir, @muscleDirFile);
}








### Read arguments
sub readArguments {

        my $mandatoryParameters = 'true';

        if (!$opts{alignmentDir}){
                $opts{alignmentDir} = '';
        }

        if (!$opts{familyFile}){
                $opts{familyFile} = '';
        }

        if (!$opts{fastqFile}){
                $opts{fastqFile} = '';
        }

        if (!$opts{rgslFile}){
                $opts{rgslFile} = '';
        }

        if (!$opts{rawRgDir}){
                $opts{rawRgDir} = '';
        }

        if (!$opts{binDir}){
                $opts{binDir} = '';
        }

        if (!$opts{kmerLength}){
                $opts{kmerLength} = '';
        }

        if (!$opts{rgId}){
                $opts{rgId} = '';
        }

        if (!$opts{readsDir}){
                $opts{readsDir} = '';
        }

        if (!$opts{familyDir}){
                $opts{familyDir} = '';
        }


        if (!$opts{minCountFamily}){
                $opts{minCountFamily} = '';
        }

        if (!$opts{muscleDir}){
                $opts{muscleDir} = '';
        }

        if (!$opts{maxFamily}){
                $opts{maxFamily} = '';
        }

        ### Mandatory parameters
        if ($opts{alignmentDir} eq ''){
                print ("Needs the -alignmentDir parameter \n");
                $mandatoryParameters = 'false';
        }

        if ($opts{familyFile} eq ''){
                print ("Needs the -familyFile parameter \n");
                $mandatoryParameters = 'false';
        }

        if ($opts{fastqFile} eq ''){
                print ("Needs the -fastqFile parameter \n");
                $mandatoryParameters = 'false';
        }
        if ($opts{rgslFile} eq ''){
                print ("Needs the -rgslFile parameter \n");
                $mandatoryParameters = 'false';
        }
        if ($opts{rawRgDir} eq ''){
                print ("Needs the -rawRgDir parameter \n");
                $mandatoryParameters = 'false';
        }

        if ($opts{binDir} eq ''){
                print ("Needs the -binDir parameter \n");
                $mandatoryParameters = 'false';
        }

        if ($opts{kmerLength} eq ''){
                print ("Needs the -kmerLength parameter \n");
                $mandatoryParameters = 'false';
        }

        if ($opts{rgId} eq ''){
                print ("Needs the -rgId parameter \n");
                $mandatoryParameters = 'false';
        }
 
        if ($opts{readsDir} eq ''){
                print ("Needs the -readsDir parameter \n");
                $mandatoryParameters = 'false';
        }
        if ($opts{familyDir} eq ''){
                print ("Needs the -familyDir parameter \n");
                $mandatoryParameters = 'false';
        }

        if ($opts{minCountFamily} eq ''){
                print ("Needs the -minCountFamily parameter \n");
                $mandatoryParameters = 'false';
        }

        if ($opts{muscleDir} eq ''){
                print ("Needs the -muscleDir parameter \n");
                $mandatoryParameters = 'false';
        }

        if ($opts{maxFamily} eq ''){
                print ("Needs the -maxFamily parameter \n");
                $mandatoryParameters = 'false';
        }


        if ($mandatoryParameters eq 'false'){
                exit(0);
        }
}
