#!/usr/bin/perl -w

use strict;
use LWP::Simple;
use File::Copy;
use Getopt::Long;
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
### NAME : catAlignment_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Generates and filters Read Families (using the minCountFamily and maxFamily parameters) to be aligned against the Reference Genome. Concatenates growing alignments.
###   
### OUTPUT : The alignment fasta file per Read Family. 
###     
###
### USAGE : perl catAlignment_v1.0.pl -readsFile /url/reads/RGSL#_#_#_F#.fasta -sequence ATGC -pmWalkId RGSL#_#_#_F#.fasta -rawRgFile /url/RG#.txt -fastqFile /url/QG#.fastq -familyFile /url/family/RGSL#_#_#_F#.fasta -minCountFamily # -alignmentFile /url/alignment/RGSL#_#_#_F#.fasta -numberWalk # -seqQueryAdd ATGC -kmerLength # -maxFamily #
###      
### OPTIONS :
###   -readsFile : File containing reads attracted by the kmer. 
###   -sequence : Kmer used to attract sequence reads.
###   -pmWalkId : Unique identifier of growing alignment fasta file. 
###   -rawRgFile : The file containing the raw Reference Genome sequence per chromosome in txt format. 
###   -fastqFile :  The file containing the Query Genome sequence reads in a single fastq file. 
###   -familyFile : The output file for the generated Read families. 
###   -minCountFamily : Minimum number of identical cut reads per Family.
###   -alignmentFile : The output file containing the growing Query Genome and Reference Genome sequences to be aligned using the MUSCLE tools. 
###   -numberWalk : Number of the current alignment extension. 
###   -seqQueryAdd : Query Genome sequence to be concatenated.
###   -kmerLength : Length of kmer.
###   -maxFamily : Maximum number of different Read Families.
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
    'readsFile=s', 
    'sequence=s', 
    'pmWalkId=s', 
    'rawRgFile=s', 
    'fastqFile=s', 
    'familyFile=s', 
    'minCountFamily=s', 
    'alignmentFile=s', 
    'numberWalk=i', 
    'seqQueryAdd=s', 
    'kmerLength=i', 
    'maxFamily=s'
    );


&readArguments();
&main();

sub main { 

    my $rawRgFile = $opts{rawRgFile};
    my $pmWalkId = $opts{pmWalkId};
    my $readsFile =  $opts{readsFile};
    my $sequence = $opts{sequence};
    my $fastqFile = $opts{fastqFile};
    my $familyFile = $opts{familyFile};
    my $minCountFamily = $opts{minCountFamily};
    my $alignmentFile = $opts{alignmentFile};
    my $numberWalk = $opts{numberWalk};
    my $seqQueryAdd = $opts{seqQueryAdd};
    my $kmerLength =  $opts{kmerLength};
    my $maxFamily =  $opts{maxFamily};

	$fastqFile =~ s/\/.*\///;
	$fastqFile =~ s/\.fastq//;
    $pmWalkId  =~ /^.*_\d+_(\d+).*/;
    my $posGenome = $1;

	print "\n\nLoading:  ${readsFile}\n";

	my @structureReadsFile;
	open(INPUT_READS_FILE, $readsFile) or die("Can't read ${readsFile} file");
                @structureReadsFile = <INPUT_READS_FILE>;
    close(INPUT_READS_FILE);

	my $toMuscleFile = $familyFile;
	$toMuscleFile =~ s/.txt//;
    if ($toMuscleFile =~ /_\d+$/){
        $toMuscleFile = $toMuscleFile."_";
    }

    my $querySequenceFile = $alignmentFile;
    $querySequenceFile =~ s/.txt//;
    if ($querySequenceFile =~ /_\d+$/){
        $querySequenceFile = $querySequenceFile."_";
    }

    my $killFile = $querySequenceFile;
    if ($querySequenceFile !~ /_$/){
        $killFile = $killFile."_";
    }

    my $familyFileName = $toMuscleFile;
    if ($familyFileName !~ /_$/){
        $familyFileName = $familyFileName."_";
    }

	my %familyFound = ();
    my $numberNucleotides = 25;
    foreach my $structureReadsFileIndex (0 .. $#structureReadsFile) {
		my $currentLine = $structureReadsFile[$structureReadsFileIndex];
		chop $currentLine;
		if ($currentLine !~ /^#/ && $currentLine !~ /^<id>/ && $currentLine !~ /^<sequence>/ && $currentLine !~ /^<sequence_n-1>/) {
      	    if ($currentLine =~ /([ATGC]{$numberNucleotides})$sequence/i){
				if ($familyFound{$1}){
                    $familyFound{$1} = $familyFound{$1} + 1;
                }else{
                    $familyFound{$1} = 1;
                }
        	}
		}
    }

    my $posGenomeAbsolutoLeft = 0;

    if ($numberWalk == 1){
        $posGenomeAbsolutoLeft = ($posGenome - ($numberWalk*$numberNucleotides))-1;
    }
    else{
        $posGenomeAbsolutoLeft = ($posGenome - (($numberWalk*$numberNucleotides)- (10*($numberWalk-1))))-1;
    }

    my $posGenomeAbsolutoRight = $posGenome + ($kmerLength-1);
    my $referenceSeq="";

    my $countFamily = 1;

    if($posGenomeAbsolutoLeft > 0) {

        open(INPUT_RAW_FILE,$rawRgFile) || die "Sorry, I couldn't open the file $rawRgFile\n";
                seek INPUT_RAW_FILE, $posGenomeAbsolutoLeft, 1;
                read INPUT_RAW_FILE, $referenceSeq, $posGenomeAbsolutoRight-$posGenomeAbsolutoLeft;
                seek INPUT_RAW_FILE, 0, 0;
        close(INPUT_RAW_FILE);

        ### Creating all Read Families file
        open(OUTPUT_FAMILY_FILE_TXT, "> ${familyFileName}family.txt") or die("Can't create file");

        	print OUTPUT_FAMILY_FILE_TXT "Query Genome: ".$fastqFile."\n";
        	print OUTPUT_FAMILY_FILE_TXT "Recovery String ID: ".$opts{pmWalkId}."\n";
        	print OUTPUT_FAMILY_FILE_TXT "Recovery String Sequence: ".$sequence."\n";
        	print OUTPUT_FAMILY_FILE_TXT "Reference Sequence: ".$referenceSeq."\n\n";
        	print OUTPUT_FAMILY_FILE_TXT "Sequence Reads Families \n\n";

            if ($seqQueryAdd ne ""){
                $sequence = $sequence.$seqQueryAdd;
            }
              
            my @validFamily;

        	foreach my $querySequence (reverse sort { $familyFound{$a} <=> $familyFound{$b} } keys %familyFound) {
        		if ($familyFound{$querySequence} >= $minCountFamily){
                    push(@validFamily, $querySequence.$sequence);

                	print OUTPUT_FAMILY_FILE_TXT $querySequence.$sequence."\t".$familyFound{$querySequence}."\n\n";
                    print OUTPUT_FAMILY_FILE_TXT ">Reference\n".$referenceSeq."\n";
                    print OUTPUT_FAMILY_FILE_TXT ">Query\n".$querySequence.$sequence."\n\n";
                    print OUTPUT_FAMILY_FILE_TXT "\n";

                    $countFamily++;
               	} 
                last if ($countFamily > $maxFamily);
            }  
        close(OUTPUT_FAMILY_FILE_TXT);  

        my $numberOfValidFamily = @validFamily;

        if ($numberOfValidFamily < $maxFamily){
            foreach my $index (0 .. $#validFamily) {  
                ### Creating alignment file
                my $numberOfFamily = $index+1;
                open(OUTPUT_TO_MUSCLE_FILE, "> ${toMuscleFile}F".$numberOfFamily.".fasta") or die("Can't create ${toMuscleFile}F".$numberOfFamily.".fasta file");
                    print OUTPUT_TO_MUSCLE_FILE ">Reference\n".$referenceSeq."\n";
                    print OUTPUT_TO_MUSCLE_FILE ">Query\n".$validFamily[$index];
                close(OUTPUT_TO_MUSCLE_FILE);             
                ### Creating file containing growing Query Genome Sequence
                open(OUTPUT_QUERY_SEQUENCE_FILE, "> ${querySequenceFile}F".$numberOfFamily."_Query.fasta") or die("Can't create ${querySequenceFile}F".$numberOfFamily.".fasta file");  
                    print OUTPUT_QUERY_SEQUENCE_FILE $validFamily[$index]."\n";
                close(OUTPUT_QUERY_SEQUENCE_FILE);  
            }
        }else{
                ### Creating kill file for exceeds Read Families case
                open(OUTPUT_KILL_EXCEED_FAMILIES_FILE, "> ${killFile}killExceedFamilies.txt") or die("Can't create file");  
                    print OUTPUT_KILL_EXCEED_FAMILIES_FILE "kill \n";
                close(OUTPUT_KILL_EXCEED_FAMILIES_FILE);  
        }
    }
    if ($countFamily == 1){
        ### Creating kill file for no Read Families case
        open(OUTPUT_KILL_NO_FAMILIES_FILE, "> ${killFile}killNoFamilies.txt") or die("Can't create file");  
            print OUTPUT_KILL_NO_FAMILIES_FILE "kill \n";
        close(OUTPUT_KILL_NO_FAMILIES_FILE);
    }
}

### Read arguments
sub readArguments {
	my $mandatoryParameters = 'true';

	if (!$opts{sequence}){
        $opts{sequence} = '';
    }
	if (!$opts{readsFile}){
		$opts{readsFile} = '';
	}
	if (!$opts{pmWalkId}){
        $opts{pmWalkId} = '';
    }
	if (!$opts{rawRgFile}){
        $opts{rawRgFile} = '';
    }
	if (!$opts{fastqFile}){
        $opts{fastqFile} = '';
    }
	if (!$opts{familyFile}){
        $opts{familyFile} = '';
    }
	if (!$opts{minCountFamily}){
        $opts{minCountFamily} = '';
    }
    if (!$opts{alignmentFile}){
        $opts{alignmentFile} = '';
    }
    if (!$opts{numberWalk}){
        $opts{numberWalk} = '';
    }
    if (!$opts{kmerLength}){
        $opts{kmerLength} = '';
    }
    if (!$opts{seqQueryAdd}){
        $opts{seqQueryAdd} = '';
    }
    if (!$opts{maxFamily}){
        $opts{maxFamily} = '';
    }

    ### Mandatory parameters
	if ($opts{sequence} eq ''){
        print ("Needs the -sequence parameter \n");
        $mandatoryParameters = 'false';
    }
	if ($opts{readsFile} eq ''){
		print ("Needs the -readsFile parameter \n");
		$mandatoryParameters = 'false';
	}
	if ($opts{pmWalkId} eq ''){
        print ("Needs the -pmWalkId parameter \n");
        $mandatoryParameters = 'false';
    }
	if ($opts{rawRgFile} eq ''){
        print ("Needs the -rawRgFile parameter \n");
        $mandatoryParameters = 'false';
    }
	if ($opts{fastqFile} eq ''){
        print ("Needs the -fastqFile parameter \n");
        $mandatoryParameters = 'false';
    }
    if ($opts{familyFile} eq ''){
        print ("Needs the -familyFile parameter \n");
        $mandatoryParameters = 'false';
    }
    if ($opts{minCountFamily} eq ''){
        print ("Needs the -minCountFamily parameter \n");
        $mandatoryParameters = 'false';
    }
    if ($opts{alignmentFile} eq ''){
        print ("Needs the -alignmentFile parameter \n");
        $mandatoryParameters = 'false';
    }
    if ($opts{numberWalk} eq ''){
        print ("Needs the -numberWalk parameter \n");
        $mandatoryParameters = 'false';
    }
    if ($opts{kmerLength} eq ''){
        print ("Needs the -kmerLength parameter \n");
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
