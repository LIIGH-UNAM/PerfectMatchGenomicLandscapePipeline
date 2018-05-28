#!/usr/bin/perl

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
### NAME : jointRefKmer_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Concatenates the RGSL with the jellyfish output (kmer count in the Query Genome).
###   
### OUTPUT : The file containing the RGSL and the PM column of the PMGL to be generated.
###     
###
### USAGE : perl jointRefKmer_v1.0.pl -referenceFile /url/LN/ -kmerCovFile /url/kmerCovPlot/ -outputFile /url/
###      
### OPTIONS :
###   -referenceFile : The chromosome RGSL file.
###   -kmerCovFile : The chromosome jellyfish output. 
###   -outputFile : Path of the output file. 
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

use Getopt::Long;
my %opts = ();

### Parameters 
GetOptions (\%opts,
   'referenceFile=s',
   'kmerCovFile=s',
   'outputFile=s',
   'h|help'
   );

&readArguments();
&main();

sub main {

  my $referenceFile =  $opts{referenceFile};
  my $kmerCovFile =  $opts{kmerCovFile};
  my $outputFile =  $opts{outputFile};

  my @structureReferenceFile;
  my @structureKmerCovFile;

  open(INPUT_REFERENCE_FILE, $referenceFile) or die("Can't read ${referenceFile} file");
    @structureReferenceFile = <INPUT_REFERENCE_FILE>;
  close(INPUT_REFERENCE_FILE);

  foreach my $structureReferenceFileIndex (0 .. $#structureReferenceFile) {
    my $currentLine = $structureReferenceFile[$structureReferenceFileIndex];
    chop $currentLine;
    my ($id,$cr,$sequence) = split("\t", $currentLine);
    $id =~ /([\w_]+)_(\d+)$/;
    $rgslData{$2} = "$id\t$cr\t$sequence";
  }

  open(INPUT_KMERCOV, $kmerCovFile) or die("Can't read ${kmerCovFile} file");
    @structureKmerCovFile = <INPUT_KMERCOV>;
  close(INPUT_KMERCOV);

 open(OUTPUT_FILE,">${outputFile}") or die("Can't create ${outputFile} file");
  foreach my $structureKmerCovFileIndex (0 .. $#structureKmerCovFile) {
    my $currentLine = $structureKmerCovFile[$structureKmerCovFileIndex];
    chop $currentLine;
    my($idCoverage,$coverage) = split("\t", $currentLine);

    if($idCoverage =~ /^\d+$/){
      print OUTPUT_FILE "$rgslData{$idCoverage}\t$coverage\n";
    }
  }
  close(OUTPUT_FILE);
}

### Read arguments
sub readArguments {
  my $mandatoryParameters = 'true';

  if (!$opts{referenceFile}){
    $opts{referenceFile} = '';
  }
  if (!$opts{kmerCovFile}){
    $opts{kmerCovFile} = '';
  }
  if (!$opts{urlFastqFile}){
    $opts{urlFastqFile} = '';
  }

  ### Mandatory parameters
  if ($opts{referenceFile} eq ''){
    print ("Needs the -referenceFile parameter \n");
    $mandatoryParameters = 'false';
  }
  if ($opts{kmerCovFile} eq ''){
    print ("Needs the -kmerCovFile parameter \n");
    $mandatoryParameters = 'false';
  }
    if ($opts{outputFile} eq ''){
    print ("Needs the -outputFile parameter \n");
    $mandatoryParameters = 'false';
  }

  if ($mandatoryParameters eq 'false'){
      exit(0);
  }
}

