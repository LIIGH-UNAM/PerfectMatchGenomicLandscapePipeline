#!/usr/bin/perl

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
### NAME : makeRGSLFile_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Generates the RGSL file using the bowtie output. 
###   
### OUTPUT : The file containing the ID, CR, SEQ, and IDF (optional) columns of the RGSL. 
###     
###
### USAGE : perl makeRGSLFile_v1.0.pl -inputFile /url/bowtie/RGSL#_#_strings_kbp.out  -outputFile /url/bowtie/RGSL#_#_strings_kbp.landscape -idFamily 
###      
### OPTIONS :
###   -inputFile : The file containing the bowtie output. 
###   -outputFile : Path of the output file. 
###   -idFamily (optional) : Allows generating the IDF column.
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

my %opts = ();

### Parameters 
GetOptions (\%opts,
        'inputFile|i=s',   
        'outputFile|o=s',  
        'idFamily|idF',    
        'h|help');

&readArguments();
&main();

sub main {

  my $inputFile =  $opts{inputFile};
  my $outputFile =  $opts{outputFile};

  my @structureInputFile;
  my %sequence = ();
  my %sequenceATCG = ();


  print "\n\nLoading:  ${inputFile}\n";

  open(INPUT_FILE, $inputFile) or die("Can't read ${inputFile} file");
    @structureInputFile = <INPUT_FILE>;
  close(INPUT_FILE);

  foreach my $index (0 .. $#structureInputFile) {
    my $currentLine = $structureInputFile[$index];
    chop $currentLine;
    my @currentLineArray = split("\t", $currentLine);

    my $positionMatch = $currentLineArray[3] + 1;   
    push @{ $sequence{$currentLineArray[0]} },"$currentLineArray[2]_$positionMatch";
    $sequenceATCG{$currentLineArray[0]} = $currentLineArray[4];
  }

  open(OUTPUT_FILE,">$outputFile") or die("Can't create ${outputFile} file");

  foreach my $key (keys %sequence){

     $num = scalar(@{$sequence{$key}});
     $key =~ /_(\d+)$/;
     if($opts{idFamily}){
        print OUTPUT_FILE "$1\t$key\t$num\t$sequenceATCG{$key}\t";
        my @sortedIdFamily = sort { $a <=> $b } @{$sequence{$key}};
        foreach my $idFamily (@sortedIdFamily){
           print OUTPUT_FILE "$idFamily,";
        }
        print OUTPUT_FILE "\n";
     }else{
        print OUTPUT_FILE "$1\t$key\t$num\t$sequenceATCG{$key}\n";
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
  if (!$opts{outputFile}){
    $opts{outputFile} = '';
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

  if ($mandatoryParameters eq 'false'){
    exit(0);
  } 
}

