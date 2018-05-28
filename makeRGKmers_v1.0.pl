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
###     Delfino Garcia Alonso, send comments to delfinog\@ccg.unam.mx
###     Jair Santiago Garcia Sotelo, send comments to jsgarcia\@liigh.unam.mx
###     Kim Palacios Flores, send comments to kimpalaciosflores\@gmail.com
###
### NAME : makeRGKmers_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Decomposes the Reference Genome into kmers of a specified length.
###    
###   
### OUTPUT : Fasta file containing each Reference String's sequence uniquely identified through its starting position within a chromosome. 
###     
###
### USAGE : perl makeRGKmers_v1.0.pl -inputFile /url/RG#.fasta -outputFile /url/Kmers/RG#_#_strings_kbp.fna -kmerLength # -rgslid RGSL#
###      
### OPTIONS :
###   -inputFile : The RG chromosome fasta file. 
###   -outputFile : Path of output directory.
###   -kmerLength : Length of kmer.
###   -rgslid : Reference Genome Self Landscape unique identifier (RGSL#, the # indicates the numeric component of the RGSL unique identifier).
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
    'inputFile|i=s',
    'outputFile|o=s',
    'kmerLength|k=i',
    'rgslid|p=s',
    'h|help');

my %dataHash = ();

&readArguments();
&main();

sub main {

  &fillHash($opts{inputFile});
  &printSeq;

}

sub printSeq{
   if($opts{outputFile}){
      open(OUTPUT_FILE,">$opts{outputFile}");
   }
   foreach my $key (sort keys %dataHash){
      for(my $index=0; $index <= length($dataHash{$key}); $index++){
      $addseq = substr($dataHash{$key},$index,$opts{kmerLength});
	 if($addseq =~ /^[ATCGatcg]{$opts{kmerLength},$opts{kmerLength}}$/){
	    $posSeq = $index+1;
     	    if($opts{rgslid}){
	       print OUTPUT_FILE ">${opts{rgslid}}_$posSeq\n$addseq\n";
     	    }else{
	       print OUTPUT_FILE ">${key}_$posSeq\n$addseq\n";
     	    }
	 }
      }
   }
   close(OUTPUT_FILE);
}

sub fillHash{
   my $inFile = shift;

   $/ = "\n>";

   open(FNA,"$inFile") || die "Could not open: $opts{inputFile}\n";
   while(<FNA>){
      my $query = $_;
      $query =~ s/\>//g;	
      my ($datos,@dna) = split(/\n/,$query);
      my $vSequence = "";
      foreach my $key (@dna){
         $vSequence .= $key; 
      }
      $dataHash{$datos} = $vSequence;
   }
   close(FNA);
}


### Read arguments
sub readArguments {
  my $mandatoryParameters = 'true';

  if (!$opts{inputFile}){
    $opts{inputFile} = '';
  }
  if (!$opts{kmerLength}){
    $opts{kmerLength} = '';
  }
  if (!$opts{rgslid}){
    $opts{rgslid} = '';
  }
  if (!$opts{outputFile}){
        $opts{outputFile} = '';
    }

  ### Mandatory parameters
  if ($opts{inputFile} eq ''){
    print ("Needs the -inputFile parameter \n");
    $mandatoryParameters = 'false';
  }
  if ($opts{kmerLength} eq ''){
    print ("Needs the -kmerLength parameter \n");
    $mandatoryParameters = 'false';
  } 
  if ($opts{rgslid} eq ''){
    print ("Needs the -rgslid parameter \n");
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
