#!/usr/bin/perl

use Getopt::Long;
my %opts = ();

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
###     Kim Palacios Flores, send comments to kimpalaciosflores\@gmail.com
###
### NAME : makeSGE_v1.0.pl
###
### VERSION : version 1.0
###
### DESCRIPTION : Generates the SGE file for executing the bash file in the cluster. 
###   
### OUTPUT : The SGE file.
###     
###
### USAGE : perl makeSGE_v1.0.pl -inputFile /url/JOBS_.bash -emos -memory #  
###      
### OPTIONS :
###   -inputFile : The bash file. 
###   -queue : The cluster queue. 
###   -emos (optional): Allows loading the amos and jellyfish modules. 
###   -memory : RAM memory to be used for executing the queue.
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

### Parameters 
GetOptions (\%opts,
      'inputFile|i=s',  
      'queue|q=s',            
      'emos|e',           
      'memory|m=i',       
      'help|h');

&readArguments();
&main();

sub main {
  my $inputFile =  $opts{inputFile};
  my $queue =  $opts{queue};
  my $emos =  $opts{emos};
  my $memory =  $opts{memory};


  $inputFile =~ /(.*\/)(.*)\.(\w{4,4})$/;
  $path = $1;
  $name = $2;

  open(OUTPUT_FILE,">$path/${name}.sge") or die("Can't create ${name}.sge file");

    print OUTPUT_FILE "\#\$ -S /bin/bash\n";
    print OUTPUT_FILE "\#\$ -cwd\n";
    if($queue){
       print OUTPUT_FILE "\#\$ -l qname=${queue}\n";
    }else{
       print OUTPUT_FILE "\#\$ -l qname=all.q\n";
    }

    if($emos){
       print OUTPUT_FILE ". /etc/profile.d/modules.sh\n";
       print OUTPUT_FILE "module load amos/3.1.0\n"; 
       print OUTPUT_FILE "module load jellyfish/1.1.10\n"; 
    }
    print OUTPUT_FILE "\#\$ -wd ${path}\n";
    print OUTPUT_FILE "\#\$ -e ${path}${name}.error\n";
    print OUTPUT_FILE "\#\$ -o ${path}${name}.salida\n";
    print OUTPUT_FILE "\#\$ -N $name\n";
    print OUTPUT_FILE "\#\$ -l virtual_free=${memory}G\n";
    print OUTPUT_FILE "source /etc/bashrc\n";
    print OUTPUT_FILE "SEEDFILE=$path${name}.bash\n";
    print OUTPUT_FILE "SEED=\$(cat \$SEEDFILE | head -n \$SGE_TASK_ID | tail -n 1)\n";
    print OUTPUT_FILE "\$SEED\n";

  close(OUTPUT_FILE);
}

### Read arguments
sub readArguments {

  my $mandatoryParameters = 'true';

  if (!$opts{inputFile}){
    $opts{inputFile} = '';
  }
  if (!$opts{memory}){
    $opts{memory} = '';
  }

  ### Mandatory parameters
  if ($opts{inputFile} eq ''){
    print ("Needs the -inputFile parameter \n");
    $mandatoryParameters = 'false';
  }
  if ($opts{memory} eq ''){
    print ("Needs the -memory parameter \n");
    $mandatoryParameters = 'false';
  }
  if ($mandatoryParameters eq 'false'){
    exit(0);
  }
}
