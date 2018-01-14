use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(cutAlignment numberOfLines waitingForTheCluster);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

###
### AUTHOR :
###     Jair Santiago Garcia Sotelo, send comments to jsgarcia\@liigh.unam.mx
###     Kim Palacios Flores, send comments to kimpalaciosflores\@gmail.com
###
### NAME : library.pm
###
### VERSION : version 1.0
###     
### DATE : 01/10/2017
###
###


### Trims a MUSCLE alignment
sub cutAlignment {

    my ($muscleDir, $alignmentDir, @muscleDirFile) = @_;
    
    foreach my $muscleDirFileIndex (0 .. $#muscleDirFile) {  

        my @structureMuscleDirFile;

        print "\n\nLoading:  $muscleDirFile[$muscleDirFileIndex]\n";  
        open(IN, ${muscleDir}.$muscleDirFile[$muscleDirFileIndex]) or die("Can't read $muscleDirFile[$muscleDirFileIndex] file");
            @structureMuscleDirFile = <IN>;
        close(IN);

        my $contendFile = "";
        my $reference = "";
        my $query = "";
        my $comparison = "";

        foreach my $structureMuscleDirFileIndex (1 .. $#structureMuscleDirFile) {  
            my $previousLine = $structureMuscleDirFile[$structureMuscleDirFileIndex-1];
            my $currentLine = $structureMuscleDirFile[$structureMuscleDirFileIndex];


            chop $currentLine;    
            if ($currentLine =~ /Reference/) {
                my $referenceAux  = $currentLine;
                $referenceAux =~ s/^Reference       //;
                $reference = $reference.$referenceAux;
            }
            if ($currentLine =~ /Query/) {
                my $queryAux  = $currentLine;                
                $queryAux =~ s/^Query           //;
                $query = $query.$queryAux;
            }
            #if ($currentLine =~ /\*/ ) {
            if ($previousLine =~ /Query/ ) {
                my $comparisonAux  = $currentLine;                
                $comparisonAux =~ s/^                //;
                $comparison = $comparison.$comparisonAux;
            }
            $contendFile =  $contendFile.$currentLine;
        }
        my @referenceArray = split("", $reference);
        my @queryArray = split("", $query);
        my @comparisonArray = split("", $comparison);

        my $cut = "false";
        my $new_reference = "";
        my $new_query = "";
        my $new_comparison = "";

        foreach my $referenceArrayIndex (0 .. $#referenceArray) {  
            if ($referenceArray[$referenceArrayIndex] =~ /A|T|G|C/i and $queryArray[$referenceArrayIndex] =~ /A|T|G|C/i){
                $cut = "true";
            }
            ##Guardar String 
            if ($cut eq "true"){
                $new_reference = $new_reference.$referenceArray[$referenceArrayIndex];
                $new_query = $new_query.$queryArray[$referenceArrayIndex];
                if (defined $comparisonArray[$referenceArrayIndex] ){
                    $new_comparison = $new_comparison.$comparisonArray[$referenceArrayIndex];
                }
                
            }
        }

        my $trimmedAlignment = $muscleDirFile[$muscleDirFileIndex];
        $trimmedAlignment =~ s/.aln-clustalw.clw//;
        
        open(OUT_FILE_ALIGNMNET, ">".$alignmentDir.$trimmedAlignment) or die("Can't create file");
            print OUT_FILE_ALIGNMNET ">Reference\n";
            print OUT_FILE_ALIGNMNET $new_reference."\n";
            print OUT_FILE_ALIGNMNET ">Query\n";  
            print OUT_FILE_ALIGNMNET $new_query."\n";
            print OUT_FILE_ALIGNMNET ">Comparison\n";  
            print OUT_FILE_ALIGNMNET $new_comparison."\n";
        close(OUT_FILE_ALIGNMNET);
    }

    return "true";
}

### Obtains number of lines in file
sub numberOfLines {
    my $inputFile = shift;
    open my $fh, q[<], $inputFile or return -1;
    my @currentLine = <$fh>;
    close $fh;
    
    return scalar @currentLine;
}

### Waits for JOB completion in cluster
sub waitingForTheCluster {
    my $qsub = shift;

    print "\tWaiting job ".$qsub." !!!\n";

    while (1) {
        my $endTime = time() + 20;
        my $toSleep = $endTime - time();
        sleep $toSleep;

        my $qstat = `qstat`;

        last if ($qstat !~ /.*${qsub}.*/);
    }

}
