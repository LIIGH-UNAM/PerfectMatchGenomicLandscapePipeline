#!/bin/bash

###
### AUTHOR :
###     Jair Santiago Garcia Sotelo, send comments to jsgarcia\@liigh.unam.mx
###     Kim Palacios Flores, send comments to kimpalaciosflores\@gmail.com
###
### NAME : findReads_v1.0.sh
###
### DESCRIPTION :  Searches the Downstream Recovery String in the Query Genome sequence reads and executes the catAlignment_v1.0.pl script.
###
### VERSION : version 1.0
###     
### DATE : 01/10/2017
###
###	$1	Kmer used to attract sequence reads.
###	$2	Path of file containing the Query Genome sequence reads in a single fastq file. 
###	$3	File containing growing alignment fasta file. 
###	$4	Path of file containing the RGSL.
###	$5	File containing the raw RG sequence per chromosome in txt format. 
###	$6	Path of directory containing reads attracted by kmers from successive alignment extensions. 
###	$7	Path of directory containing generated Read Families.
###	$8	Path of directory containing all PMGL pipeline scripts.
###	$9	Minimum number of identical cut reads per Family.
###	$10	Path of directory containing alignments in fasta format. 
###	$11	Length of kmer.
### $12	Number of the current alignment extension. 
### $13	Maximum number of different Read Families.
### $14	Query Genome sequence to be concatenated.
###
###

grep -i "$1" $2  > $6$3_forward.txt

var1="$(echo $1 | rev )"
var1="$(echo $var1 | tr [:lower:] [:upper:])"
var1="$(echo $var1 | tr ATGC TACG)"

grep -i "$var1" $2  > $6$3_reverse.txt

cat $6$3_reverse.txt | perl -pe 'chomp;tr/ACGTacgt/TGCAtgca/;$_=reverse."\n"' > $6$3_reverseToForward.txt

cat $6$3_forward.txt > $6$3.txt
cat $6$3_reverseToForward.txt >> $6$3.txt

##### Executing the catAlignment_v1.0.pl script.

perl $8catAlignment_v2.0.pl -readsFile $6$3.txt -fastqFile $2 -sequence $1 -pmWalkId $3 -rawRgFile $5 -familyFile $7$3.txt -minCountFamily $9 -alignmentFile ${10}$3.txt -kmerLength ${11} -numberWalk ${12} -maxFamily ${13} -seqQueryAdd ${14} 
