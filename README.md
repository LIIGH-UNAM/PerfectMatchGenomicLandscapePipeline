# PerfectMatchGenomicLandscapePipeline
Precise detection of genome variation using the Perfect Match Genomic Landscape Pipeline

-------------------------------------------------
Citation
-------------------------------------------------

A Perfect Match Genomic Landscape Provides a Unified Framework for the Precise Detection of Variation in Natural and Synthetic Haploid Genomes

Kim Palacios-Flores, Jair García-Sotelo, Alejandra Castillo, Carina Uribe, Luis Aguilar, Lucía Morales, Laura Gómez-Romero, José Reyes, Alejandro Garciarubio, Margareta Boege and Guillermo Dávila

GENETICS April 1, 2018 vol. 208 no. 4 1631-1641

https://doi.org/10.1534/genetics.117.300589

PMID: 29367403 

-------------------------------------------------
AUTHORS
-------------------------------------------------

Kim Palacios-Flores
Jair Garcia-Sotelo

Contact information:

Kim Palacios Flores, kimpalaciosflores@gmail.com
Jair Santiago Garcia Sotelo, jsgarcia@liigh.unam.mx

-------------------------------------------------
PMGL PIPELINE Versions
-------------------------------------------------
Version 1.0 PMGL Pipeline description:

This pipeline is used to reveal signatures of variation embedded in a Perfect Match Genomic Landscape (PMGL), discover the underlying variants, generate a customized Reference Genome, and validate the precise location and nature of the introgressed variants. 

Version 2.0 PMGL Pipeline description:

The following corrections / improvements were made to the PMGL pipeline:

Correction:

1)	The extent of Reference Genome sequence extracted to be aligned with the Query Genome sequence exceeded the size of the Query Genome sequence for alignment extensions number 2 and onwards. This has been corrected in the current PMGL pipeline version such that both sequences are the same length.

Improvements:

1)	Chromosomes in the RGSL structure were incorrectly ordered for reference genomes with more than 19 chromosomes. This has now been corrected and all chromosomes are organized numerically. Also, a regular expression was modified such that more than 99 chromosomes can be scanned for signatures of variation.

2)	A validation step was added to the customization function such that the customization script is halted if no signatures of variation were solved in the previous PMGL pipeline module (module 5: Orchestrates the interpretation and extension of alignments. Discovers the nature of variants). 

3)	The PMGL pipeline was modified such that signature of variation sites that yield a perfect alignment (no underlying variation) can be handled. 


-------------------------------------------------
PMGL PIPELINE
-------------------------------------------------

This pipeline is used to reveal signatures of variation embedded in a Perfect Match Genomic Landscape (PMGL), discover the underlying variants, generate a customized Reference Genome, and validate the precise location and nature of the introgressed variants. The PMGL pipeline is divided into six modules:

1) Generation of a Reference Genome Self Landscape (RGSL).
2) Generation of a Perfect Match Genomic Landscape (PMGL).
3) Scanning of the PMGL using the zero-trail scan. Pinpoints the position of variants. 
4) Generation of the first alignment at each signature of variation using the MUSCLE Multiple Sequence Alignment tool.
5) Orchestrates the interpretation and extension of alignments. Discovers the nature of variants.
6) Generation of a customized Reference Genome. The customized Reference Genome sequence is validated by performing steps 1), 2), and 3) using the customized Reference Genome and the original Query Genome sequence reads.

---------------------
Cluster configuration
---------------------

System of Cluster Administration: Bright Cluster Manager 7.1
S.O. Centos 7 x86_64
Resources Scheduler: SGE (Sun Grid Engine) 2011.11p1
Environment Modules: 3.2.10
Internal Management Network: 1GbE
Internal Network: Infiniband 40Gbps

------------
Requirements
------------

R version 3.3.2
R libraries:
- stringr
- seqinr
- stringi
- readtext
- optparse

Perl version 5.16.3
Perl libraries:
- Getopt::Long;
- LWP::Simple;
- File::Copy;
- Math::Round;
- Math::BigFloat;

Bowtie 0.12.7
Jellyfish 1.1.10

Internet access is required to execute the Client of MUSCLE Multiple Sequence Alignment tool.
IMPORTANT NOTE: The Client of MUSCLE Multiple Sequence Alignment tool was modified by Jair Santiago Garcia Sotelo and Kim Palacios Flores on 02/11/2017 to re-execute MUSCLE using a new JOB ID when it is not responding. The modified Client of MUSCLE Multiple Sequence Alignment tool is integrated into the PMGL pipeline scripts. 

-------------
Generate RGSL
-------------

DESCRIPTION

Generates a Reference Genome Self Landscape (RGSL) from a Reference Genome sequence in fasta format. 

USAGE

The only script that needs to be executed by the user for this module is the following:

perl generateRGSL_v1.0.pl -binDir /url/ -fastaDir /url/RG#_#.fasta -bowtieDir /url/ -outputDir /url/ -kmerLength # -rgslId RGSL# -memory #
   
OUTPUT

RGSL#.tab file. The RGSL reports each Reference String's unique identifier (ID), the number of times its sequence is present in the entire Reference Genome (CR), its DNA sequence (SEQ), and the unique identifiers of all Reference Strings in the entire Reference Genome that share the same sequence (IDF). 

Generates four types of output directories:

bowtie: Contains the comparison files between the bowtie database and the chromosome kmers fasta files.
fasta: Contains Reference Genome fasta file.
kmers: Contains the chromosome kmers fasta files. 
LN: Contains the chromosome landscapes.
RGSL: Contains the Reference Genome Self Landscape.

Scripts comprising the Generate RGSL module:

generateRGSL_v1.0.pl
makeSGE_v1.0.pl
makeRGKmers_v1.0.pl
makeRGSLFile_v1.0.pl
orderFamily_v1.0.pl

-------------
Generate PMGL
-------------

DESCRIPTION

Generates a Perfect Match Genomic Landscape (PMGL) from a Reference Genome Self Landscape (RGSL) and a Query Genome (sequence reads in fastq format).

USAGE

The only script that needs to be executed by the user for this module is the following:

perl generatePMGL_v1.0.pl -binDir /url/ -fastqFile /url/QG#.fastq -fastaRgDir /url/RG#_#.fasta -lnRgDir /url/ -outputDir /url/ -kmerLength # -pmglId PMGL# -memory #
   
OUTPUT

PMGL#_RG#_all_PMnCR_SV.tab file. The PMGL reports the RGSL (except the IDF column) plus the number of perfect match occurrences of each Reference String in 
the Read Strings Dataset (PM), each Reference String's PM normalized by its CR (PMnCR), each Reference String's signature value (SV).

Generates four types of output directories:

BD: Contains the jellyfish database.
kmerCovPlot: Contains the comparison files between the Reference Genome chromosome fasta files and the jellyfish database.
landscape: Contains the chromosome landscapes.
PMGL: Contains the Perfect Match Genomic Landscape.

Scripts comprising the Generate PMGL module:

generatePMGL_v1.0.pl
makeSGE_v1.0.pl
jointRefKmer-cov_v1.0.pl
normalizedByCountReference_v1.0.pl
signatureValue_v1.0.pl

-----------------------------
Generate PMGL Zero Trail Scan
-----------------------------

DESCRIPTION

Locates signatures of variation along a Perfect Match Genomic Landscape (PMGL) using the zero-trail scan.

USAGE

The only script that needs to be executed by the user for this module is the following:

perl PMGLZeroTrailScan.pl -inputFile PMGL#_RG#_all_PMnCR_SV.tab -minPMn # -maxPMn_1 # -CRn_1 # -lowComplexity -zeroTrail

OUTPUT

PMGL#_RG#_all_PMnCR_SV_ZeroTrailScan.tab File. Contains the PMGL rows corresponding to the Downstream Recovery String of each signature of variation plus its zero-trail length (last column).

Scripts comprising the PMGL Zero Trail Scan module:

PMGLZeroTrailScan_V0.3.pl

------------------------
Generate First Alignment
------------------------

DESCRIPTION 

For each signature of variation located using the zero-trail scan, generates an alignment anchored by the Downstream Recovery String. 

USAGE

The only script that needs to be executed by the user for this module is the following:

perl generateAligment_v1.0.pl -binDir /url/ -outputDir /url/ -zeroTrailScanFile PMGL#_RG#_all_PMnCR_SV_ZeroTrailScan.tab -rgslFile /URL/RGSL#.tab -fastqFile /url/QG#.fastq -rawRgDir /url/ -rgId RG# -minCountFamily # -memory # -kmerLength # -maxFamily #
   
OUTPUT 

Generates four types of output directories:

alignment: Contains the formatted alignment files.
bash: Contains the script files for JOB execution.
family: Contains the different Query Genome sequences (Read Families) that contain a perfect match with the Downstream Recovery String, and that have a minimum number of occurrences in the Query Genome.
muscle: Contains the alignment files generated by the MUSCLE Multiple Sequence Alignment tool.
reads: Contains the Query Genome sequence reads containing a perfect match with the Downstream Recovery String. 

Scripts comprising the Generate Alignment module:

generateAlignment_v1.0.pl
makeSGE_v1.0.pl
findReads_v1.0.sh
catAlignment_v1.0.pl
muscle_lwp.pl

------------------------------------------
Interpretation and Extension of Alignments
------------------------------------------

DESCRIPTION 

Orchestrates iterative alignment extensions at signatures of variation to uncover the underlying variants. Signatures of variation are classified as either solved or unsolved.

USAGE

The only script that needs to be executed by the user for this module is the following:

Rscript --vanilla /url/organizer_v1.0.R --binDir /url/ --rawRgDir /url/ --readsDir /url/ --familyDir /url/ --muscleDir /url/ --alignmentDir /url/ --organizerDir /url/ --rgslFile /url/RGSL#.tab --fastqFile /url/QG#.fastq --organizerName organizer --rgId RG# --kmerLength # --minCountFamily # --maxFamily # --minAnchor # --maxPmWalk #

OUTPUT

1) A file reporting the final status (solved or unsolved plus additional information) of each signature of variation.
2) An R object containing the final status report for all signatures of variation, the final alignment report for solved signatures of variation, and a partial report for unsolved signatures of variation.

Scripts comprising the Interpretation and Extension of Alignments module:

organizer_v1.0.R

Note: It is recommended to execute this script in a computer node.

-------------
Customization
-------------

DESCRIPTION

Generates a new, customized Reference Genome using the Common Variation Motifs found for solved signatures of variation by the organizer_v1.0.R program. 

USAGE

The only script that needs to be executed by the user for this module is the following:

Rscript --vanilla /url/customization_v1.0.R --rawRgDir /url/ --custRgDir /url/ --organizerFile /url/organizer.rds --kmerLength # --rgId RG# --custRgId RG#
   
OUTPUT 

1) A file reporting all sequence changes performed by the customization process decomposed into single nucleotide changes, deletions, or insertions. 
   Sequence changes are reported in the order that they were effectuated, from the more downstream positions to the more upstream ones within each chromosome. NOTE: In the case where an alignment has invaded another alignment located further upstream, and each alignment instructs a different change to be made on the same nucleotide(s) from the Reference Genome to match the Query Genome, such nucleotide(s) is not customized. Additionally, a file reporting all uncustomized Reference Genome positions due to conflicting alignment changes is generated.

2) The new, customized Reference Genome sequence per chromosome in fasta and txt format. 

Scripts comprising the Customization module:

customization_v1.0.R

--------------------
Directory structure
--------------------
.
├── bin
│   ├── catAlignment_v1.0.pl
│   ├── customization_v1.0.R
│   ├── findReads_v1.0.sh
│   ├── generateAlignment_v1.0.pl
│   ├── generateBashAligment_v1.0.pl
│   ├── generatePMGL_v1.0.pl
│   ├── generatePMWalk_v1.0.pl
│   ├── generateRGSL_v1.0.pl
│   ├── jointRefKmer-cov_v1.0.pl
│   ├── library.pm
│   ├── makeRGKmers_v1.0.pl
│   ├── makeRGSLFile_v1.0.pl
│   ├── makeSGE_v1.0.pl
│   ├── muscle_lwp.pl
│   ├── normalizedByCountReference_v1.0.pl
│   ├── orderFamily_v1.0.pl
│   ├── organizer_v1.0.R
│   ├── PMGLZeroTrailScan_V0.3.pl
│   └── signatureValue_v1.0.pl
├── PMGL
│   └── PMGL#
│       ├── BD
│       │   ├── JOBS_makeDB_Jellifish.bash
│       │   ├── JOBS_makeDB_Jellifish.error
│       │   ├── JOBS_makeDB_Jellifish.salida
│       │   ├── JOBS_makeDB_Jellifish.sge
│       │   └── PMGL#_BD_0
│       ├── kmerCovPlot
│       │   ├── JOBS_Final_LN.bash
│       │   ├── JOBS_Final_LN.error
│       │   ├── JOBS_Final_LN.salida
│       │   ├── JOBS_Final_LN.sge
│       │   ├── JOBS_kmer-cov-plot.bash
│       │   ├── JOBS_kmer-cov-plot.salida
│       │   └── PMGL#_RG#_#.kmer-cov-plot
│       ├── landscape
│       │   └── PMGL#_RG#_#_landscape.tab
│       └── PMGL
│           ├── JOBS_normalized_string.bash
│           ├── JOBS_normalized_string.error
│           ├── JOBS_normalized_string.salida
│           ├── JOBS_normalized_string.sge
│           ├── JOBS_signature_value.bash
│           ├── JOBS_signature_value.error
│           ├── JOBS_signature_value.salida
│           ├── JOBS_signature_value.sge
│           ├── PMGL#_RG#_all_PMnCR_SV.tab
│           ├── PMGL#_RG#_all_PMnCR.tab
│           └── PMGL#_RG#_all.tab
├── PMGL_SCAN
│   └── PMGL#_RG#_all_PMnCR_SV_ZeroTrailScan.tab
├── PMGL_VARIANTS
│   └── PMGL#
│       ├── Organizer1_PMGL#_IDs_status_report.txt
│       ├── Organizer1_PMGL#.rds
│       ├── alignment
│       │   ├── RGSL#_#_#_F#.fasta
│       │   ├── RGSL#_#_#_F#_Query.fasta
│       │   ├── RGSL#_#_#_F#_W#F#.fasta
│       │   ├── RGSL#_#_#_F#_W#F#_Query.fasta
│       ├── bash
│       │   ├── RG1_alignment.bash
│       │   ├── RG1_alignment.error
│       │   ├── RG1_alignment.salida
│       │   └── RG1_alignment.sge
│       ├── family
│       │   ├── RGSL#_#_#_F#.fasta
│       │   ├── RGSL#_#_#_family.txt
│       │   ├── RGSL#_#_#_F#_W#F#.fasta
│       │   └── RGSL#_#_#_F#_W#_family.txt
│       ├── muscle
│       │   ├── RGSL#_#_#_F#.fasta.aln-clustalw.clw
│       │   ├── RGSL#_#_#_F#.fasta.out.txt
│       │   ├── RGSL#_#_#_F#.fasta.phylotree.ph
│       │   ├── RGSL#_#_#_F#.fasta.pim.pim
│       │   ├── RGSL#_#_#_F#.fasta.sequence.txt
│       │   ├── RGSL#_#_#_F#_W#F#.fasta.aln-clustalw.clw
│       │   ├── RGSL#_#_#_F#_W#F#.fasta.out.txt
│       │   ├── RGSL#_#_#_F#_W#F#.fasta.phylotree.ph
│       │   ├── RGSL#_#_#_F#_W#F#.fasta.pim.pim
│       │   └── RGSL#_#_#_F#_W#F#.fasta.sequence.txt
│       └── reads
│           ├── RGSL#_#_#_F#_W#_forward.txt
│           ├── RGSL#_#_#_F#_W#_reverseToForward.txt
│           ├── RGSL#_#_#_F#_W#_reverse.txt
│           ├── RGSL#_#_#_F#_W#.txt
│           ├── RGSL#_#_#_forward.txt
│           ├── RGSL#_#_#_reverseToForward.txt
│           ├── RGSL#_#_#_reverse.txt
│           └── RGSL#_#_#.txt
├── QG
│   └── QG#
│       └── QG#.fastq
├── RG
│   └── RG#
│    	├── RG#_to_RG#_customization_report.txt
│       ├── RG#_#.fasta
│       └── RG#_#.txt
├── RGSL
│   └── RGSL#
│       ├── bowtie
│       │   ├── JOBS_makeLandscape.bash
│       │   ├── JOBS_makeLandscape.error
│       │   ├── JOBS_makeLandscape.salida
│       │   ├── JOBS_makeLandscape.sge
│       │   ├── JOBS_Order_family.bash
│       │   ├── JOBS_Order_family.error
│       │   ├── JOBS_Order_family.salida
│       │   ├── JOBS_Order_family.sge
│       │   ├── JOBS_sortLandscape.bash
│       │   ├── RGSL#_#_strings_kbp.landscape
│       │   ├── RGSL#_#_strings_kbp.landscape_sort
│       │   └── RGSL#_#_strings_kbp.out
│       ├── fasta
│       │   ├── JOBS_bowtie-build.bash
│       │   ├── JOBS_bowtie-build.error
│       │   ├── JOBS_bowtie-build.salida
│       │   ├── JOBS_bowtie-build.sge
│       │   ├── RG#_allChr.fasta
│       │   ├── RG#_allChr.fasta.1.ebwt
│       │   ├── RG#_allChr.fasta.2.ebwt
│       │   ├── RG#_allChr.fasta.3.ebwt
│       │   ├── RG#_allChr.fasta.4.ebwt
│       │   ├── RG#_allChr.fasta.rev.1.ebwt
│       │   └── RG#_allChr.fasta.rev.2.ebwt
│       ├── Kmers
│       │   ├── JOBS_bowtieALL_RGSL.bash
│       │   ├── JOBS_bowtieALL_RGSL.error
│       │   ├── JOBS_bowtieALL_RGSL.salida
│       │   ├── JOBS_bowtieALL_RGSL.sge
│       │   └── RG#_#_strings_kbp.fna
│       ├── LN
│       │   └── RGSL#_#_landscape.tab
│       └── RGSL
│           └── RGSL#.tab
└── tree.txt

