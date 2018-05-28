#!/usr/bin/env Rscript

### AUTHOR :
###     Jair Santiago Garcia Sotelo, send comments to jsgarcia@liigh.unam.mx
###     Kim Palacios Flores, send comments to kimpalaciosflores@gmail.com
###
### NAME : organizer_v1.0.R
###
### VERSION : version 1.0
###
### DESCRIPTION : Orchestrates iterative alignment extensions at signatures of variation to uncover the underlying variants. Signatures of variation are classified as either solved or unsolved.
###   
### OUTPUT : 
### 1) A file reporting the final status (solved or unsolved plus additional information) of each signature of variation.
### 2) An R object containing the final status report for all signatures of variation, the final alignment report for solved signatures of variation, and a partial report for unsolved signatures of variation.
###
### USAGE : Rscript --vanilla /url/organizer_v1.0.R --binDir /url/ --rawRgDir /url/ --readsDir /url/ --familyDir /url/ --muscleDir /url/ --alignmentDir /url/ --organizerDir /url/ --rgslFile /url/RGSL#.tab --fastqFile /url/QG#.fastq --organizerName organizer --rgId RG# --kmerLength # --minCountFamily # --maxFamily # --minAnchor # --maxPmWalk #
###      
### DATE : 11/16/2017
###
### Requirements : 
###     R version 3.3.2
###     R libraries:
###       -stringr
###       -seqinr
###       -stringi
###       -optparse


library(stringr)
library(seqinr)
library(stringi)
library(optparse)


option_list = list(
  make_option(c("-b","--binDir"),type="character",default=NULL,
              help="Path of directory containing all PMGL pipeline scripts",metavar="character"),
  
  make_option(c("-d","--rawRgDir"),type="character",default=NULL,
              help="Path of directory containing the raw RG sequence per chromosome in txt format",metavar="character"),
  
  make_option(c("-e","--readsDir"),type="character",default=NULL,
              help="Path of directory containing reads attracted by kmers from successive alignment extensions",metavar="character"),
  
  make_option(c("-i","--familyDir"),type="character",default=NULL,
              help="Path of directory containing generated Read Families",metavar="character"),
  
  make_option(c("-m","--muscleDir"),type="character",default=NULL,
              help="Path of directory containing MUSCLE alignments",metavar="character"),
  
  make_option(c("-j", "--alignmentDir"), type="character", default=NULL, 
              help="Path of directory containing alignments in fasta format", metavar="character"),
  
  make_option(c("-l","--organizerDir"),type="character",default=NULL,
              help="Path of directory where organizer R object will be saved",metavar="character"),
  
  make_option(c("-n","--rgslFile"),type="character",default=NULL,
              help="Path of file containing the RGSL",metavar="character"),
  
  make_option(c("-q","--fastqFile"),type="character",default=NULL,
              help="Path of file containing the sequencing reads in fastq format",metavar="character"),
  
  make_option(c("-o","--organizerName"),type="character",default=NULL,
              help="Name of organizer",metavar="character"),
  
  make_option(c("-r","--rgId"),type="character",default=NULL,
              help="Reference Genome unique identifier",metavar="character"),
  
  make_option(c("-k","--kmerLength"),type="integer",default=NULL,
              help="Length of kmer",metavar="integer"),
  
  make_option(c("-c","--minCountFamily"),type="integer",default=NULL,
              help="Minimum number of identical cut reads per Family",metavar="integer"),
  
  make_option(c("-f","--maxFamily"),type="integer",default=10,
              help="Maximum number of different Read Families [default= %default]",metavar="integer"),
  
  make_option(c("-a","--minAnchor"),type="integer",default=20,
              help="Minimum length for perfect match anchor on extending side of alignment [default= %default]",metavar="integer"),
  
  make_option(c("-w","--maxPmWalk"),type="integer",default=3,
              help="Maximum number of alignment extensions [default= %default]",metavar="integer")
);

opt_parser=OptionParser(option_list=option_list);
opt=parse_args(opt_parser);


missingMandatoryParameter<-FALSE

if (is.null(opt$binDir)){
  print("binDir argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$binDir))){
  print(paste0("--binDir Directory does not exist: ",opt$binDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if (is.null(opt$rawRgDir)){
  print("rawRgDir argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$rawRgDir))){
  print(paste0("--rawRgDir Directory does not exist: ",opt$rawRgDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if (is.null(opt$readsDir)){
  print("readsDir argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$readsDir))){
  print(paste0("--readsDir Directory does not exist: ",opt$readsDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}


if (is.null(opt$familyDir)){
  print("familyDir argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$familyDir))){
  print(paste0("--familyDir Directory does not exist: ",opt$familyDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if (is.null(opt$muscleDir)){
  print("muscleDir argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$muscleDir))){
  print(paste0("--muscleDir Directory does not exist: ",opt$muscleDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}


if (is.null(opt$alignmentDir)){
  print("alignmentDir argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$alignmentDir))){
  print(paste0("--alignmentDir Directory does not exist: ",opt$alignmentDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if (is.null(opt$organizerDir)){
  print("organizerDir argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$organizerDir))){
  print(paste0("--organizerDir Directory does not exist: ",opt$organizerDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if (is.null(opt$rgslFile)){
  print("rgslFile argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(file_test("-f",opt$rgslFile))){
  print(paste0("--rgslFile File does not exist: ",opt$rgslFile), quote=FALSE)
  missingMandatoryParameter<-TRUE
}else{
  fileName<-basename(opt$rgslFile)
  if(!(grepl("^RGSL\\d+.tab$",fileName))){
    print("--rgslFile requires the following nomenclature: RGSL#.tab", quote=FALSE)
    missingMandatoryParameter<-TRUE
  }
}

if (is.null(opt$fastqFile)){
  print("fastqFile argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(file_test("-f",opt$fastqFile))){
  print(paste0("--fastqFile File does not exist: ",opt$fastqFile), quote=FALSE)
  missingMandatoryParameter<-TRUE
}else{
  fileName<-basename(opt$fastqFile)
  if(!(grepl("^QG\\d+.fastq$",fileName))){
    print("--fastqFile requires the following nomenclature: QG#.fastq", quote=FALSE)
    missingMandatoryParameter<-TRUE
  }
}


if (is.null(opt$organizerName)){
  print("organizerName argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if (is.null(opt$rgId)){
  print("rgId argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(grepl("^RG\\d+$",opt$rgId))){
  print("--rgId requires the following nomenclature: RG#")
  missingMandatoryParameter<-TRUE
}

if (is.null(opt$kmerLength)){
  print("kmerLength argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if (is.null(opt$minCountFamily)){
  print("minCountFamily argument must be supplied", quote=FALSE)
  missingMandatoryParameter<-TRUE
}


if(missingMandatoryParameter){
  print_help(opt_parser)
  stop("Missing mandatory parameter")
}


### DESCRIPTION: For all initial or extended alignments generated for a signature of variation, determines if:
### 1) They share a Common Variation Motif relative to the corresponding region of the Reference Genome.
### 2) A single alignment contains the Query Genome Sequence that differs from the corresponding region of the Reference Genome only at the Common Variation Motif. 
### 3) If conditions 1 and 2 are met, determines if the chosen alignment is anchored, both upstream and downstream, by perfect matches to the Reference Genome.


ChooseAlignment<-function(alignmentDir,pmWalkId,kmerLength,minAnchor){
  filenamesAlignments<-list.files(path=alignmentDir,pattern=paste0("^",paste0(pmWalkId,"\\d+\\.fasta$")))
  numAlignments<-length(filenamesAlignments)
  allAlignments<-vector("list",numAlignments)
  lengthMaxAlignment<-0
  for(i in seq(1,numAlignments)){
    alignment<-read.alignment(file=paste0(alignmentDir,filenamesAlignments[i]),format="fasta")
    
    alignment$seq[[1]]<-stri_reverse(alignment$seq[[1]])
    alignment$seq[[2]]<-stri_reverse(alignment$seq[[2]])
    alignment$seq[[3]]<-stri_reverse(alignment$seq[[3]])
    
    allAlignments[[i]]<-alignment
    lengthAlignment<-nchar(alignment$seq[[1]])
    if(lengthAlignment>lengthMaxAlignment){
      lengthMaxAlignment<-lengthAlignment
    }
  }
  referenceReport<-matrix(data="NA",nrow=numAlignments,ncol=lengthMaxAlignment)
  queryReport<-matrix(data="NA",nrow=numAlignments,ncol=lengthMaxAlignment)
  for(i in seq(1,numAlignments)){
    alignment<-allAlignments[[i]]
    mismatches<-gregexpr(pattern=" {1}",alignment$seq[[3]])
    mismatchPositions<-unlist(mismatches)
    if(!(-1 %in% mismatchPositions)){
      for(j in seq(1,length(mismatchPositions))){
        referenceChar<-substr(alignment$seq[[1]],mismatchPositions[j],mismatchPositions[j])
        queryChar<-substr(alignment$seq[[2]],mismatchPositions[j],mismatchPositions[j])
        referenceReport[i,mismatchPositions[j]]<-referenceChar
        queryReport[i,mismatchPositions[j]]<-queryChar
      }
    }
  }
  
  commonPositions<-numeric()
  commonReference<-character()
  commonQuery<-character()
  
  flagStopCommon<-0
  for(i in seq(1,lengthMaxAlignment)){
    queryCharsPerPosition<-unique(queryReport[,i])
    if(length(queryCharsPerPosition)==1){
      if(queryCharsPerPosition!="NA"){
        referenceCharsPerPosition<-unique(referenceReport[,i])
        if(length(referenceCharsPerPosition)==1){
          if(referenceCharsPerPosition!="NA"){
            if(flagStopCommon==0){
              commonPositions<-append(commonPositions,i)
              commonReference<-append(commonReference,referenceCharsPerPosition)
              commonQuery<-append(commonQuery,queryCharsPerPosition)
            }
          }
        }
      }
    }else{
      flagStopCommon<-1
    }
  }
  
  if(length(commonPositions)>0){
    print("A common variation motif has been found",quote=FALSE)
    commonVariationMotif<-list(commonPositions,commonReference,commonQuery)
    names(commonVariationMotif)<-c("position","referenceNucleotide","queryNucleotide")
  }else{
    print("A common variation motif across alignments has not been found",quote=FALSE)
    report<-list(referenceReport,queryReport)
    names(report)<-c("referenceReport","queryReport")
    return(report)
  }
  
  chooseQuerySeq<-rep("NA",lengthMaxAlignment)
  chooseQuerySeq[commonPositions]<-commonQuery
  referenceSeq<-rep("NA",lengthMaxAlignment)
  referenceSeq[commonPositions]<-commonReference
  
  chooseAlignment<-numeric()
  for(i in seq(1,numAlignments)){
    if(all(queryReport[i,]==chooseQuerySeq)){
      if(all(referenceReport[i,]==referenceSeq)){
        chooseAlignment<-append(chooseAlignment,i)
      }
    }
  }
  
  
  if(length(chooseAlignment)==1){
    print("A unique alignment containing only the common variation motif has been found",quote=FALSE)
  }else{
    print("A unique alignment containing only the common variation motif has not been found",quote=FALSE)
    report<-list(referenceReport,queryReport,commonVariationMotif)
    names(report)<-c("referenceReport","queryReport","commonVariationMotif")
    return(report)
  }
  
  comparison<-strsplit(allAlignments[[chooseAlignment]]$seq[[3]],split="")[[1]]
  anchorLength<-length(comparison)-max(which(comparison==" "))
  drsAnchorLength<-min(which(comparison==" "))-1
  if(anchorLength>=minAnchor & drsAnchorLength==kmerLength){
    missingAnchor<-FALSE
  }else{
    missingAnchor<-TRUE
  }
  
  filenameAlignment<-list.files(path=alignmentDir,pattern=paste0("^",paste0(pmWalkId,paste0(chooseAlignment,"\\.fasta$"))))
  alignment<-allAlignments[[chooseAlignment]]
  
  alignmentReport<-list(filenameAlignment,alignment$seq)
  names(alignmentReport)<-c("filenameAlignment","alignment")
  
  report<-list(alignmentReport,referenceReport,queryReport,commonVariationMotif,missingAnchor)
  names(report)<-c("alignmentReport","referenceReport","queryReport","commonVariationMotif","missingAnchor")
  return(report)
}


### DESCRIPTION: For all signatures of variation, performs the requested number of alignment extensions. 
### Using the output from the ChooseAlignment function and the number of generated Read Families at each iteration, classifies signatures of variation as solved, unsolved, or transiting. 
### Only the alignments corresponding to transiting signatures of variation enter the next round of alignment extension. If no solution is found before the maximum number of alignment extensions is
### reached, the transiting signatures of variation are classified as unsolved. 
### For solved signatures of variation, the final alignment is reported (the sequence orientation is reversed such that the alignment begins with the Downstream Recovery String), as well as the Common Variation Motif (which is used by the customization_v1.R script to customize the Reference Genome).
### For unsolved signatures of variation, the result generated by the ChooseAlignment function in the iteration prior to the signature of variation being classified as unsolved is reported. 

Organizer<-function(alignmentDir,kmerLength,minAnchor,maxPmWalk,statusIds,transitingIds,solvedIds,unsolvedIds,alignmentDirString,maxFamilyString,fastqFileString,rgslFileString,rawRgDirString,kmerLengthString,binDirString,rgIdString,readsDirString,familyDirString,minCountFamilyString,muscleDirString){
  print("IDs Status Report:",quote=FALSE)
  print(statusIds)
  if(sum(statusIds$solved)+sum(statusIds$unsolved)==dim(statusIds)[1]){
    report<-list(statusIds,solvedIds,unsolvedIds)
    names(report)<-c("statusIds","solvedIds","unsolvedIds")
    return(report)
  }else{
    for(i in seq(1,dim(statusIds)[1])){
      if(statusIds$initial[i]){
        killExceedFamilies<-list.files(path=alignmentDir,pattern=paste0(paste0("^",rownames(statusIds)[i]),"_killExceedFamilies.txt$"))
        killNoFamilies<-list.files(path=alignmentDir,pattern=paste0(paste0("^",rownames(statusIds)[i]),"_killNoFamilies.txt$"))
        if(length(killExceedFamilies)==0 & length(killNoFamilies)==0){
          message<-sprintf("Iteration %d for ID %s",statusIds$pmWalk[i],rownames(statusIds)[i])
          print(message,quote=FALSE)
          idReport<-ChooseAlignment(alignmentDir,paste0(rownames(statusIds)[i],"_F"),kmerLength,minAnchor)
          if("alignmentReport" %in% names(idReport)){
            statusIds$commonVariationMotif[i]<-TRUE
            statusIds$commonVariationMotifOnly[i]<-TRUE
            if(!(idReport$missingAnchor)){
              statusIds$anchor[i]<-TRUE
              statusIds$solved[i]<-TRUE
              solvedIds[[length(solvedIds)+1]]<-idReport
              names(solvedIds)[length(solvedIds)]<-rownames(statusIds)[i]
              message<-sprintf("ID %s has been solved",rownames(statusIds)[i])
              print(message,quote=FALSE)
            }else{
              transitingIds[[i]]<-idReport
            }
          }else if("commonVariationMotif" %in% names(idReport)){
            statusIds$commonVariationMotif[i]<-TRUE
            transitingIds[[i]]<-idReport
          }else{
            transitingIds[[i]]<-idReport
          }
        }else{
          statusIds$unsolved[i]<-TRUE
          unsolvedIds[[length(unsolvedIds)+1]]<-NA
          names(unsolvedIds)[length(unsolvedIds)]<-rownames(statusIds)[i]
          if(length(killExceedFamilies)!=0){
            statusIds$exceedsMaxReadFamilies[i]<-TRUE
            message<-sprintf("ID %s has not been solved because too many Read Families have been generated",rownames(statusIds)[i])
            print(message,quote=FALSE)
          }else if(length(killNoFamilies)!=0){
            statusIds$noReadFamilies[i]<-TRUE
            message<-sprintf("ID %s has not been solved because no Read Families have been generated",rownames(statusIds)[i])
            print(message,quote=FALSE)
          }
        }
        statusIds$initial[i]<-FALSE
      }
      else if(all(!(statusIds$unsolved[i]) & !(statusIds$solved[i]))){
        if(statusIds$pmWalk[i]<maxPmWalk){
          flagKillExceedFamilies<-FALSE
          flagKillNoFamilies<-FALSE
          filenamesFamilies<-list.files(path=alignmentDir,pattern=paste(paste0(paste0(paste0("^",rownames(statusIds)[i]),"_F\\d+"),paste0(paste(rep("_W\\d+F\\d+",statusIds$pmWalk[i]),collapse=""),"\\.fasta$"))))
          for(j in seq(1,length(filenamesFamilies))){
            familyFileString<-paste("-familyFile",filenamesFamilies[j],sep=" ")
            cmd<-paste("perl",paste0(opt$binDir,"generatePMWalk_v1.0.pl"),familyFileString,alignmentDirString,maxFamilyString,fastqFileString,rgslFileString,rawRgDirString,kmerLengthString,binDirString,rgIdString,readsDirString,familyDirString,minCountFamilyString,muscleDirString)
            system(cmd)
            killExceedFamilies<-list.files(path=alignmentDir,pattern=paste0(paste0(paste0(paste0("^",rownames(statusIds)[i]),"_F\\d+"),paste(rep("_W\\d+F\\d+",statusIds$pmWalk[i]),collapse="")),"_W\\d+_killExceedFamilies.txt$"))
            killNoFamilies<-list.files(path=alignmentDir,pattern=paste0(paste0(paste0(paste0("^",rownames(statusIds)[i]),"_F\\d+"),paste(rep("_W\\d+F\\d+",statusIds$pmWalk[i]),collapse="")),"_W\\d+_killNoFamilies.txt$"))
            if(length(killExceedFamilies)!=0){
              flagKillExceedFamilies<-TRUE
              break
            }else if(length(killNoFamilies)!=0){
              flagKillNoFamilies<-TRUE
              break
            }
          }
          if(!flagKillExceedFamilies & !flagKillNoFamilies){
            statusIds$pmWalk[i]<-statusIds$pmWalk[i]+1
            message<-sprintf("Iteration %d for ID %s",statusIds$pmWalk[i],rownames(statusIds)[i])
            print(message,quote=FALSE)
            pmWalkId<-paste0(paste0(rownames(statusIds)[i],"_F\\d+"),paste0(paste(rep("_W\\d+F\\d+",statusIds$pmWalk[i]-1),collapse=""),"_W\\d+F"))
            idReport<-ChooseAlignment(alignmentDir,pmWalkId,kmerLength,minAnchor)
            if("alignmentReport" %in% names(idReport)){
              statusIds$commonVariationMotif[i]<-TRUE
              statusIds$commonVariationMotifOnly[i]<-TRUE
              if(!(idReport$missingAnchor)){
                statusIds$anchor[i]<-TRUE
                statusIds$solved[i]<-TRUE
                solvedIds[[length(solvedIds)+1]]<-idReport
                names(solvedIds)[length(solvedIds)]<-rownames(statusIds)[i]
                message<-sprintf("ID %s has been solved",rownames(statusIds)[i])
                print(message,quote=FALSE)
              }else{
                transitingIds[[i]]<-idReport
              }
            }else if("commonVariationMotif" %in% names(idReport)){
              statusIds$commonVariationMotif[i]<-TRUE
              statusIds$commonVariationMotifOnly[i]<-FALSE
              transitingIds[[i]]<-idReport
            }else{
              statusIds$commonVariationMotif[i]<-FALSE
              statusIds$commonVariationMotifOnly[i]<-FALSE
              transitingIds[[i]]<-idReport
            }
          }else{
            statusIds$unsolved[i]<-TRUE
            unsolvedIds[[length(unsolvedIds)+1]]<-transitingIds[[i]]
            names(unsolvedIds)[length(unsolvedIds)]<-rownames(statusIds)[i]
            if(flagKillExceedFamilies){
              statusIds$exceedsMaxReadFamilies[i]<-TRUE
              message<-sprintf("ID %s has not been solved because too many Read Families have been generated",rownames(statusIds)[i])
              print(message,quote=FALSE)
            }else if(flagKillNoFamilies){
              statusIds$noReadFamilies[i]<-TRUE
              message<-sprintf("ID %s has not been solved because no Read Families have been generated",rownames(statusIds)[i])
              print(message,quote=FALSE)
            }  
          }
        }else{
          statusIds$unsolved[i]<-TRUE
          unsolvedIds[[length(unsolvedIds)+1]]<-transitingIds[[i]]
          names(unsolvedIds)[length(unsolvedIds)]<-rownames(statusIds)[i]
          message<-sprintf("ID %s has not been solved",rownames(statusIds)[i])
          print(message,quote=FALSE)
        }
      }
    }
    Organizer(alignmentDir,kmerLength,minAnchor,maxPmWalk,statusIds,transitingIds,solvedIds,unsolvedIds,alignmentDirString,maxFamilyString,fastqFileString,rgslFileString,rawRgDirString,kmerLengthString,binDirString,rgIdString,readsDirString,familyDirString,minCountFamilyString,muscleDirString)
  }
}


system(paste("rm -rf",paste0(opt$alignmentDir,"*_W*")))
system(paste("rm -rf",paste0(opt$familyDir,"*_W*")))
system(paste("rm -rf",paste0(opt$muscleDir,"*_W*")))
system(paste("rm -rf",paste0(opt$readsDir,"*_W*")))

IDs_alignment<-unique(str_extract(list.files(path=opt$alignmentDir,pattern="^[A-z0-9]+_[0-9]+_[0-9]+_F\\d+\\.fasta$"),"^[A-z0-9]+_[0-9]+_[0-9]+"))
IDs_killExceedFamilies<-str_extract(list.files(path=opt$alignmentDir,pattern="^[A-z0-9]+_[0-9]+_[0-9]+_killExceedFamilies.txt$"),"^[A-z0-9]+_[0-9]+_[0-9]+")
IDs_killNoFamilies<-str_extract(list.files(path=opt$alignmentDir,pattern="^[A-z0-9]+_[0-9]+_[0-9]+_killNoFamilies.txt$"),"^[A-z0-9]+_[0-9]+_[0-9]+")
IDs<-c(IDs_alignment,IDs_killExceedFamilies,IDs_killNoFamilies)
if(length(IDs)==0){
  stop("No alignments have been generated for any ID")
}
statusIds<-data.frame("initial"=rep(TRUE,length(IDs)),"commonVariationMotif"=rep(FALSE,length(IDs)),"commonVariationMotifOnly"=rep(FALSE,length(IDs)),"anchor"=rep(FALSE,length(IDs)),"solved"=rep(FALSE,length(IDs)),"unsolved"=rep(FALSE,length(IDs)),"pmWalk"=rep(0,length(IDs)),"exceedsMaxReadFamilies"=rep(FALSE,length(IDs)),"noReadFamilies"=rep(FALSE,length(IDs)),row.names=IDs,stringsAsFactors=FALSE)
solvedIds<-list()
unsolvedIds<-list()
transitingIds<-vector("list",length(IDs))


organizerReport<-Organizer(opt$alignmentDir,opt$kmerLength,opt$minAnchor,opt$maxPmWalk,statusIds,transitingIds,solvedIds,unsolvedIds,paste("-alignmentDir",opt$alignmentDir,sep=" "),paste("-maxFamily",opt$maxFamily,sep=" "),paste("-fastqFile",opt$fastqFile,sep=" "),paste("-rgslFile",opt$rgslFile,sep=" "),paste("-rawRgDir",opt$rawRgDir,sep=" "),paste("-kmerLength",opt$kmerLength,sep=" "),paste("-binDir",opt$binDir,sep=" "),paste("-rgId",opt$rgId,sep=" "),paste("-readsDir",opt$readsDir,sep=" "),paste("-familyDir",opt$familyDir,sep=" "),paste("-minCountFamily",opt$minCountFamily,sep=" "),paste("-muscleDir",opt$muscleDir,sep=" "))
write.table(organizerReport$statusIds,file=paste0(opt$organizerDir,paste0(opt$organizerName,"_IDs_status_report.txt")),sep="\t",quote=FALSE)
print(paste(paste0(opt$organizerName,"_IDs_status_report.txt"),"has been generated"),quote=FALSE)
saveRDS(organizerReport,file=paste0(opt$organizerDir,paste0(opt$organizerName,".rds")))
print(paste(paste0(opt$organizerName,".rds"),"has been saved"),quote=FALSE)