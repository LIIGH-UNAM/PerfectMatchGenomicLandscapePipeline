#!/usr/bin/env Rscript


### AUTHOR :
###     Jair Santiago Garcia Sotelo, send comments to jsgarcia@liigh.unam.mx
###     Kim Palacios Flores, send comments to kimpalaciosflores@gmail.com
###
### NAME : customization_v1.0.R
###
### VERSION : version 1.0
###
### DESCRIPTION : Generates a new, customized Reference Genome using the Common Variation Motifs found for solved signatures of variation by the organizer_v1.0.R program. 
###   
### OUTPUT : 
### 1) A file reporting all sequence changes performed by the customization process decomposed into single nucleotide changes, deletions, or insertions. 
###    Sequence changes are reported in the order that they were effectuated, from the more downstream positions to the more upstream ones within each chromosome. 
###    NOTE: In the case where an alignment has invaded another alignment located further upstream, and each alignment instructs a different change to be made
###    on the same nucleotide(s) from the Reference Genome to match the Query Genome, such nucleotide(s) is not customized. Additionally, a file reporting all uncustomized 
###    Reference Genome positions due to conflicting alignment changes is generated.
### 2) The new, customized Reference Genome sequence per chromosome in fasta and txt format. 
###
### USAGE : Rscript --vanilla /url/customization_v1.0.R --rawRgDir /url/ --custRgDir /url/ --organizerFile /url/organizer.rds --kmerLength # --rgId RG# --custRgId RG#
###      
### DATE : 11/16/2017
###
### Requirements : 
###     R version 3.3.2
###     R libraries:
###       -stringr
###       -stringi
###       -readtext
###       -seqinr
###       -optparse



library(stringr)
library(stringi)
library(readtext)
library(seqinr)
library(optparse)

option_list = list(
  make_option(c("-d","--rawRgDir"),type="character",default=NULL,
              help="Path of directory containing the corresponding raw RG sequence per chromosome in txt format",metavar="character"),
  
  make_option(c("-o", "--custRgDir"), type="character", default=NULL, 
              help="Path of directory for resulting customization report and customized RG", metavar="character"),
  
  make_option(c("-i","--organizerFile"),type="character",default=NULL,
              help="Path of organizer R object file",metavar="character"),
  
  make_option(c("-k","--kmerLength"),type="integer",default=NULL,
              help="Length of kmer",metavar="integer"),
  
  make_option(c("-r","--rgId"),type="character",default=NULL,
              help="RG ID",metavar="character"),
  
  make_option(c("-c", "--custRgId"), type="character", default=NULL, 
              help="RG ID for resulting customized RG", metavar="character")
);

opt_parser=OptionParser(option_list=option_list);
opt=parse_args(opt_parser);


missingMandatoryParameter<-FALSE

if(is.null(opt$rawRgDir)){
  print("rawRgDir argument must be supplied",quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$rawRgDir))){
  print(paste0("--rawRgDir Directory does not exist: ",opt$rawRgDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if(is.null(opt$custRgDir)){
  print("custRgDir argument must be supplied",quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(dir.exists(opt$custRgDir))){
  print(paste0("--custRgDir Directory does not exist: ",opt$custRgDir), quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if(is.null(opt$organizerFile)){
  print("organizerFile argument must be supplied",quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(file_test("-f",opt$organizerFile))){
  print(paste0("--organizerFile File does not exist: ",opt$organizerFile), quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if(is.null(opt$kmerLength)){
  print("kmerLength argument must be supplied",quote=FALSE)
  missingMandatoryParameter<-TRUE
}

if(is.null(opt$rgId)){
  print("rgId argument must be supplied",quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(grepl("^RG\\d+$",opt$rgId))){
  print("--rgId requires the following nomenclature: RG#")
  missingMandatoryParameter<-TRUE
}

if(is.null(opt$custRgId)){
  print("custRgId argument must be supplied",quote=FALSE)
  missingMandatoryParameter<-TRUE
}else if(!(grepl("^RG\\d+$",opt$custRgId))){
  print("--custRgId requires the following nomenclature: RG#")
  missingMandatoryParameter<-TRUE
}


if(missingMandatoryParameter){
  print_help(opt_parser)
  stop("Missing mandatory parameter")
}


### DESCRIPTION: Generates a table of sequence changes to be performed by the Customize function.
ObtainRgChanges<-function(organizer,kmerLength,rgId,custRgId,custRgDir){
  if(length(organizer$solvedIds)==0){
    stop("There are no solved signatures of variation. The Reference Genome customization will not be performed.", call. = FALSE)
  }else{
    solved<-organizer$solvedIds
    allCommonVariationMotifs<-lapply(solved,function(Id){
      return(Id$commonVariationMotif)
    })
    positions<-lapply(allCommonVariationMotifs,function(commonVariationMotif){
      return(commonVariationMotif$position)
    })
    reference<-lapply(allCommonVariationMotifs,function(commonVariationMotif){
      return(commonVariationMotif$referenceNucleotide)
    })
    chromosomes<-character()
    ord<-numeric()
    for(i in seq(1,length(positions))){
      IdPosition<-as.numeric(strsplit(names(positions)[i],split="_")[[1]][3])
      chromosomes<-append(chromosomes,rep(strsplit(names(positions)[i],split="_")[[1]][2],length(positions[[i]])))
      if(!("-" %in% reference[[i]])){
        ord<-append(ord,rep(1,length(positions[[i]])))
      }else{
        substract<-numeric(length(positions[[i]]))
        for(j in seq(1,length(positions[[i]]))){
          substract[j]<-length(which(reference[[i]][1:j]=="-"))
        }
        positions[[i]]<-positions[[i]]-substract
        for(j in seq(1,length(unique(positions[[i]])))){
          countSamePos<-length(which(positions[[i]]==unique(positions[[i]])[j]))
          ord<-append(ord,seq(1,countSamePos))
        }
      }
      positions[[i]]<-IdPosition-(positions[[i]]-kmerLength)
    }
    positions<-unlist(unname(positions))
    reference<-unlist(unname(reference))
    query<-lapply(allCommonVariationMotifs,function(commonVariationMotif){
      return(commonVariationMotif$queryNucleotide)
    })
    query<-unlist(unname(query))
    
    rgChanges<-data.frame("chr"=chromosomes,"RG_pos"=positions,"ord"=ord,"Reference"=reference,"Query"=query,stringsAsFactors=FALSE)
    rgChanges<-rgChanges[order(as.numeric(rgChanges$chr),rgChanges$RG_pos,-rgChanges$ord,decreasing=TRUE),]
    rgChanges<-unique(rgChanges)
    rgChangesPos<-unique(rgChanges$RG_pos)
    rgChangesPerPos<-lapply(rgChangesPos,function(pos){
      return(rgChanges[which(rgChanges$RG_pos==pos),])
    })
    names(rgChangesPerPos)<-rgChangesPos
    conflictingChanges<-vapply(rgChangesPerPos,function(posChanges){
      if(dim(unique(posChanges[,c("chr","RG_pos","ord","Reference")]))[1]==1 & length(unique(posChanges$Query))>1){
        return(FALSE)
      }else{
        return(TRUE)
      }
    },logical(1),USE.NAMES=TRUE)
    eliminatePos<-names(conflictingChanges)[which(!(conflictingChanges))]
    conflictingChangesFilepath<-paste0(custRgDir,paste0(paste(rgId,custRgId,sep="_to_"),"_conflicting_changes_report.txt"))
    if(file.exists(conflictingChangesFilepath)){
      file.remove(conflictingChangesFilepath)
    }
    if(length(eliminatePos)>0){
      for(i in seq(1,length(eliminatePos))){
        cat(paste("Due to conflicting alignment changes, the following RG position will not be customized: ",eliminatePos[i]),file=conflictingChangesFilepath,sep="\n",append=TRUE)
      }
      print(paste(paste0(paste(rgId,custRgId,sep="_to_"),"_conflicting_changes_report.txt"),"has been generated"))
    }
    rgChanges<-rgChanges[which(!(rgChanges$RG_pos %in% eliminatePos)),]
    return(rgChanges)
  }
}


### DESCRIPTION: Introgresses the sequence changes to generate a customized Reference Genome. 
Customize<-function(rgChanges,rawRgDir,rgId,custRgId,custRgDir){
  nucleotides<-c("A","a","G","g","T","t","C","c")
  chromsToBeCustomized<-unique(rgChanges$chr)
  customizedChromsSeq<-lapply(chromsToBeCustomized,function(chrom){
    chromRaw<-readtext(file=paste0(rawRgDir,paste0(paste(rgId,chrom,sep="_"),".txt")))
    chromSeq<-chromRaw$text
    chromRgChanges<-rgChanges[which(rgChanges$chr==chrom),]
    message<-sprintf("Customizing chromosome %s",chrom)
    print(message,quote=FALSE)
    for(i in seq(1,dim(chromRgChanges)[1])){
      if(chromRgChanges$Reference[i] %in% nucleotides & chromRgChanges$Query[i] %in% nucleotides){
        if(toupper(substr(chromSeq,chromRgChanges$RG_pos[i],chromRgChanges$RG_pos[i]))==toupper(chromRgChanges$Reference[i])){
          stri_sub(chromSeq,chromRgChanges$RG_pos[i],chromRgChanges$RG_pos[i])<-chromRgChanges$Query[i]
        }else{
          stop("Unexpected nucleotide in raw RG sequence")
        }
      }else if(chromRgChanges$Reference[i] %in% nucleotides & chromRgChanges$Query[i]=="-"){
        if(toupper(substr(chromSeq,chromRgChanges$RG_pos[i],chromRgChanges$RG_pos[i]))==toupper(chromRgChanges$Reference[i])){
          stri_sub(chromSeq,chromRgChanges$RG_pos[i],chromRgChanges$RG_pos[i])<-""
        }else{
          stop("Unexpected nucleotide in raw RG sequence")
        }
      }else if(chromRgChanges$Reference[i]=="-" & chromRgChanges$Query[i] %in% nucleotides){
        stri_sub(chromSeq,chromRgChanges$RG_pos[i],chromRgChanges$RG_pos[i]-1)<-chromRgChanges$Query[i]
      }
    }
    chromSeq<-toupper(chromSeq)
    cat(chromSeq,file=paste0(custRgDir,paste0(paste(custRgId,chrom,sep="_"),".txt")))
    write.fasta(chromSeq,as.string=TRUE,names=paste(paste(custRgId,chrom,sep="_"),paste("Customized from",rgId)),file.out=paste0(custRgDir,paste0(paste(custRgId,chrom,sep="_"),".fasta")))
    return(chromSeq)
  })
  names(customizedChromsSeq)<-chromsToBeCustomized
  allRgChroms<-str_extract(str_extract(list.files(path=rawRgDir,pattern=paste0("^",paste0(rgId,"_\\d+\\.txt$"))),"_\\d+\\."),"\\d+")
  intactChroms<-allRgChroms[which(!(allRgChroms %in% rgChanges$chr))]
  if(length(intactChroms)!=0){
    intactRgFiles<-paste0(paste(rgId,intactChroms,sep="_"),".txt")
    file.copy(paste0(rawRgDir,intactRgFiles),paste0(custRgDir,paste0(paste(custRgId,intactChroms,sep="_"),".txt")))
    for(i in seq(1,length(intactChroms))){
      intactChromRaw<-readtext(file=paste0(rawRgDir,intactRgFiles[i]))
      intactChromSeq<-intactChromRaw$text
      intactChromSeq<-toupper(intactChromSeq)
      write.fasta(intactChromSeq,as.string=TRUE,names=paste(paste(custRgId,intactChroms[i],sep="_"),paste("Customized from",rgId)),file.out=paste0(custRgDir,paste0(paste(custRgId,intactChroms[i],sep="_"),".fasta")))
    }
  }
  write.table(rgChanges,file=paste0(custRgDir,paste0(paste(rgId,custRgId,sep="_to_")),"_customization_report.txt"),sep="\t",quote=FALSE,row.names=FALSE)
  print(paste(paste0(paste(rgId,custRgId,sep="_to_"),"_customization_report.txt"),"has been generated"))
  return(customizedChromsSeq)
}



organizer<-readRDS(file=opt$organizerFile)
rgChanges<-ObtainRgChanges(organizer,opt$kmerLength,opt$rgId,opt$custRgId,opt$custRgDir)
customizedChromsSeq<-Customize(rgChanges,opt$rawRgDir,opt$rgId,opt$custRgId,opt$custRgDir)

