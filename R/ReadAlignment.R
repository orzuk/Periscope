#This alphabet contains all letters from the Single-Letter Amino Acid Code (see ?AMINO_ACID_CODE) plus "*" (the stop letter), "-" (the gap letter), "+" (the hard masking letter), and "." (the not a letter or not available letter). It is stored in the AA_ALPHABET predefined constant (character vector).
library("entropy")
library("Biostrings")
setwd("D:\\Yair\\Proteins\\")
Alignments=readAAMultipleAlignment(filepath="CPD_ARATH_O04147_1-181_e-4_m30_PF07823_jackhmmer_e-4_m30_complete_run.fa")
singlescount=consensusMatrix(Alignments,as.prob=TRUE,baseOnly=TRUE)
methods=c("MI","FisherExact","CHISQ","FISHERZAMIR","CHISQZAMIR")
reps=100
Zamir=1

apply.tests= function (pairscount) 
{
  nr.methods = length(methods)
  stats = rep(NA, nr.methods)
  stats.index=1
  stats[stats.index]=mi.empirical(pairscount)
  stats.index=stats.index+1
  stats[stats.index]=fisher.test(pairscount,simulate.p.value=TRUE,B=reps)$p.value
  stats.index=stats.index+1
  stats[stats.index]=chisq.test(pairscount,simulate.p.value=TRUE,B=reps)$p.value
  stats.index=stats.index+1
  stats[stats.index]=fisher.test(pairscount+Zamir,simulate.p.value=TRUE,B=reps)$p.value
  stats.index=stats.index+1
  stats[stats.index]=chisq.test(pairscount+Zamir,simulate.p.value=TRUE,B=reps)$p.value
  return(stats)
}
# Remove lower case
badcolumns=NULL
AlignmentMatrix=as.matrix(Alignments)
for (j in 1:dim(AlignmentMatrix)[2]){
  if ( (AlignmentMatrix[1,j] != toupper(AlignmentMatrix[1,j]) ) |  (AlignmentMatrix[1,j] ==".") )
    badcolumns=c(badcolumns,j)
}
AlignmentMatrix=AlignmentMatrix[,-badcolumns]

d=dim(AlignmentMatrix)[2]
#d=2 for debug 
omat = matrix(NA, nrow = (d-1)*d/2, ncol = 2+length(methods))
omatnogaps = matrix(NA, nrow = (d-1)*d/2, ncol = 2+length(methods))
colnames(omat)=c("index1","index2",methods)
colnames(omatnogaps)=c("index1","index2",methods)
count=0
for (i in 1:(d-1)){
  print(i)
  for (j in (i+1):d){
    count=count+1
    pairscount=matrix(table(interaction(AlignmentMatrix[,i],AlignmentMatrix[,j])),nrow=length(as.vector(table(AlignmentMatrix[,i]))), ncol=length(as.vector(table(AlignmentMatrix[,j]))))
    omat[count,]=c(i,j,apply.tests(pairscount))
    # Also perform tests without gaps and without rows and columns that became empty after gaps were removed
    gaps=unlist(strsplit(names(table(interaction(AlignmentMatrix[,i],AlignmentMatrix[,j])))[1],""))
    if (gaps[3]=="-"){
      pairscount=pairscount[-1,]
      badcolumns=NULL
      for (k in 1:dim(pairscount)[2]){
        if (all(pairscount[,k]==0))
          badcolumns=c(badcolumns,k)
      }
      if (length(badcolumns)>0)
        pairscount=pairscount[,-badcolumns]
    }
    if (gaps[1]=="-"){
      pairscount=pairscount[,-1]
      badrows=NULL
      for (k in 1:dim(pairscount)[1]){
        if (all(pairscount[k,]==0))
          badrows=c(badrows,k)
      }
      if (length(badrows)>0)
        pairscount=pairscount[-badrows,]
    }
    omatnogaps[count,]=c(i,j,apply.tests(pairscount))
  }
}



