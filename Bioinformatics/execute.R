#####
ProfileMostProbable<-function(Text,k,profile)
{
  t<-nchar(Text)-k+1
  prob<-rep(1,t)
  for (i in 1:t)
  {
    temp<-substring(Text,i,i+k-1)
    for (j in 1:k)
    {
      ttemp<-substring(temp,j,j)
      if(ttemp=="A")
        prob[i]<-prob[i]*profile[1,j]
      if(ttemp=="C")
        prob[i]<-prob[i]*profile[2,j]
      if(ttemp=="G")
        prob[i]<-prob[i]*profile[3,j]
      if(ttemp=="T")
        prob[i]<-prob[i]*profile[4,j]
    }
  }
  tttemp<-which.max(prob)
  return(substring(Text,tttemp,tttemp+k-1))
}
#####
library(stringr)
ScoreConstruction<-function(Motifs)
{
  ScoreMatrix<-matrix(NA,nrow=4,ncol=nchar(Motifs[1]))
  temp<-strsplit(Motifs[1],"")[[1]]
  t<-length(Motifs)
  for (i in 2:t)
  {
    temp<-paste(temp,strsplit(Motifs[i],"")[[1]],sep="")
  }
  tt<-nchar(Motifs[1])
  for (j in 1:tt)
  {
    ttemp<-nchar(temp[j])
    ScoreMatrix[1,j]<-ttemp-str_count(temp[j],"A")
    ScoreMatrix[2,j]<-ttemp-str_count(temp[j],"C")
    ScoreMatrix[3,j]<-ttemp-str_count(temp[j],"G")
    ScoreMatrix[4,j]<-ttemp-str_count(temp[j],"T")
  }
  return(apply(ScoreMatrix,2,min))
}
#####
ProfileConstruction<-function(Motifs)
{
  profile<-matrix(NA,nrow=4,ncol=nchar(Motifs[1]))
  temp<-strsplit(Motifs[1],"")[[1]]
  t<-length(Motifs)
  for (i in 2:t)
  {
    temp<-paste(temp,strsplit(Motifs[i],"")[[1]],sep="")
  }
  tt<-nchar(Motifs[1])
  for (j in 1:tt)
  {
    ttemp<-nchar(temp[j])+1
    profile[1,j]<-(str_count(temp[j],"A")+1)/ttemp
    profile[2,j]<-(str_count(temp[j],"C")+1)/ttemp
    profile[3,j]<-(str_count(temp[j],"G")+1)/ttemp
    profile[4,j]<-(str_count(temp[j],"T")+1)/ttemp
  }
  return(profile)
}
#####
MotifsFinding<-function(Profile,Dna,k)
{
  t<-length(Dna)
  Motifs<-rep(NA,t)
  for (i in 1:t)
  {
    Motifs[i]<-ProfileMostProbable(Dna[i],k,Profile)
  }
  return(Motifs)
}
#####
RandomizedMotifSearch<-function(Dna,k)
{
  t<-length(Dna)
  Motifs<-rep(NA,t)
  for (i in 1:t)
  {
    temp<-sample(1:(nchar(Dna[i])-k+1),1)
    Motifs[i]<-substring(Dna[i],temp,temp+k-1)
  }
  BestMotifs<-Motifs
  while (TRUE)
  {
    Profile<-ProfileConstruction(Motifs)
    Motifs<-MotifsFinding(Profile,Dna,k)
    if (sum(ScoreConstruction(Motifs))<sum(ScoreConstruction(BestMotifs)))
      BestMotifs<-Motifs
    else return(BestMotifs)
  }
}
setwd("C:/Users/Nan Sun/Desktop")
Text<-paste(readLines("dataset_161_5.txt"),collapse="")
n<-20
k<-15
Dna<-rep(NA,n)
for (i in 1:n)
{
  Dna[i]<-substring(Text,1+(i-1)*nchar(Text)/n,i*nchar(Text)/n)
}
#####
f<-function(i){RandomizedMotifSearch(Dna,k)}