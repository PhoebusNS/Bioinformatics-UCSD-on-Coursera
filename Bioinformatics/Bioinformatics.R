#####PatternCount
PatternCount<-function(Text,Pattern) 
{
  count<-0;
  for (i in 1:(nchar(Text)-nchar(Pattern)+1)) 
    {
    if (substring(Text,i,i+nchar(Pattern)-1)==Pattern) 
      {
      count<-count+1
      }
     }
  return(count)
}
Text<-"CGGGAAA"
PatternCount(Text,"CG")
#####FrequentWords
FrequentWords<-function(Text,k)
{
  FrequentPatterns<-NULL
  Count<-rep(0,length=nchar(Text)-k+1)
  for (i in 1:(nchar(Text)-k+1))
  {
    Pattern<-substring(Text,i,i+k-1)
    Count[i]<-PatternCount(Text,Pattern)
    maxCount<-max(Count)
  }
  for (i in 1:(nchar(Text)-k+1))
  {
    if (Count[i]==maxCount)
      FrequentPatterns=cbind(FrequentPatterns,substring(Text,i,i+k-1))
  }
  FrequentPatterns<-unique(unlist(strsplit(FrequentPatterns, " ")))
  return(FrequentPatterns)
}
FrequentWords(Text,2)
#####ReverseComplement
temp<-chartr("ACGT","TGCA",Text)
result<-paste(rev(unlist(strsplit(temp,split=""))),collapse="")
#####PatternMatching
PatternMatching<-function(Text,Pattern)
{
  temp<-integer(length=0)
  for(i in 1:(nchar(Text)-nchar(Pattern)+1))
  {
    if(substring(Text,i,i+nchar(Pattern)-1)==Pattern)
      temp<-cbind(temp,i-1)
  }
  return(temp)
}
fileName<-'dataset_3_5 (5).txt'
Text<-readChar(fileName,file.info(fileName)$size)
Pattern<-"CTTGATCAT"
result<-PatternMatching(Text,Pattern)
result<-toString(result[1,])
result<-gsub(",","",result)
print(cat(result,sep='\n'))
#####PatternToNumber
PatternToNumber<-function(Pattern)
{
  symbols="A|T|C|G"
  if(grepl(symbols,Pattern)==FALSE)
    return (0)
  symbol<-substring(Pattern,nchar(Pattern),nchar(Pattern))
  Prefix<-substring(Pattern,1,nchar(Pattern)-1)
  return(4*PatternToNumber(Prefix)+as.numeric(chartr("ACGT","0123",symbol)))
}
PatternToNumber("CGAGTACGGTTCCTCTC")
#####NumberToPattern
NumberToPattern<-function(index,k)
{
  if (k==1)
    return(chartr("0123","ACGT",toString(index)))
  prefixIndex<-index%/%4
  r<-index%%4
  symbol<-chartr("0123","ACGT",toString(r))
  PrefixPattern<-NumberToPattern(prefixIndex,k-1)
  return (paste(PrefixPattern,symbol,sep=""))
}
NumberToPattern(8323,8)
#####ComputingFrequencies
ComputingFrequencies<-function(Text,k)
{
  FrequencyArray<-rep(0,length=4^k)
  for (i in 1:(nchar(Text)-k+1))
  {
  Pattern<-substring(Text,i,i+k-1)
  j<-PatternToNumber(Pattern)
  FrequencyArray[j+1]<-FrequencyArray[j+1]+1
  }
  return(FrequencyArray)
}
result<-ComputingFrequencies("CGTCTTCGCGTGATACGGTATCTTCCCAGCCGGCCATTAGCCCACAATTCCACGGTCCCCGGCTGTTCTGAAGACCCCTACGTAATGACACGGTGCGGTAGATATCGCGCACTGCCCCTAAAACCCCATTTCTAACACGTTCACAGACTGCATCCGTGTTACGAACTAGGAGATCCTGAGATCCCAGGGAAGTCACGCTTTTGAGCTGTGTCATGTACCCGATTATGCAGCTGTAGCGGTCAAGTCTTAAGCGATGGAAGGGACTGTGGGCTTGCGGTATCGGCCCGACCAGCTCGTCTCCAAGCGTTTCGATCTGTACCGTTTGTTATTTTACTCTAGTTTGGCTCATATAGATCCTTTCTTTCGCCACAGATTCGGATTACGAGGTCAGTGTGAGTGAAATCTATGACGTAGTTCTATATACAATGAATGCTCCCCTCAAGGAAGGAATCAGAAGTTACGCCAGGATCTGCATTTCACTCGTACTTCACAAGGTTTCGATGGCCCAGTAGAGCTCACGTTGTTACAAAATGGCTGATAACTTTGCACTCTCTCCGTTAAAGTCAATGGGTTCTTTTCTGATATCCCTCACTTACCTATCTTAGCGAAACTTTTGTACCCCGTAGGCGCAATGGGCCATGGTGGACGTGTGGTGGTCTTCGCGATGGTTCGACCTAATGCGCTGATCGAGACTATGGAAAAACTCACGTGAACTACAGGCCTGGCCAGAAGCATCTCGTGAGTGTCTTCGCGGG",5)
result<-toString(result)
result<-gsub(",","",result)
print(cat(result,sep='\n'))
#####ClumpFinding
ClumpFinding<-function(Genome,k,t,L)
{
  FrequentPatterns<-NULL
  Clump<-rep(0,length=(4^k))
  Text<-substring(Genome,1,L)
  FrequencyArray<-ComputingFrequencies(Text,k)
  for (i in 1:(4^k))
  {
    if (FrequencyArray[i]>=t)
      Clump[i]<-1
  }
  for (i in 2:(nchar(Genome)-L+1))
  {
    FirstPattern<-substring(Genome,i-1,i-1+k-1)
    index<-PatternToNumber(FirstPattern)
    FrequencyArray[index+1]<-FrequencyArray[index+1]-1
    LastPattern<-substring(Genome,i+L-k,i+L-1)
    index<-PatternToNumber(LastPattern)
    FrequencyArray[index+1]<-FrequencyArray[index+1]+1
    if(FrequencyArray[index+1]>=t)
    {
      Clump[index+1]<-1
    }
  }
  for (i in 1:(4^k))
  {
    if (Clump[i]==1)
      {
      Pattern<-NumberToPattern(i-1,k)
      FrequentPatterns<-paste(FrequentPatterns,Pattern,sep=" ")
      }
  }
  return(FrequentPatterns)
}
ClumpFinding(singleString,k=9,t=3,L=500)
setwd("C:/Users/Nan Sun/Desktop")
singleString<-paste(readLines("test.txt"),collapse="")
#####MinimumSkew
MinimumSkew<-function(Genome)
{
  skew<-c(NA,length=(nchar(Genome)+1))
  skew[1]<-0
  for (i in 1:nchar(Genome))
  {
    if (substring(Genome,i,i)=="G")
      skew[i+1]<-skew[i]+1
    else if (substring(Genome,i,i)=="C")
      skew[i+1]<-skew[i]-1
    else
      skew[i+1]<-skew[i]
  }
  return(which(skew==min(skew))-1)
}
Genome<-"GCATACACTTCCCAGTAGGTACTG"
MinimumSkew(Genome)
setwd("C:/Users/Nan Sun/Desktop")
Genome<-paste(readLines("dataset_7_6 (1).txt"),collapse="")
#####HammingDistance
HammingDistance<-function(Genome1,Genome2)
{
  count<-0
  for (i in 1:nchar(Genome1))
  {
    if(substring(Genome1,i,i)!=substring(Genome2,i,i))
      count<-count+1
    else count<-count
  }
  return(count)
}
Genome1<-"CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG"
Genome2<-"ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT"
HammingDistance(Genome1,Genome2)
#####ApproximatePatternMatching
ApproximatePatternMatching<-function(Text,Pattern,d)
{
  temp<-integer(length=0)
  for(i in 1:(nchar(Text)-nchar(Pattern)+1))
  {
    if(HammingDistance(substring(Text,i,i+nchar(Pattern)-1),Pattern)<=d)
      temp<-cbind(temp,i-1)
  }
  return(temp)
}
setwd("C:/Users/Nan Sun/Desktop")
fileName<-'dataset_9_4.txt'
Text<-readChar(fileName,file.info(fileName)$size)
Pattern<-"TTGTGTCCCGCG"
d<-5
result<-ApproximatePatternMatching(Text,Pattern,d)
result<-toString(result[1,])
result<-gsub(",","",result)
print(cat(result,sep='\n'))
#####ApproximatePatternCount
ApproximatePatternCount<-function(Text,Pattern,d) 
{
  count<-0;
  for (i in 1:(nchar(Text)-nchar(Pattern)+1)) 
  {
    if (HammingDistance(substring(Text,i,i+nchar(Pattern)-1),Pattern)<=d) 
    {
      count<-count+1
    }
  }
  return(count)
}
Text<-"CGTGACAGTGTATGGGCATCTTT"
ApproximatePatternCount(Text,"TGT",1)
#####Neighbors
Neighbors<-function(Pattern,d)
{
  if (d==0)
    return(Pattern)
  if (nchar(Pattern)==1)
    return(c("A","T","C","G"))
  else
    Neighborhood<-NULL
    SuffixNeighbors<-Neighbors(substring(Pattern,2,nchar(Pattern)),d)
    for(i in 1:length(SuffixNeighbors))
      if(HammingDistance(substring(Pattern,2,nchar(Pattern)),SuffixNeighbors[i])<d)
        Neighborhood=cbind(Neighborhood,
                           paste("A",SuffixNeighbors[i],sep=""),
                           paste("T",SuffixNeighbors[i],sep=""),
                           paste("C",SuffixNeighbors[i],sep=""),
                           paste("G",SuffixNeighbors[i],sep=""))
      else
        Neighborhood=cbind(Neighborhood,paste(substring(Pattern,1,1),SuffixNeighbors[i],sep=""))
    return(Neighborhood)
}
Pattern<-"ACGT"
temp<-Neighbors(Pattern,1)
temp<-toString(temp)
temp<-gsub(",","",temp)
print(cat(temp,sep=' '))
#####FrequentWordsWithMismatches
FrequentWordsWithMismatches<-function(Text,k,d)
{
  FrequentPatterns<-NULL
  Close<-rep(0,length=(4^k))
  FrequencyArray<-rep(0,length=(4^k))
  for(i in 1:(nchar(Text)-k+1))
  {
    Neighborhood<-Neighbors(substring(Text,i,i+k-1),d)
    for(j in 1:length(Neighborhood))
    {
      index<-PatternToNumber(Neighborhood[j])
      Close[index]<-1
    }
  }
  for (i in 1:(4^k))
  {
    if(Close[i]==1)
    {
      Pattern<-NumberToPattern(i,k)
      FrequencyArray[i]<-ApproximatePatternCount(Text,Pattern,d)
    } 
  }
  maxCount<-max(FrequencyArray)
  for(i in 1:(4^k))
  {
    if(FrequencyArray[i]==maxCount)
    {
      Pattern<-NumberToPattern(i,k)
      FrequentPatterns<-paste(FrequentPatterns,Pattern,sep=" ")
    } 
  }
  return(FrequentPatterns)
}
Text<-"TGAGAAACTGAGGGCGGCTGTGGCGGCAAACTGAGTGTTGAGTGAGTGAGAAACGCCGGCGCCTGTGCCGGCGCCGCCAAACGCCGCCGGCGGCGGCGGCGCCGGCTGTGCCTGAGTGAGTGAGAAACTGAGTGTTGAGTGAGGGCGCCGCCTGAGAAACGGCAAACAAACGGCTGTTGTGCCAAACTGTTGTTGTGGCAAACTGAGGGCGCCAAACTGAGAAACTGAGTGTAAACAAACTGAGGGCAAACTGTGCCGCCGCCGCCGCCTGTGGCGCCTGTTGAGGGCAAACTGTTGTAAACGGCAAACTGAGGCCTGTTGTGCCTGTTGAGGCCGGCGGCTGTTGTTGTGCCAAACGGCTGTTGAGAAACAAAC"
FrequentWordsWithMismatches(Text,5,2)
#####FrequentWordsWithMismatchesAndReverseComplements
FrequentWordsWithMismatchesAndReverseComplements<-function(Text,k,d)
{
  FrequentPatterns<-NULL
  Close<-rep(0,length=(4^k))
  FrequencyArray<-rep(0,length=(4^k))
  for(i in 1:(nchar(Text)-k+1))
  {
    Neighborhood<-Neighbors(substring(Text,i,i+k-1),d)
    for(j in 1:length(Neighborhood))
    {
      index<-PatternToNumber(Neighborhood[j])
      Close[index]<-1
    }
  }
  for (i in 1:(4^k))
  {
    if(Close[i]==1)
    {
      Pattern<-NumberToPattern(i,k)
      Pattern_rc<-paste(rev(unlist(strsplit(chartr("ACGT","TGCA",Pattern),split=""))),collapse="")
      FrequencyArray[i]<-ApproximatePatternCount(Text,Pattern,d)+ApproximatePatternCount(Text,Pattern_rc,d)
    } 
  }
  maxCount<-max(FrequencyArray)
  for(i in 1:(4^k))
  {
    if(FrequencyArray[i]==maxCount)
    {
      Pattern<-NumberToPattern(i,k)
      FrequentPatterns<-paste(FrequentPatterns,Pattern,sep=" ")
    } 
  }
  return(FrequentPatterns)
}
setwd("C:/Users/Nan Sun/Desktop")
Text<-paste(readLines("dataset_9_8 (1).txt"),collapse="")
MinimumSkew(Text)
FrequentWordsWithMismatchesAndReverseComplements(substring(Text,3764856-500,3764856+500),9,1)
#####MotifEnumeration
MotifEnumeration<-function(Dna,k,d)
{
  Patterns<-NULL
  for(i in 1:(nchar(Dna[1])-k+1))
  {
    temp<-Neighbors(substring(Dna[1],i,k+i-1),d)
    for(j in 1:length(temp))
    {
      count<-0
      for(l in 2:length(Dna))
      {
        count<-count+ifelse(ApproximatePatternCount(Dna[l],temp[j],d)>0,1,0)
      }
      if(count>=(length(Dna)-1))
        Patterns=cbind(Patterns,temp[j])
    }
  }
  Patterns<-unique(unlist(strsplit(Patterns, " ")))
  return(Patterns)
}
Dna<-c("CGCGGACAACTCCACCTATGAGGGC","TCGTCACCTCTGAAGAGCGCGTCGC","TCGCCAGGGCCCAGGTCAGTCAACT","TTGGGAACCGGGAGCTGACCGATAT","AGCGCTGAGCCGTTAACGTTGGAAC","AGGGCTAAGGTACTATCGACCTGCA")
temp<-MotifEnumeration(Dna,5,2)
temp<-toString(temp)
temp<-gsub(",","",temp)
print(cat(temp,sep=' '))
#####Entropy
8*0.1*log2(0.1)+4*0.2*log2(0.2)+2*0.3*log2(0.3)+3*0.4*log2(0.4)+1*0.5*log2(0.5)+2*0.6*log2(0.6)+2*0.7*log2(0.7)+1*0.8*log2(0.8)+3*0.9*log2(0.9)
#####DistanceBetweenPatternAndStrings
DistanceBetweenPatternAndStrings<-function(Pattern,Dna)
{
  k<-nchar(Pattern)
  distance<-0
  for (i in 1:length(Dna))
  {
    HD<-Inf
    for (j in 1:(nchar(Dna[i])-k+1))
    {
      if(HD>HammingDistance(Pattern,substring(Dna[i],j,j+k-1)))
        HD<-HammingDistance(Pattern,substring(Dna[i],j,j+k-1))
    }
    distance<-distance+HD
  }
  return(distance)
}
#####MedianString
MedianString<-function(Dna,k)
{
  distance<-Inf
  for (i in 1:(4^k))
  {
    Pattern<-NumberToPattern(i-1,k)
    if (distance>DistanceBetweenPatternAndStrings(Pattern,Dna))
    {
      distance<-DistanceBetweenPatternAndStrings(Pattern,Dna)
      Median<-Pattern
    }
  }
  return(Median)
}
Dna<-c("CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCCGCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTCGGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG")
MedianString(Dna,7)
#####ProfileMostProbable
Text<-c("ATGCGCTGCACGGATCTTCGCGGGATAGTAGCAACAGAACAGAATGGGATAGCAATGAGTGATAGTGGTTTACATATGGAGTACCATTACTCTCAGAACGGGTTGTAGTGAAGGTCCTATCAGATTCGCAACTAAAAGGATTGAGTATCGTGAAAGTCCGACCCCATTGTAACCTCCGACCTAATGTGTGCACTGATCAACCTACCCGACTCCATGTTCGAGTGGTGCTACAACCTAGGTGCGGATCCGATTAAAGACCCCTGGTAATTGAATCTCACTTACGGGAGGTGTATCATCGCCTTCCCTCAACAAGTACACATTGCAGTATAGCAGGAGCTTTCAGCGGGTTTGACGGTAGGGAAATTACTGAAACACTCGACCCACTCCGTACGGATCCTTTTACCAGGCTTGTCGAAGTCCCCCGGTTAGAGATACAGTGTGTGGTAGTCATCTCGACCAATTCGGGTTGTGTGAGATTGACACAAACGAATCCTAGTTTAATAGAATACGAATAACTTTGAAGATTTTATGCTGAGAAGGGAGTTATGTAAATGGTCCACGTCCGTAAGGTAAGGCTCAAGCTTCGGTAGGGTTTAGTGCGAAGACCGATTCAGCGCTCGGCAGTTCGGGGCCCTAAAATAGATTACGGGACTTGATGAGTAGCCGGGGGCGCACGCACACTGCGATAAAGGTTCCCTGGTAGATGTGACCTTCTGTAAAAATCCGGTGGGTTCCGTAGGTCATATGAGTTCTGATTCGAGGGGGACCTGTGTTGCGCATGATGAGGACAAAATAGCAGGATAGAGGTTGCAGATACATCAGTTGCCCGCGCTTGCATATTCTCCGCAGCACGAATCATATCACCTTTTATTATGTAGTGACTCTTGAGGCACTTGTCTATACTGGAAGAACCTAGCGCGAGGGCCCATAGCTAGGCTCAACAACAGCCGGGACTCTAGACTCTTGCGAGGAACAATTGGTAGATCATCCGGAGCTTCCA")
k<-13
profile<-read.table("C:/Users/Nan Sun/Desktop/profile.txt")
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
ProfileMostProbable(Text,k,profile)
#####ProfileConstruction
library(stringr)
ProfileConstruction<-function(Motifs)
{
  profile<-matrix(NA,nrow=4,ncol=nchar(Motifs[1]))
  temp<-strsplit(Motifs[1],"")[[1]]
  for (i in 2:length(Motifs))
  {
    temp<-paste(temp,strsplit(Motifs[i],"")[[1]],sep="")
  }
  for (j in 1:nchar(Motifs[1]))
  {
   profile[1,j]<-str_count(temp[j],"A")/nchar(temp[j])
   profile[2,j]<-str_count(temp[j],"C")/nchar(temp[j])
   profile[3,j]<-str_count(temp[j],"G")/nchar(temp[j])
   profile[4,j]<-str_count(temp[j],"T")/nchar(temp[j])
  }
  return(profile)
}
Motifs<-c("TGA","GTT","GAA","TGT")
Dna<-c("TGACGTTC","TAAGAGTT","GGACGAAA","CTGTTCGC")
profile<-ProfileConstruction(Motifs)
#####Score
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
ScoreConstruction(Motifs)
#####GreedyMotifSearch
Text<-paste(readLines("dataset_159_5 (3).txt"),collapse="")
n<-25
k<-12
Dna<-rep(NA,n)
for (i in 1:n)
{
  Dna[i]<-substring(Text,1+(i-1)*nchar(Text)/n,i*nchar(Text)/n)
}
GreedyMotifSearch<-function(Dna,k)
{
  BestMotifs<-NULL
  for (i in 1:length(Dna))
  {
    BestMotifs<-cbind(BestMotifs,substring(Dna[i],1,1+k-1))
  }
  for (i in 1:(nchar(Dna[1])-k+1))
  {
    Motifs<-substring(Dna[1],i,i+k-1)
    temp<-NULL
    for (j in 2:length(Dna))
    {
      profile<-ProfileConstruction(cbind(temp,Motifs[j-1]))
      temp<-cbind(temp,Motifs[j-1])
      Motifs<-cbind(Motifs,ProfileMostProbable(Dna[j],k,profile))
    }
    if(sum(ScoreConstruction(Motifs))<=sum(ScoreConstruction(BestMotifs)))
      BestMotifs<-Motifs
  }
  return(BestMotifs)
}
temp<-GreedyMotifSearch(Dna,k)
temp<-toString(temp)
temp<-gsub(",","",temp)
print(cat(temp,sep=' '))
#####Using Laplace's Rule
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
profile<-ProfileConstruction(Motifs)
#####MotifsFinding
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
k<-3
MotifsFinding(profile,Dna,k)
#####RandomizedMotifSearch
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
Dna<-c("AATTGGCACATCATTATCGATAACGATTCGCCGCATTGCC","GGTTAACATCGAATAACTGACACCTGCTCTGGCACCGCTC","AATTGGCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT","GGTTAATGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG","AATTGGACGGCAACTACGGTTACAACGCAGCAAGAATATT","GGTTAACTGTTGTTGCTAACACCGTTAAGCGACGGCAACT","AATTGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG","GGTTAAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA")
k<-6
f<-function(i)
{
  RandomizedMotifSearch(Dna,k)
}
result<-lapply(1:1000,f)
temp<-rep(NA,1000)
for (i in 1:1000)
{
  temp[i]<-sum(ScoreConstruction(result[[i]]))
}
result[[which.min(temp)]]
#####
Text<-paste(readLines("dataset_161_5.txt"),collapse="")
n<-20
k<-15
Dna<-rep(NA,n)
for (i in 1:n)
{
  Dna[i]<-substring(Text,1+(i-1)*nchar(Text)/n,i*nchar(Text)/n)
}
f<-function(i)
{
  RandomizedMotifSearch(Dna,k)
}
result<-lapply(1:1000,f)
temp<-rep(NA,1000)
for (i in 1:1000)
{
  temp[i]<-sum(ScoreConstruction(result[[i]]))
}
temp<-result[[which.min(temp)]]
temp<-toString(temp)
temp<-gsub(",","",temp)
print(cat(temp,sep='\n'))
#####
p<-(600-15)/(600-15+1)
1-p^{10}
((1-p)^{2})*(p^{8})*45
#####GibbsSampler
Text<-c("ATGCGCTGCACGGATCTTCGCGGGATAGTAGCAACAGAACAGAATGGGATAGCAATGAGTGATAGTGGTTTACATATGGAGTACCATTACTCTCAGAACGGGTTGTAGTGAAGGTCCTATCAGATTCGCAACTAAAAGGATTGAGTATCGTGAAAGTCCGACCCCATTGTAACCTCCGACCTAATGTGTGCACTGATCAACCTACCCGACTCCATGTTCGAGTGGTGCTACAACCTAGGTGCGGATCCGATTAAAGACCCCTGGTAATTGAATCTCACTTACGGGAGGTGTATCATCGCCTTCCCTCAACAAGTACACATTGCAGTATAGCAGGAGCTTTCAGCGGGTTTGACGGTAGGGAAATTACTGAAACACTCGACCCACTCCGTACGGATCCTTTTACCAGGCTTGTCGAAGTCCCCCGGTTAGAGATACAGTGTGTGGTAGTCATCTCGACCAATTCGGGTTGTGTGAGATTGACACAAACGAATCCTAGTTTAATAGAATACGAATAACTTTGAAGATTTTATGCTGAGAAGGGAGTTATGTAAATGGTCCACGTCCGTAAGGTAAGGCTCAAGCTTCGGTAGGGTTTAGTGCGAAGACCGATTCAGCGCTCGGCAGTTCGGGGCCCTAAAATAGATTACGGGACTTGATGAGTAGCCGGGGGCGCACGCACACTGCGATAAAGGTTCCCTGGTAGATGTGACCTTCTGTAAAAATCCGGTGGGTTCCGTAGGTCATATGAGTTCTGATTCGAGGGGGACCTGTGTTGCGCATGATGAGGACAAAATAGCAGGATAGAGGTTGCAGATACATCAGTTGCCCGCGCTTGCATATTCTCCGCAGCACGAATCATATCACCTTTTATTATGTAGTGACTCTTGAGGCACTTGTCTATACTGGAAGAACCTAGCGCGAGGGCCCATAGCTAGGCTCAACAACAGCCGGGACTCTAGACTCTTGCGAGGAACAATTGGTAGATCATCCGGAGCTTCCA")
k<-13
profile<-read.table("C:/Users/Nan Sun/Desktop/profile.txt")
ProbConstruction<-function(Text,k,profile)
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
  return(prob)
}
ProbConstruction(Text,k,profile)
#####
GibbsSampler<-function(Dna,k,N)
{
  t<-length(Dna)
  Motifs<-rep(NA,t)
  for (i in 1:t)
  {
    temp<-sample(1:(nchar(Dna[i])-k+1),1)
    Motifs[i]<-substring(Dna[i],temp,temp+k-1)
  }
  BestMotifs<-Motifs
  for (j in 1:N)
  {
    i<-sample(1:t,1)
    profile<-ProfileConstruction(Motifs[-i])
    myprob<-ProbConstruction(Dna[i],k,profile)
    myloc<-sample(1:(nchar(Dna[i])-k+1),1,prob=(myprob/sum(myprob)))
    Motifs[i]<-substring(Dna[i],myloc,myloc+k-1)
    if (sum(ScoreConstruction(Motifs))<sum(ScoreConstruction(BestMotifs)))
      BestMotifs<-Motifs
  }
  return(BestMotifs)
}
N<-100
k<-8
Text<-paste(readLines("dataset_161_5.txt"),collapse="")
n<-5
Dna<-rep(NA,n)
for (i in 1:n)
{
  Dna[i]<-substring(Text,1+(i-1)*nchar(Text)/n,i*nchar(Text)/n)
}
ff<-function(i)
{
  GibbsSampler(Dna,k,N)
}
result<-lapply(1:20,ff)
temp<-rep(NA,20)
for (i in 1:20)
{
  temp[i]<-sum(ScoreConstruction(result[[i]]))
}
temp<-result[[which.min(temp)]]
print(cat(temp,sep='\n'))
###################
StringComposition<-function(k,Text)
{
  Composition<-rep(NA,nchar(Text)-k+1)
  for (i in 1:(nchar(Text)-k+1))
  {
    Composition[i]<-substring(Text,i,i+k-1)
  }
  return(Composition)  
}
setwd("C:/Users/Nan Sun/Desktop")
Text<-paste(readLines("dataset_197_3.txt"),collapse="")
k<-100
result<-StringComposition(k,Text)
cat(result,file="outfile.txt",sep='\n')
#####
n<-length(readLines("dataset_198_3.txt"))
PatternOriginal<-paste(readLines("dataset_198_3.txt"),collapse="")
Pattern<-rep(NA,n)
Pattern[1]<-c("ATTCATACGTACATGCACTACCTAG")
k<-nchar(Pattern[1])
for (i in 2:n)
{
  Pattern[i]<-substring(PatternOriginal,(i-1)*k+1,i*k)
}
ReconstructFromGenomePath<-function(Pattern,k,n)
{
  Text<-Pattern[1]
  for (i in 2:n)
  {
    Text<-paste(Text,substring(Pattern[i],k,k),sep="")
  }
  return(Text)
}
ReconstructFromGenomePath(Pattern,k,n)
#####
n<-length(readLines("dataset_198_3.txt"))
PatternOriginal<-paste(readLines("dataset_198_3.txt"),collapse="")
Pattern<-rep(NA,n)
Pattern[1]<-c("ATACTACTAGCAAGAGGTCCGTGGC")
k<-nchar(Pattern[1])
for (i in 2:n)
{
  Pattern[i]<-substring(PatternOriginal,(i-1)*k+1,i*k)
}
OverlapGraph<-function(Pattern,k,n)
{
  Pre_Suf<-matrix(NA,nrow=n,ncol=2)
  result<-matrix(NA,ncol=2)
  for (i in 1:n)
  {
    Pre_Suf[i,1]<-substring(Pattern[i],1,k-1)
    Pre_Suf[i,2]<-substring(Pattern[i],2,k)
  }
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      if (Pre_Suf[i,2]==Pre_Suf[j,1] & Pattern[i]!=Pattern[j])
      {
        temp<-cbind(Pattern[i],Pattern[j])
        result<-rbind(result,temp)
      }
    }
  }
  return(result[-1,])
}
result<-OverlapGraph(Pattern,k,n)
temp<-matrix(NA,nrow=nrow(result),ncol=1)
for (i in 1:nrow(result))
{
  temp[i]<-paste(result[i,],collapse=' -> ')
}
cat(temp,file="outfile.txt",sep='\n')
#####
Text<-paste(readLines("dataset_203_1.txt"),collapse="")
k<-3
Nodes<-StringComposition(k-1,Text)
n<-nchar(Text)
DeBruijnGraph<-function(Nodes,k,n)
{
  result<-matrix(NA,ncol=2)
  for (i in 1:(n-k+1))
  {
    temp<-cbind(Nodes[i],Nodes[i+1])
    result<-rbind(result,temp)
  }
  return(result[-1,])
}
result<-DeBruijnGraph(Nodes,k,n)
result<-data.frame(result)
result<-aggregate(X2~X1,data=result,toString)
temp<-matrix(NA,nrow=nrow(result),ncol=1)
for (i in 1:nrow(result))
{
  temp[i]<-paste(result[i,]$X1,result[i,]$X2,sep=' -> ')
}
temp<-gsub(", ",",",temp)
cat(temp,file="outfile.txt",sep='\n')
#####
n<-length(readLines("dataset_198_3.txt"))
PatternOriginal<-paste(readLines("dataset_198_3.txt"),collapse="")
Pattern<-rep(NA,n)
Pattern[1]<-c("AAAT")
k<-nchar(Pattern[1])
for (i in 2:n)
{
  Pattern[i]<-substring(PatternOriginal,(i-1)*k+1,i*k)
}
DeBruijnGraphFromKmer<-function(Pattern,k,n)
{
  result<-matrix(NA,ncol=2)
  for (i in 1:n)
  {
    temp<-cbind(substring(Pattern[i],1,k-1),substring(Pattern[i],2,k))
    result<-rbind(result,temp)
  }
  return(result[-1,])
}
result<-DeBruijnGraphFromKmer(Pattern,k,n)
result<-data.frame(result)
result<-aggregate(X2~X1,data=result,toString)
temp<-matrix(NA,nrow=nrow(result),ncol=1)
for (i in 1:nrow(result))
{
  temp[i]<-paste(result[i,]$X1,result[i,]$X2,sep=' -> ')
}
temp<-gsub(", ",",",temp)
cat(temp,file="outfile.txt",sep='\n')
#####
library(stringr)
setwd("C:/Users/Nan Sun/Desktop")
AdjacencyList<-readLines("outfile.txt")
for (i in 1:length(AdjacencyList))
{
  AdjacencyList[i]<-sub(' -> ',',',AdjacencyList[i])
}
Nodes<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Nodes[i]<-strsplit(AdjacencyList[i],",")[[1]][1]
}
Neighbors<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Neighbors[i]<-list(strsplit(AdjacencyList[i],",")[[1]][-1])
}
Outdegree<-rep(NA,length(AdjacencyList))
for (i in 1:length(Nodes))
{
  Outdegree[i]<-length(Neighbors[[i]])
}
Indegree<-rep(NA,length(AdjacencyList))
useful<-unlist(Neighbors,",")
for (i in 1:length(Nodes))
{
  Indegree[i]<-length(which(useful==Nodes[i]))
}
#####
Circuit<-NULL
Stack<-NULL
CurrentNode<-sample(Nodes,1)
while (length(Neighbors[which(Nodes==CurrentNode)][[1]])>0 || length(Stack)!=0)
{
  if(length(Neighbors[which(Nodes==CurrentNode)][[1]])==0)
  {
    Circuit<-append(Circuit,CurrentNode)
    CurrentNode<-Stack[length(Stack)]
    Stack<-Stack[-length(Stack)]
  }
  if(length(Neighbors[which(Nodes==CurrentNode)][[1]])>0)
  {
    Stack<-append(Stack,CurrentNode)
    CurrentNodeNeighbors<-unlist(Neighbors[which(Nodes==CurrentNode)][[1]])
    toy<-sample(CurrentNodeNeighbors,1)
    Neighbors[which(Nodes==CurrentNode)][[1]]<-CurrentNodeNeighbors[CurrentNodeNeighbors!=toy]
    CurrentNode<-toy
  }
}
temp<-rev(Circuit)
cat(temp,file="outfile.txt",sep='->')
##########
#####
library(stringr)
setwd("C:/Users/Nan Sun/Desktop")
AdjacencyList<-readLines("outfile.txt")
for (i in 1:length(AdjacencyList))
{
  AdjacencyList[i]<-sub(' -> ',',',AdjacencyList[i])
}
Nodes<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Nodes[i]<-strsplit(AdjacencyList[i],",")[[1]][1]
}
Neighbors<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Neighbors[i]<-list(strsplit(AdjacencyList[i],",")[[1]][-1])
}
Outdegree<-rep(NA,length(AdjacencyList))
for (i in 1:length(Nodes))
{
  Outdegree[i]<-length(Neighbors[[i]])
}
Indegree<-rep(NA,length(AdjacencyList))
useful<-unlist(Neighbors,",")
for (i in 1:length(Nodes))
{
  Indegree[i]<-length(which(useful==Nodes[i]))
}
#####
Circuit<-NULL
Stack<-NULL
CurrentNode<-Nodes[which(Outdegree-Indegree==1)]
x<-Nodes
y<-unique(unlist(Neighbors))
Nodes[length(Nodes)+1]<-y[!y%in%x]
Neighbors[length(Neighbors)+1]<-NULL
while (length(Neighbors[which(Nodes==CurrentNode)][[1]])>0 || length(Stack)!=0)
{
  if(length(Neighbors[which(Nodes==CurrentNode)][[1]])==0)
  {
    Circuit<-append(Circuit,CurrentNode)
    CurrentNode<-Stack[length(Stack)]
    Stack<-Stack[-length(Stack)]
  }
  if(length(Neighbors[which(Nodes==CurrentNode)][[1]])>0)
  {
    Stack<-append(Stack,CurrentNode)
    CurrentNodeNeighbors<-unlist(Neighbors[which(Nodes==CurrentNode)][[1]])
    toy<-sample(CurrentNodeNeighbors,1)
    Neighbors[which(Nodes==CurrentNode)][[1]]<-CurrentNodeNeighbors[CurrentNodeNeighbors!=toy]
    CurrentNode<-toy
  }
}
temp<-rev(Circuit)
result<-temp[1]
for (i in 2:length(temp))
{
  result<-cat(result,substring(temp[i],nchar(temp[i]),nchar(temp[i])),sep="")
}
cat(temp,file="outfile.txt",sep='->')
#####
n<-length(readLines("dataset_203_1.txt"))
PatternOriginal<-paste(readLines("dataset_203_1.txt"),collapse="")
Pattern<-rep(NA,n)
Pattern[1]<-c("AGAAGACGTCCCCCACAGGTAAACG")
k<-nchar(Pattern[1])
for (i in 2:n)
{
  Pattern[i]<-substring(PatternOriginal,(i-1)*k+1,i*k)
}
result<-DeBruijnGraphFromKmer(Pattern,k,n)
result<-data.frame(result)
result<-aggregate(X2~X1,data=result,toString)
temp<-matrix(NA,nrow=nrow(result),ncol=1)
for (i in 1:nrow(result))
{
  temp[i]<-paste(result[i,]$X1,result[i,]$X2,sep=' -> ')
}
temp<-gsub(", ",",",temp)
cat(temp,file="outfile.txt",sep='\n')
AdjacencyList<-readLines("outfile.txt")
for (i in 1:length(AdjacencyList))
{
  AdjacencyList[i]<-sub(' -> ',',',AdjacencyList[i])
}
Nodes<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Nodes[i]<-strsplit(AdjacencyList[i],",")[[1]][1]
}
Neighbors<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Neighbors[i]<-list(strsplit(AdjacencyList[i],",")[[1]][-1])
}
Outdegree<-rep(NA,length(AdjacencyList))
for (i in 1:length(Nodes))
{
  Outdegree[i]<-length(Neighbors[[i]])
}
Indegree<-rep(NA,length(AdjacencyList))
useful<-unlist(Neighbors,",")
for (i in 1:length(Nodes))
{
  Indegree[i]<-length(which(useful==Nodes[i]))
}
#####
Circuit<-NULL
Stack<-NULL
CurrentNode<-Nodes[which(Outdegree-Indegree==1)]
x<-Nodes
y<-unique(unlist(Neighbors))
Nodes[length(Nodes)+1]<-y[!y%in%x]
Neighbors[length(Neighbors)+1]<-NULL
while (length(Neighbors[which(Nodes==CurrentNode)][[1]])>0 || length(Stack)!=0)
{
  if(length(Neighbors[which(Nodes==CurrentNode)][[1]])==0)
  {
    Circuit<-append(Circuit,CurrentNode)
    CurrentNode<-Stack[length(Stack)]
    Stack<-Stack[-length(Stack)]
  }
  if(length(Neighbors[which(Nodes==CurrentNode)][[1]])>0)
  {
    Stack<-append(Stack,CurrentNode)
    CurrentNodeNeighbors<-unlist(Neighbors[which(Nodes==CurrentNode)][[1]])
    toy<-sample(CurrentNodeNeighbors,1)
    Neighbors[which(Nodes==CurrentNode)][[1]]<-CurrentNodeNeighbors[CurrentNodeNeighbors!=toy]
    CurrentNode<-toy
  }
}
temp<-rev(Circuit)
result<-temp[1]
for (i in 2:length(temp))
{
  result<-cat(result,substring(temp[i],nchar(temp[i]),nchar(temp[i])),sep="")
}
###############
AdjacencyList<-t(expand.grid(replicate(9,0:1,simplify=FALSE)))
AdjacencyList<-as.matrix(AdjacencyList)
cat(AdjacencyList,file="outfile.txt",sep="")
n<-ncol(AdjacencyList)
PatternOriginal<-readLines("outfile.txt")
Pattern<-rep(NA,n)
Pattern[1]<-c("000000000")
k<-nchar(Pattern[1])
for (i in 2:n)
{
  Pattern[i]<-substring(PatternOriginal,(i-1)*k+1,i*k)
}
result<-DeBruijnGraphFromKmer(Pattern,k,n)
result<-data.frame(result)
result<-aggregate(X2~X1,data=result,toString)
temp<-matrix(NA,nrow=nrow(result),ncol=1)
for (i in 1:nrow(result))
{
  temp[i]<-paste(result[i,]$X1,result[i,]$X2,sep=' -> ')
}
temp<-gsub(", ",",",temp)
cat(temp,file="outfile1.txt",sep='\n')
AdjacencyList<-readLines("outfile1.txt")
for (i in 1:length(AdjacencyList))
{
  AdjacencyList[i]<-sub(' -> ',',',AdjacencyList[i])
}
Nodes<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Nodes[i]<-strsplit(AdjacencyList[i],",")[[1]][1]
}
Neighbors<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Neighbors[i]<-list(strsplit(AdjacencyList[i],",")[[1]][-1])
}
Circuit<-NULL
Stack<-NULL
CurrentNode<-sample(Nodes,1)
while (length(Neighbors[which(Nodes==CurrentNode)][[1]])>0 || length(Stack)!=0)
{
  if(length(Neighbors[which(Nodes==CurrentNode)][[1]])==0)
  {
    Circuit<-append(Circuit,CurrentNode)
    CurrentNode<-Stack[length(Stack)]
    Stack<-Stack[-length(Stack)]
  }
  if(length(Neighbors[which(Nodes==CurrentNode)][[1]])>0)
  {
    Stack<-append(Stack,CurrentNode)
    CurrentNodeNeighbors<-unlist(Neighbors[which(Nodes==CurrentNode)][[1]])
    toy<-sample(CurrentNodeNeighbors,1)
    Neighbors[which(Nodes==CurrentNode)][[1]]<-CurrentNodeNeighbors[CurrentNodeNeighbors!=toy]
    CurrentNode<-toy
  }
}
temp<-rev(Circuit)
result<-temp[1]
for (i in 2:length(temp))
{
  result<-cat(result,substring(temp[i],nchar(temp[i]),nchar(temp[i])),sep="")
}
##########
setwd("C:/Users/Nan Sun/Desktop")
PairedReads<-readLines("dataset_203_1.txt")
length(PairedReads)
Prefix<-function(PairedReads)
{
  temp<-rep(NA,length(PairedReads))
  for (i in 1:length(PairedReads))
  {
    temp[i]<-paste(substring(PairedReads[i],1,(nchar(PairedReads[1])-1)/2-1),substring(PairedReads[i],(nchar(PairedReads[1])-1)/2+2,nchar(PairedReads[1])-1),sep='|')
  }
  return(temp)
}
Prefix<-Prefix(PairedReads)
Suffix<-function(PairedReads)
{
  temp<-rep(NA,length(PairedReads))
  for (i in 1:length(PairedReads))
  {
    temp[i]<-paste(substring(PairedReads[i],2,(nchar(PairedReads[1])-1)/2),substring(PairedReads[i],(nchar(PairedReads[1])-1)/2+3,nchar(PairedReads[1])),sep='|')
  }
  return(temp)
}
Suffix<-Suffix(PairedReads)
OrderedSequence<-function(Prefix,Suffix,PairedReads)
{
  temp<-rep(NA,length(PairedReads))
  temp[1]<-PairedReads[!Prefix%in%Suffix]
  for (i in 2:length(PairedReads))
  {
   temp[i]<-PairedReads[which(Prefix==Suffix[which(PairedReads==temp[i-1])])]
  }
  return(temp)
}
GappedPatterns<-OrderedSequence(A,B,PairedReads)
#####
result1<-substring(GappedPatterns[1],1,(nchar(PairedReads[1])-1)/2)
for (i in 2:length(GappedPatterns))
{
  result1<-paste(result1,substring(GappedPatterns[i],(nchar(GappedPatterns[i])-1)/2,(nchar(GappedPatterns[i])-1)/2),sep="")
}
result2<-substring(GappedPatterns[1],(nchar(PairedReads[1])-1)/2+2,nchar(PairedReads[1]))
for (i in 2:length(GappedPatterns))
{
  result2<-paste(result2,substring(GappedPatterns[i],nchar(GappedPatterns[i]),nchar(GappedPatterns[i])),sep="")
}
k<-50
d<-200
temp1<-substring(result1,k+d+1,nchar(result1))
temp2<-substring(result2,1,nchar(result2)-k-d)
result<-paste(result1,substring(result2,nchar(temp1)+1,nchar(result2)),sep="")
result
###############
setwd("C:/Users/Nan Sun/Desktop")
Patterns<-readLines("dataset_203_1.txt")
k<-nchar(Patterns[1])
n<-length(Patterns)
result<-DeBruijnGraphFromKmer(Patterns,k,n)
result<-data.frame(result)
result<-aggregate(X2~X1,data=result,toString)
temp<-matrix(NA,nrow=nrow(result),ncol=1)
for (i in 1:nrow(result))
{
  temp[i]<-paste(result[i,]$X1,result[i,]$X2,sep=' -> ')
}
temp<-gsub(", ",",",temp)
cat(temp,file="outfile1.txt",sep='\n')
AdjacencyList<-readLines("outfile1.txt")
for (i in 1:length(AdjacencyList))
{
  AdjacencyList[i]<-sub(' -> ',',',AdjacencyList[i])
}
Nodes<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Nodes[i]<-strsplit(AdjacencyList[i],",")[[1]][1]
}
Neighbors<-rep(NA,length(AdjacencyList))
for (i in 1:length(AdjacencyList))
{
  Neighbors[i]<-list(strsplit(AdjacencyList[i],",")[[1]][-1])
}
Outdegree<-rep(NA,length(AdjacencyList))
for (i in 1:length(Nodes))
{
  Outdegree[i]<-length(Neighbors[[i]])
}
Indegree<-rep(NA,length(AdjacencyList))
useful<-unlist(Neighbors,",")
for (i in 1:length(Nodes))
{
  Indegree[i]<-length(which(useful==Nodes[i]))
}
Location<-which(!(Indegree==1 & Outdegree==1) & Outdegree>0)
IsoLocation<-which(Indegree==1 & Outdegree==1)
Paths<-NULL
for (i in 1:length(Location))
{
  for (j in 1:length(Neighbors[Location[i]][[1]]))
  {
    NonBranchingPath<-paste(Nodes[Location[i]],substring(Neighbors[Location[i]][[1]][j],nchar(Neighbors[Location[i]][[1]][j]),nchar(Neighbors[Location[i]][[1]][j])),sep="")
    w<-unlist(Neighbors[Location[i]][[1]][j])
    while (length(intersect(which(Nodes==w),IsoLocation))!=0)
    {
      u<-unlist(Neighbors[Nodes==w])
      NonBranchingPath<-paste(NonBranchingPath,substring(u,nchar(u),nchar(u)),sep="")
      w<-u
    }
    Paths<-append(Paths,NonBranchingPath)
  }
}
cat(Paths,file="outfile2.txt",sep='\n')
##################
setwd("C:/Users/Nan Sun/Desktop")
GeneticCode<-readLines("RNA_codon_table_1.txt")
Pattern<-readLines("paste.txt")
number<-nchar(Pattern)/3
temp<-NULL
for (i in 1:number)
{
  temp<-paste(temp,substring(GeneticCode,5,5)[which(substring(GeneticCode,1,3)==substring(Pattern,1+(i-1)*3,i*3))],sep="")
}
#####
Pattern<-readLines("paste.txt")
Pattern<-chartr("ACGT","ACGU",Pattern)
goal<-"PAGQIEYSQ"
numberforpattern<-nchar(Pattern)-nchar(goal)*3+1
TempForPattern<-rep(NA,numberforpattern)
for (i in 1:numberforpattern)
{
  temp<-substring(Pattern,i,i+nchar(goal)*3-1)
  for(j in 1:nchar(goal))
  {
    TempForPattern[i]<-paste(TempForPattern[i],substring(GeneticCode,5,5)[which(substring(GeneticCode,1,3)==substring(temp,1+(j-1)*3,j*3))],sep="")
  }
}
for(i in 1:numberforpattern)
{
  TempForPattern[i]<-substring(TempForPattern[i],3,nchar(TempForPattern[i]))
}
numberforreversecomlement<-nchar(Pattern)-nchar(goal)*3+1
TempForReverseComplement<-rep(NA,numberforreversecomlement)
for (i in 1:numberforreversecomlement)
{
  temp<-substring(Pattern,i,i+nchar(goal)*3-1)
  temp<-chartr("ACGU","UGCA",temp)
  temp<-paste(rev(unlist(strsplit(temp,split=""))),collapse="")
  for(j in 1:nchar(goal))
  {
    TempForReverseComplement[i]<-paste(TempForReverseComplement[i],substring(GeneticCode,5,5)[which(substring(GeneticCode,1,3)==substring(temp,1+(j-1)*3,j*3))],sep="")
  }
}
for(i in 1:numberforreversecomlement)
{
  TempForReverseComplement[i]<-substring(TempForReverseComplement[i],3,nchar(TempForReverseComplement[i]))
}
vector<-rep(0,numberforreversecomlement)
vector[which(TempForPattern==goal)]<-1
vector[which(TempForReverseComplement==goal)]<-1
location<-which(vector==1)
result<-NULL
for (i in 1:length(location))
{
  result<-append(result,substring(Pattern,location[i],location[i]+3*nchar(goal)-1))
}
result<-chartr("ACGU","ACGT",result)
cat(result,file="outfile.txt",sep='\n')
##########
setwd("C:/Users/Nan Sun/Desktop")
IntegerMass<-readLines("integer_mass_table (1).txt")
Peptide<-"NQEL"
PrefixMassFunction<-function(Peptide,IntegerMass)
{
  PrefixMass<-rep(NA,nchar(Peptide)+1)
  PrefixMass[1]<-0
  for (i in 1:nchar(Peptide))
  {
    for (j in 1:20)
    {
      if (substring(IntegerMass[j],1,1)==substring(Peptide,i,i))
        PrefixMass[i+1]<-PrefixMass[i]+as.numeric(substring(IntegerMass[j],3,nchar(IntegerMass[j])))
    }
  }
  return(PrefixMass)
}
PrefixMass<-PrefixMassFunction(Peptide,IntegerMass)
CyclicSpectrumFunction<-function(Peptide,PrefixMass)
{
  CyclicSpectrum<-NULL
  for (i in 1:nchar(Peptide))
  {
    for (j in (i+1):(nchar(Peptide)+1))
    {
      CyclicSpectrum<-append(CyclicSpectrum,PrefixMass[j]-PrefixMass[i])
      if (i>1 & j<(nchar(Peptide)+1))
      {CyclicSpectrum<-append(CyclicSpectrum,PrefixMass[1+nchar(Peptide)]-(PrefixMass[j]-PrefixMass[i]))}
    }
  }
  return(sort(append(CyclicSpectrum,0)))
}
CyclicSpectrum<-CyclicSpectrumFunction(Peptide,PrefixMass)
#####
LinearSpectrumFunction<-function(Peptide,PrefixMass)
{
  LinearSpectrum<-NULL
  for (i in 1:nchar(Peptide))
  {
    for (j in (i+1):(nchar(Peptide)+1))
    {
      LinearSpectrum<-append(LinearSpectrum,PrefixMass[j]-PrefixMass[i])
      if (i>1 & j<(nchar(Peptide)+1))
      {LinearSpectrum<-append(LinearSpectrum,PrefixMass[j]-PrefixMass[i])}
    }
  }
  return(sort(append(LinearSpectrum,0)))
}
LinearSpectrum<-LinearSpectrumFunction(Peptide,PrefixMass)
##########
Spectrum<-readLines("dataset_100_6 (5).txt")
Spectrum<-as.numeric(strsplit(Spectrum,split=" ")[[1]])
result<-NULL
Peptides<-substring(IntegerMass,1,1)
while(length(Peptides)!=0)
{
  Peptides<-as.vector(outer(Peptides,substring(IntegerMass,1,1),paste,sep=""))
  for (j in 1:length(Peptides))
  {
    MassPeptide<-PrefixMassFunction(Peptides[j],IntegerMass)[1+nchar(Peptides[j])]
    CyclicSpectrumPeptide<-CyclicSpectrumFunction(Peptides[j],PrefixMassFunction(Peptides[j],IntegerMass))
    LinearSpectrumPeptide<-LinearSpectrumFunction(Peptides[j],PrefixMassFunction(Peptides[j],IntegerMass))
    if (MassPeptide==max(Spectrum))
    {
      if (all(length(CyclicSpectrumPeptide)==length(Spectrum)) && all(CyclicSpectrumPeptide==Spectrum))
      {
        result<-append(result,Peptides[j])
      }
      Peptides[j]<-NA
    }
    else if (!all(LinearSpectrumPeptide%in%Spectrum))
    {
      Peptides[j]<-NA
    }
  }
  Peptides<-Peptides[!is.na(Peptides)]
}
newresult<-strsplit(result,split=" ")
for (i in 1:length(result))
{
  for (j in 1:nchar(result[i]))
  {
    temp<-IntegerMass[which(substring(IntegerMass,1,1)==substring(result[i],j,j))]
    newresult[[i]][j+1]<-paste("-",substring(temp,3,nchar(temp)),sep="")
  }
}
for (i in 1:length(newresult))
{
  newresult[[i]]<-paste(unlist(strsplit(newresult[[i]],split=" ")),collapse="")
}
newresult<-unlist(newresult)
for (i in 1:length(newresult))
{
  newresult[i]<-substring(newresult[i],nchar(result[1])+2,nchar(newresult[i]))
}
newresult<-unique(newresult)
cat(newresult,file="outfile.txt",sep=' ')
###############
setwd("C:/Users/Nan Sun/Desktop")
IntegerMass<-readLines("integer_mass_table (1).txt")
Peptide<-"KWMCEHRYPFYREEQIDRVTVKPIYYETALCPWEV"
PrefixMassFunction<-function(Peptide,IntegerMass)
{
  PrefixMass<-rep(NA,nchar(Peptide)+1)
  PrefixMass[1]<-0
  for (i in 1:nchar(Peptide))
  {
    for (j in 1:length(IntegerMass))
    {
      if (substring(IntegerMass[j],1,1)==substring(Peptide,i,i))
        PrefixMass[i+1]<-PrefixMass[i]+as.numeric(substring(IntegerMass[j],3,nchar(IntegerMass[j])))
    }
  }
  return(PrefixMass)
}
PrefixMass<-PrefixMassFunction(Peptide,IntegerMass)
CyclicSpectrumFunction<-function(Peptide,PrefixMass)
{
  CyclicSpectrum<-NULL
  for (i in 1:nchar(Peptide))
  {
    for (j in (i+1):(nchar(Peptide)+1))
    {
      CyclicSpectrum<-append(CyclicSpectrum,PrefixMass[j]-PrefixMass[i])
      if (i>1 & j<(nchar(Peptide)+1))
      {CyclicSpectrum<-append(CyclicSpectrum,PrefixMass[1+nchar(Peptide)]-(PrefixMass[j]-PrefixMass[i]))}
    }
  }
  return(sort(append(CyclicSpectrum,0)))
}
CyclicSpectrum<-CyclicSpectrumFunction(Peptide,PrefixMass)
Spectrum<-readLines("Spectrum.txt")
Spectrum<-as.numeric(strsplit(Spectrum,split=" ")[[1]])
UniqueSpectrum<-unique(Spectrum)
ScoreCompute<-function(Spectrum,UniqueSpectrum,CyclicSpectrum)
{
  score<-0
  for (i in 1:length(UniqueSpectrum))
  {
    score<-score+min(length(Spectrum[Spectrum==UniqueSpectrum[i]]),length(CyclicSpectrum[CyclicSpectrum==UniqueSpectrum[i]]))
  }
  return(score)
}
ScoreCompute(Spectrum,UniqueSpectrum,CyclicSpectrum)
##########
setwd("C:/Users/Nan Sun/Desktop")
IntegerMass<-readLines("integer_mass_table (1).txt")
Peptide<-"PEEP"
Spectrum<-readLines("Spectrum.txt")
Spectrum<-as.numeric(strsplit(Spectrum,split=" ")[[1]])
UniqueSpectrum<-unique(Spectrum)
PrefixMass<-PrefixMassFunction(Peptide,IntegerMass)
LinearSpectrumFunction<-function(Peptide,PrefixMass)
{
  LinearSpectrum<-NULL
  for (i in 1:nchar(Peptide))
  {
    for (j in (i+1):(nchar(Peptide)+1))
    {
      LinearSpectrum<-append(LinearSpectrum,PrefixMass[j]-PrefixMass[i])
      if (i>1 & j<(nchar(Peptide)+1))
      {LinearSpectrum<-append(LinearSpectrum,PrefixMass[j]-PrefixMass[i])}
    }
  }
  return(sort(append(LinearSpectrum,0)))
}
LinearSpectrum<-LinearSpectrumFunction(Peptide,PrefixMass)
CyclicSpectrum<-CyclicSpectrumFunction(Peptide,PrefixMass)
CyclicSpectrum<-CyclicSpectrum[(!CyclicSpectrum%in%LinearSpectrum)==FALSE]
ScoreCompute(Spectrum,UniqueSpectrum,CyclicSpectrum)
#####
setwd("C:/Users/Nan Sun/Desktop")
IntegerMass<-readLines("integer_mass_table (1).txt")
Leaderboard<-c("NQEL","QLEN","NQEL","QLEN","NQEL","QLEN","NQEL","QLEN","NQEL","QLEN")
Spectrum<-readLines("Spectrum.txt")
Spectrum<-as.numeric(strsplit(Spectrum,split=" ")[[1]])
UniqueSpectrum<-unique(Spectrum)
LinearScores<-rep(0,length(Leaderboard))
Trim<-function(Leaderboard,IntegerMass,Spectrum,UniqueSpectrum,N)
{
  for (i in 1:length(Leaderboard))
  {
    Peptide<-Leaderboard[i]
    PrefixMass<-PrefixMassFunction(Peptide,IntegerMass)
    LinearSpectrum<-LinearSpectrumFunction(Peptide,PrefixMass)
    CyclicSpectrum<-CyclicSpectrumFunction(Peptide,PrefixMass)
    CyclicSpectrum<-CyclicSpectrum[(!CyclicSpectrum%in%LinearSpectrum)==FALSE]
    LinearScores[i]<-ScoreCompute(Spectrum,UniqueSpectrum,CyclicSpectrum)
  }
  Leaderboard<-Leaderboard[order(LinearScores,decreasing=TRUE)]
  LinearScores<-LinearScores[order(LinearScores,decreasing=TRUE)]
  for (j in (N+1):length(Leaderboard))
  {
    if (LinearScores[j]<LinearScores[N])
    {
      Leaderboard[j:length(Leaderboard)]<-NA
      Leaderboard<-Leaderboard[!is.na(Leaderboard)]
      return(Leaderboard)
    }
  }
  return(Leaderboard)
}
#####
N<-5
Trim(Leaderboard,IntegerMass,Spectrum,UniqueSpectrum,N)
###############
setwd("C:/Users/Nan Sun/Desktop")
Spectrum<-readLines("dataset_102_8 (26).txt")
Spectrum<-as.numeric(strsplit(Spectrum,split=" ")[[1]])
UniqueSpectrum<-unique(Spectrum)
Leaderboard<-substring(IntegerMass,1,1)
LeaderPeptide<-sample(Leaderboard,1)
N<-50
while (nchar(LeaderPeptide)==1)
{
  Leaderboard<-as.vector(outer(Leaderboard,substring(IntegerMass,1,1),paste,sep=""))
  for (j in 1:length(Leaderboard))
  {
    MassPeptide<-PrefixMassFunction(Leaderboard[j],IntegerMass)[1+nchar(Leaderboard[j])]
    CyclicSpectrumPeptide<-CyclicSpectrumFunction(Leaderboard[j],PrefixMassFunction(Leaderboard[j],IntegerMass))
    ScorePeptide<-ScoreCompute(Spectrum,UniqueSpectrum,CyclicSpectrumPeptide)
    CyclicSpectrumLeaderPeptide<-CyclicSpectrumFunction(LeaderPeptide,PrefixMassFunction(LeaderPeptide,IntegerMass))
    ScoreLeaderPeptide<-ScoreCompute(Spectrum,UniqueSpectrum,CyclicSpectrumLeaderPeptide)
    if (MassPeptide==max(Spectrum))
    {
      if (ScorePeptide>ScoreLeaderPeptide)
      {
        LeaderPeptide<-Leaderboard[j]
      }
    }
    else if (MassPeptide>max(Spectrum))
    {
      Leaderboard[j]<-NA
    }
  }
  Leaderboard<-Leaderboard[!is.na(Leaderboard)]
  LinearScores<-rep(0,length(Leaderboard))
  Leaderboard<-Trim(Leaderboard,IntegerMass,Spectrum,UniqueSpectrum,N)
}
result<-LeaderPeptide
newresult<-strsplit(result,split=" ")
for (i in 1:length(result))
{
  for (j in 1:nchar(result[i]))
  {
    temp<-IntegerMass[which(substring(IntegerMass,1,1)==substring(result[i],j,j))]
    newresult[[i]][j+1]<-paste("-",substring(temp,3,nchar(temp)),sep="")
  }
}
for (i in 1:length(newresult))
{
  newresult[[i]]<-paste(unlist(strsplit(newresult[[i]],split=" ")),collapse="")
}
newresult<-unlist(newresult)
for (i in 1:length(newresult))
{
  newresult[i]<-substring(newresult[i],nchar(result[1])+2,nchar(newresult[i]))
}
newresult<-unique(newresult)
cat(newresult,file="outfile.txt",sep=' ')
##########
setwd("C:/Users/Nan Sun/Desktop")
Spectrum<-readLines("value.txt")
Spectrum<-as.numeric(strsplit(Spectrum,split=" ")[[1]])
UniqueSpectrum<-unique(Spectrum)
vector<-sort(Spectrum)
matrix<-abs(outer(vector,vector,'-'))
convolution<-matrix[upper.tri(matrix,diag=FALSE)]
convolution<-convolution[convolution!=0]
#cat(convolution,file="outfile.txt",sep=' ')
convolution<-convolution[order(convolution,decreasing=TRUE)]
convolution<-convolution[convolution>=57]
convolution<-convolution[convolution<=200]
M<-17
Store<-sort(table(convolution),decreasing=TRUE)
NewConvolution<-convolution[convolution%in%as.numeric(names(Store[1:max(which(Store==Store[M][[1]]))]))]
UniqueNewConvolution<-unique(NewConvolution)
IntegerMass<-paste(letters[seq(1:17)],UniqueNewConvolution[1:17],sep=" ")
IntegerMass<-paste(letters[seq(1:26)],UniqueNewConvolution[1:26],sep=" ")
IntegerMass<-append(IntegerMass,paste(LETTERS[seq(1:(length(UniqueNewConvolution)-26))],UniqueNewConvolution[27:length(UniqueNewConvolution)],sep=" "))
Leaderboard<-substring(IntegerMass,1,1)
LeaderPeptide<-sample(Leaderboard,1)
N<-60
while (length(Leaderboard)!=0)
{
  Leaderboard<-as.vector(outer(Leaderboard,substring(IntegerMass,1,1),paste,sep=""))
  for (j in 1:length(Leaderboard))
  {
    MassPeptide<-PrefixMassFunction(Leaderboard[j],IntegerMass)[1+nchar(Leaderboard[j])]
    CyclicSpectrumPeptide<-CyclicSpectrumFunction(Leaderboard[j],PrefixMassFunction(Leaderboard[j],IntegerMass))
    ScorePeptide<-ScoreCompute(Spectrum,UniqueSpectrum,CyclicSpectrumPeptide)
    CyclicSpectrumLeaderPeptide<-CyclicSpectrumFunction(LeaderPeptide,PrefixMassFunction(LeaderPeptide,IntegerMass))
    ScoreLeaderPeptide<-ScoreCompute(Spectrum,UniqueSpectrum,CyclicSpectrumLeaderPeptide)
    if (MassPeptide==max(Spectrum))
    {
      if (ScorePeptide>ScoreLeaderPeptide)
      {
        LeaderPeptide<-Leaderboard[j]
      }
    }
    else if (MassPeptide>max(Spectrum))
    {
      Leaderboard[j]<-NA
    }
  }
  Leaderboard<-Leaderboard[!is.na(Leaderboard)]
  LinearScores<-rep(0,length(Leaderboard))
  Leaderboard<-Trim(Leaderboard,IntegerMass,Spectrum,UniqueSpectrum,N)
}
###############
result<-LeaderPeptide
newresult<-strsplit(result,split=" ")
for (i in 1:length(result))
{
  for (j in 1:nchar(result[i]))
  {
    temp<-IntegerMass[which(substring(IntegerMass,1,1)==substring(result[i],j,j))]
    newresult[[i]][j+1]<-paste("-",substring(temp,3,nchar(temp)),sep="")
  }
}
for (i in 1:length(newresult))
{
  newresult[[i]]<-paste(unlist(strsplit(newresult[[i]],split=" ")),collapse="")
}
newresult<-unlist(newresult)
for (i in 1:length(newresult))
{
  newresult[i]<-substring(newresult[i],nchar(result[1])+2,nchar(newresult[i]))
}
newresult<-unique(newresult)
cat(newresult,file="outfile.txt",sep=' ')
###############
DPChange<-function(money,Coins)
{
  MinNumCoins<-rep(0,money+1)
  for (m in 1:money)
  {
    MinNumCoins[m+1]<-Inf
    for (i in 1:length(Coins))
    {
      if (m>=Coins[i])
      {
        if ((MinNumCoins[m+1-Coins[i]]+1)<MinNumCoins[m+1])
        {
          MinNumCoins[m+1]<-MinNumCoins[m+1-Coins[i]]+1
        }
      }
    }
  }
  return(MinNumCoins[money+1])
}
money<-17073
Coins<-c(9,8,5,3,1)
DPChange(money,Coins)
###############
ManhattanTourist<-function(n,m,Down,Right)
{
  S<-matrix(0,nrow=n+1,ncol=m+1)
  for (i in 1:n)
  {
    S[i+1,1]<-S[i-1+1,1]+Down[i,1]
  }
  for (j in 1:m)
  {
    S[1,j+1]<-S[1,j-1+1]+Right[1,j]
  }
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      S[i+1,j+1]<-max((S[i,j+1]+Down[i,j+1]),(S[i+1,j]+Right[i+1,j]))
    }
  }
  return(S[n+1,m+1])
}
setwd("C:/Users/Nan Sun/Desktop")
Down<-read.table("Down.txt",header=FALSE,sep=" ")
Right<-read.table("Right.txt",header=FALSE,sep=" ")
n<-dim(Down)[1]
m<-dim(Right)[2]
ManhattanTourist(n,m,Down,Right)
###############
LCSBackTrack<-function(V,W)
{
  S<-matrix(0,nrow=n+1,ncol=m+1)
  BackTrack<-matrix(NA,nrow=n,ncol=m)
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      temp<-S[i-1+1,j-1+1]
      if (substring(V,i,i)==substring(W,j,j))
      {
        temp<-1+temp
      }
      S[i+1,j+1]<-max(S[i-1+1,j+1],S[i+1,j-1+1],temp)
      if (S[i+1,j+1]==S[i-1+1,j+1])
      {
        BackTrack[i,j]<-"down"
      }
      else if (S[i+1,j+1]==S[i+1,j-1+1])
      {
        BackTrack[i,j]<-"right"
      }
      else if ((S[i+1,j+1]==S[i-1+1,j-1+1]+1) & (substring(V,i,i)==substring(W,j,j)))
      {
        BackTrack[i,j]<-"diagonal"
      }
    }
  }
  return(BackTrack)
}
setwd("C:/Users/Nan Sun/Desktop")
V<-readLines("V.txt")
W<-readLines("W.txt")
V<-"CTCGAT"
W<-"TACGTC"
n<-nchar(V)
m<-nchar(W)
BackTrack<-LCSBackTrack(V,W)
###############
OutputLCS<-function(BackTrack,V,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS(BackTrack,V,i-1,j)
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS(BackTrack,V,i,j-1)
  }
  else 
  {
    OutputLCS(BackTrack,V,i-1,j-1)
    sink('outfile.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
}
OutputLCS(BackTrack,V,n,m)
###############
setwd("C:/Users/Nan Sun/Desktop")
Graph<-read.table("Graph.txt",header=FALSE,sep="\n")
Graph<-as.matrix(Graph)
Graph<-unlist(strsplit(Graph,split=("->")))
NewGraph<-matrix(NA,nrow=length(Graph)/2,ncol=1)
for (i in 1:dim(NewGraph)[1])
{
  NewGraph[i]<-paste(Graph[2*i-1],Graph[2*i],collapse=" ")
}
NewGraph<-unlist(strsplit(NewGraph,split=(":")))
NewNewGraph<-matrix(NA,nrow=length(NewGraph)/2,ncol=1)
for (i in 1:dim(NewNewGraph)[1])
{
  NewNewGraph[i]<-paste(NewGraph[2*i-1],NewGraph[2*i],collapse=" ")
}
Nodes<-strsplit(NewNewGraph,split=(" "))
LeftNodes<-matrix(NA,nrow=length(Nodes),ncol=1)
for (i in 1:dim(LeftNodes)[1])
{
  LeftNodes[i]<-Nodes[[i]][1]
}
RightNodes<-matrix(NA,nrow=length(Nodes),ncol=1)
for (i in 1:dim(RightNodes)[1])
{
  RightNodes[i]<-Nodes[[i]][2]
}
#####
AllNodes<-sort(append(unique(as.numeric(RightNodes)),min(LeftNodes)))
for (i in 2:length(AllNodes))
{
  predecessors<-LeftNodes[which(RightNodes==AllNodes[i])]
  predecessors<-predecessors[predecessors%in%AllNodes]
  if (length(predecessors)==0)
  {
    AllNodes[which(AllNodes==AllNodes[i])]<-NA
  }
}
AllNodes<-AllNodes[!is.na(AllNodes)]
WeightMatrix<-matrix(0,nrow=length(AllNodes),ncol=length(AllNodes))
for (i in 1:length(Nodes))
{
  s<-match(Nodes[[i]][1],AllNodes)
  t<-match(Nodes[[i]][2],AllNodes)
  WeightMatrix[s,t]<-as.numeric(Nodes[[i]][3])
}
##########
LongestPathBackTracking<-function(AllNodes,sink,WeightMatrix)
{
  score<-rep(0,length(AllNodes))
  recorder<-rep(NA,length(AllNodes))
  for (i in 2:length(AllNodes))
  {
    predecessors<-LeftNodes[which(RightNodes==AllNodes[i])]
    predecessors<-predecessors[predecessors%in%AllNodes]
    if (length(predecessors)!=0)
    {
      for (j in 1:length(predecessors))
      {
        MaxOnePred<-score[which(AllNodes==predecessors[j])]+WeightMatrix[which(AllNodes==predecessors[j]),i]
        if (MaxOnePred>score[i])
        {
          score[i]<-MaxOnePred
          recorder[i]<-predecessors[j]
        }
      }
    }
  }
  return(list(score[which(AllNodes==sink)],recorder))
}
sink<-7
result<-LongestPathBackTracking(AllNodes,sink,WeightMatrix)
LongestLength<-result[[1]]
recorder<-result[[2]]
reverseoutput<-"21"
while(length(recorder)>1)
{
  reverseoutput<-paste(reverseoutput,recorder[length(recorder)],collapse=" ")
  recorder<-recorder[1:which(recorder[length(recorder)]==AllNodes)]
} 
paste(rev(strsplit(reverseoutput,split=" ")[[1]]),collapse="->")
###############
setwd("C:/Users/Nan Sun/Desktop")
ScoringMatrix<-read.table("ScoringMatrix.txt")
V<-readLines("V.txt")
W<-readLines("W.txt")
n<-nchar(V)
m<-nchar(W)
sigma<-(-5)
LCSBackTrack_ScoreMatrix<-function(V,W,ScoringMatrix,sigma)
{
  S<-matrix(0,nrow=n+1,ncol=m+1)
  S[1,]<-seq(from=0,to=sigma*(ncol(S)-1),by=sigma)
  S[,1]<-seq(from=0,to=sigma*(nrow(S)-1),by=sigma)
  BackTrack<-matrix(NA,nrow=n,ncol=m)
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      temp1<-S[i-1+1,j-1+1]
      temp1<-temp1+ScoringMatrix[which((rownames(ScoringMatrix))==substring(V,i,i)),which(colnames(ScoringMatrix)==substring(W,j,j))]
      temp2<-S[i-1+1,j+1]+sigma
      temp3<-S[i+1,j-1+1]+sigma
      S[i+1,j+1]<-max(temp1,temp2,temp3)
      if (S[i+1,j+1]==temp2)
      {
        BackTrack[i,j]<-"down"
      }
      if (S[i+1,j+1]==temp3)
      {
        BackTrack[i,j]<-"right"
      }
      if (S[i+1,j+1]==temp1)
      {
        BackTrack[i,j]<-"diagonal"
      }
    }
  }
  return(BackTrack)
}
BackTrack<-LCSBackTrack_ScoreMatrix(V,W,ScoringMatrix,sigma)
##########
OutputLCS_ScoreMatrix<-function(BackTrack,V,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i,j-1)
    sink('outfile1.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else 
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j-1)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
}
OutputLCS_ScoreMatrix(BackTrack,V,n,m)
##########
OutputLCS_ScoreMatrix<-function(BackTrack,W,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j)
    sink('outfile2.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
  else 
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
}
OutputLCS_ScoreMatrix(BackTrack,W,n,m)
##########
setwd("C:/Users/Nan Sun/Desktop")
ScoringMatrix<-read.table("ScoringMatrix.txt")
V<-readLines("V.txt")
W<-readLines("W.txt")
n<-nchar(V)
m<-nchar(W)
sigma<-(-5)
LCSBackTrack_ScoreMatrix<-function(V,W,ScoringMatrix,sigma)
{
  S<-matrix(0,nrow=n+1,ncol=m+1)
  S[1,]<-seq(from=0,to=sigma*(ncol(S)-1),by=sigma)
  S[,1]<-seq(from=0,to=sigma*(nrow(S)-1),by=sigma)
  BackTrack<-matrix(NA,nrow=n,ncol=m)
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      temp1<-S[i-1+1,j-1+1]
      temp1<-temp1+ScoringMatrix[which((rownames(ScoringMatrix))==substring(V,i,i)),which(colnames(ScoringMatrix)==substring(W,j,j))]
      temp2<-S[i-1+1,j+1]+sigma
      temp3<-S[i+1,j-1+1]+sigma
      temp4<-0
      S[i+1,j+1]<-max(temp1,temp2,temp3,temp4)
      if (S[i+1,j+1]==temp4)
      {
        BackTrack[i,j]<-"skip"
      }
      if (S[i+1,j+1]==temp2)
      {
        BackTrack[i,j]<-"down"
      }
      if (S[i+1,j+1]==temp3)
      {
        BackTrack[i,j]<-"right"
      }
      if (S[i+1,j+1]==temp1)
      {
        BackTrack[i,j]<-"diagonal"
      }
    }
  }
  S[n+1,m+1]<-max(S)
  return(BackTrack)
}
BackTrack<-BackTrack[1:(which(S==max(S),arr.ind=TRUE)[1,]-1)[1],1:(which(S==max(S),arr.ind=TRUE)[1,]-1)[2]]
##########
I<-dim(BackTrack)[1]
J<-dim(BackTrack)[2]
OutputLCS_ScoreMatrix<-function(BackTrack,V,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
  if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i,j-1)
    sink('outfile1.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  if (BackTrack[i,j]=="diagonal") 
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j-1)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
  if (BackTrack[i,j]=="skip") 
  {
    return()
  }
}
OutputLCS_ScoreMatrix(BackTrack,V,I,J)
##########
OutputLCS_ScoreMatrix<-function(BackTrack,W,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j)
    sink('outfile2.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
  if (BackTrack[i,j]=="diagonal") 
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
  if (BackTrack[i,j]=="skip") 
  {
    return()
  }
}
OutputLCS_ScoreMatrix(BackTrack,W,I,J)
###############
V<-readLines("V.txt")
W<-readLines("W.txt")
n<-nchar(V)
m<-nchar(W)
LCSBackTrack<-function(V,W)
{
  S<-matrix(0,nrow=n+1,ncol=m+1)
  S[1,]<-seq(from=0,to=(-1)*(ncol(S)-1),by=-1)
  S[,1]<-seq(from=0,to=(-1)*(nrow(S)-1),by=-1)
  BackTrack<-matrix(NA,nrow=n,ncol=m)
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      temp<-S[i-1+1,j-1+1]
      if (substring(V,i,i)==substring(W,j,j))
      {
        temp<-temp
      }
      if (substring(V,i,i)!=substring(W,j,j))
      {
        temp<-temp-1
      }
      S[i+1,j+1]<-max(S[i-1+1,j+1]-1,S[i+1,j-1+1]-1,temp)
      if (S[i+1,j+1]==S[i-1+1,j+1]-1)
      {
        BackTrack[i,j]<-"down"
      }
      else if (S[i+1,j+1]==S[i+1,j-1+1]-1)
      {
        BackTrack[i,j]<-"right"
      }
      else if ((S[i+1,j+1]==temp))
      {
        BackTrack[i,j]<-"diagonal"
      }
    }
  }
  return(BackTrack)
}
BackTrack<-LCSBackTrack_ScoreMatrix(V,W,ScoringMatrix,sigma)
##########
OutputLCS_ScoreMatrix<-function(BackTrack,V,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i,j-1)
    sink('outfile1.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else 
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j-1)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
}
OutputLCS_ScoreMatrix(BackTrack,V,n,m)
##########
OutputLCS_ScoreMatrix<-function(BackTrack,W,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j)
    sink('outfile2.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
  else 
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
}
OutputLCS_ScoreMatrix(BackTrack,W,n,m)
Genome1<-read.table("outfile1.txt")
Genome1<-toString(Genome1$V1)
Genome2<-read.table("outfile2.txt")
Genome2<-toString(Genome2$V1)
HammingDistance<-function(Genome1,Genome2)
{
  count<-0
  for (i in 1:nchar(Genome1))
  {
    if(substring(Genome1,i,i)!=substring(Genome2,i,i))
      count<-count+1
    else count<-count
  }
  return(count)
}
HammingDistance(Genome1,Genome2)
####################Overlap Alignment
V<-"GATACAGCACACTAGTACTACTTGAC"
W<-"CTA-TAGT-CTTAACGATATACGAC"
n<-nchar(V)
m<-nchar(W)
LCSBackTrack<-function(V,W)
{
  S<-matrix(0,nrow=n+1,ncol=m+1)
  S[1,]<-seq(from=0,to=(-2)*(ncol(S)-1),by=-1)
  S[,1]<-0
  BackTrack<-matrix(NA,nrow=n,ncol=m)
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      temp<-S[i-1+1,j-1+1]
      if (substring(V,i,i)==substring(W,j,j))
      {
        temp<-temp+1
      }
      if (substring(V,i,i)!=substring(W,j,j))
      {
        temp<-temp
      }
      S[i+1,j+1]<-max(S[i-1+1,j+1]-2,S[i+1,j-1+1]-2,temp)
      if (S[i+1,j+1]==S[i-1+1,j+1]-2)
      {
        BackTrack[i,j]<-"down"
      }
      else if (S[i+1,j+1]==S[i+1,j-1+1]-2)
      {
        BackTrack[i,j]<-"right"
      }
      else if ((S[i+1,j+1]==temp))
      {
        BackTrack[i,j]<-"diagonal"
      }
    }
  }
  max(S[,dim(S)[2]])
  return(BackTrack)
}
BackTrack<-BackTrack[1:((which(S[,dim(S)[2]]==max(S[,dim(S)[2]]),arr.ind=TRUE))[1]-1),1:(dim(S)[2]-1)]
I<-dim(BackTrack)[1]
J<-dim(BackTrack)[2]
OutputLCS_ScoreMatrix<-function(BackTrack,V,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i,j-1)
    sink('outfile1.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else 
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j-1)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
}
OutputLCS_ScoreMatrix(BackTrack,V,I,J)
##########
OutputLCS_ScoreMatrix<-function(BackTrack,W,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j)
    sink('outfile2.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
  else 
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
}
OutputLCS_ScoreMatrix(BackTrack,W,I,J)
####################Fitting Alignment
V<-"GTTGGATTACGAATCGATATCTGTTTG"
W<-"ACGTCG"
n<-nchar(V)
m<-nchar(W)
LCSBackTrack<-function(V,W)
{
  S<-matrix(0,nrow=n+1,ncol=m+1)
  BackTrack<-matrix(NA,nrow=n,ncol=m)
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      temp<-S[i-1+1,j-1+1]
      if (substring(V,i,i)==substring(W,j,j))
      {
        temp<-temp+1
      }
      if (substring(V,i,i)!=substring(W,j,j))
      {
        temp<-temp-1
      }
      S[i+1,j+1]<-max(S[i-1+1,j+1]-1,S[i+1,j-1+1]-1,temp)
      if (S[i+1,j+1]==S[i-1+1,j+1]-1)
      {
        BackTrack[i,j]<-"down"
      }
      else if (S[i+1,j+1]==S[i+1,j-1+1]-1)
      {
        BackTrack[i,j]<-"right"
      }
      else if ((S[i+1,j+1]==temp))
      {
        BackTrack[i,j]<-"diagonal"
      }
    }
  }
  max(S[dim(S)[1],])
  return(BackTrack)
}
BackTrack<-BackTrack[1:(dim(S)[1]-1),1:((which(S[dim(S)[1],]==max(S[dim(S)[1],]),arr.ind=TRUE))[1]-1)]
I<-dim(BackTrack)[1]
J<-dim(BackTrack)[2]
OutputLCS_ScoreMatrix<-function(BackTrack,V,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i,j-1)
    sink('outfile1.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else 
  {
    OutputLCS_ScoreMatrix(BackTrack,V,i-1,j-1)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
}
OutputLCS_ScoreMatrix(BackTrack,V,I,J)
##########
OutputLCS_ScoreMatrix<-function(BackTrack,W,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (BackTrack[i,j]=="down")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j)
    sink('outfile2.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else if (BackTrack[i,j]=="right")
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
  else 
  {
    OutputLCS_ScoreMatrix(BackTrack,W,i-1,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
}
OutputLCS_ScoreMatrix(BackTrack,W,I,J)
####################
setwd("C:/Users/Nan Sun/Desktop")
ScoringMatrix<-read.table("ScoringMatrix.txt")
V<-"GQCPTPYPGCWASDSDKGWYAECETPDVHMQLEGHHLFCKDFVACTHRWGKHGSFIHDYYDETALNYFTRIK"
W<-"GQYNTGVWASDSDCSNHGPDTGWYDECETPDVHMQLEGHMLFCKQESNKDFVACTHAWGKHGSFIHDYYDETATYYDTRIK"
n<-nchar(V)
m<-nchar(W)
Lower_S<-matrix(0,nrow=n+1,ncol=m+1)
Middle_S<-matrix(0,nrow=n+1,ncol=m+1)
Upper_S<-matrix(0,nrow=n+1,ncol=m+1)
sigma<-11
epsilon<-1
for (i in 2:(n+1))
{
  Middle_S[i,1]<--sigma-epsilon*(i-2)
}
for (j in 2:(m+1))
{
  Middle_S[1,j]<--sigma-epsilon*(j-2)
}
for (j in 1:(m+1))
{
  Lower_S[1,j]<--Inf
}
for (i in 1:(n+1))
{
  Upper_S[i,1]<--Inf
}
S<-Middle_S
#####
for (i in 1:n)
{
 for (j in 1:m)
  {
   temp<-ScoringMatrix[which((rownames(ScoringMatrix))==substring(V,i,i)),which(colnames(ScoringMatrix)==substring(W,j,j))]
   Middle_S[i+1,j+1]<-max(Lower_S[i-1+1,j-1+1],Middle_S[i-1+1,j-1+1],Upper_S[i-1+1,j-1+1])+temp
   Lower_S[i+1,j+1]<-max(Lower_S[i-1+1,j+1]-epsilon,Middle_S[i-1+1,j+1]-sigma)
   Upper_S[i+1,j+1]<-max(Upper_S[i+1,j-1+1]-epsilon,Middle_S[i+1,j-1+1]-sigma)
   S[i+1,j+1]<-max(Lower_S[i+1,j+1],Middle_S[i+1,j+1],Upper_S[i+1,j+1])
  }
}
UpperScore<-Upper_S[2:(n+1),2:(m+1)]
MiddleScore<-Middle_S[2:(n+1),2:(m+1)]
LowerScore<-Lower_S[2:(n+1),2:(m+1)]
SScore<-S[2:(n+1),2:(m+1)]
U_BackTrack<-matrix(NA,nrow=n,ncol=m)
M_BackTrack<-matrix(NA,nrow=n,ncol=m)
L_BackTrack<-matrix(NA,nrow=n,ncol=m)
for (i in 1:(dim(UpperScore)[1]))
{
  for (j in 1:(dim(UpperScore)[2]))
  {
    if(UpperScore[i,j]==SScore[i,j])
    {
      U_BackTrack[i,j]<-"right"
    }
    else
    {
      U_BackTrack[i,j]<-"middle"
    }
  }
}
for (i in 1:(dim(LowerScore)[1]))
{
  for (j in 1:(dim(LowerScore)[2]))
  {
    if(LowerScore[i,j]==SScore[i,j])
    {
      L_BackTrack[i,j]<-"down"
    }
    else
    {
      L_BackTrack[i,j]<-"middle"
    }
  }
}
for (i in 1:(dim(MiddleScore)[1]))
{
  for (j in 1:(dim(MiddleScore)[2]))
  {
    if(MiddleScore[i,j]==SScore[i,j])
    {
      M_BackTrack[i,j]<-"diagonal"
    }
    else if (LowerScore[i,j]>UpperScore[i,j])
    {
      M_BackTrack[i,j]<-"lower"
    }
    else
    {
      M_BackTrack[i,j]<-"upper"
    }
  }
}
##########
Output_ScoreMatrix_V<-function(U_BackTrack,M_BackTrack,L_BackTrack,V,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (L_BackTrack[i,j]=="down")
  {
    Output_ScoreMatrix_V(U_BackTrack,M_BackTrack,L_BackTrack,V,i-1,j)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
  else if (U_BackTrack[i,j]=="right")
  {
    Output_ScoreMatrix_V(U_BackTrack,M_BackTrack,L_BackTrack,V,i,j-1)
    sink('outfile1.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else if (M_BackTrack[i,j]=="diagonal")
  {
    Output_ScoreMatrix_V(U_BackTrack,M_BackTrack,L_BackTrack,V,i-1,j-1)
    sink('outfile1.txt',append=TRUE)
    cat(substring(V,i,i),sep="")
    sink()
  }
}  
Output_ScoreMatrix_V(U_BackTrack,M_BackTrack,L_BackTrack,V,n,m)
Output_ScoreMatrix_W<-function(U_BackTrack,M_BackTrack,L_BackTrack,W,i,j)
{
  if ((i==0) | (j==0))
  {
    return()
  }
  if (L_BackTrack[i,j]=="down")
  {
    Output_ScoreMatrix_W(U_BackTrack,M_BackTrack,L_BackTrack,W,i-1,j)
    sink('outfile2.txt',append=TRUE)
    cat("-",sep="")
    sink()
  }
  else if (U_BackTrack[i,j]=="right")
  {
    Output_ScoreMatrix_W(U_BackTrack,M_BackTrack,L_BackTrack,W,i,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
  else if (M_BackTrack[i,j]=="diagonal") 
  {
    Output_ScoreMatrix_W(U_BackTrack,M_BackTrack,L_BackTrack,W,i-1,j-1)
    sink('outfile2.txt',append=TRUE)
    cat(substring(W,j,j),sep="")
    sink()
  }
}
Output_ScoreMatrix_W(U_BackTrack,M_BackTrack,L_BackTrack,W,n,m)
####################
setwd("C:/Users/Nan Sun/Desktop")
ScoringMatrix<-read.table("ScoringMatrix.txt")
V<-readLines("V.txt")
W<-readLines("W.txt")
n<-nchar(V)
m<-nchar(W)
middle<-floor(m/2)
S<-matrix(0,nrow=n+1,ncol=middle+1)
for (i in 2:(n+1))
{
  S[i,1]<--5*(i-1)
}
for (j in 2:(middle+1))
{
  S[1,j]<--5*(j-1)
}
for (i in 1:n)
{
  for (j in 1:middle)
  {
    temp<-ScoringMatrix[which((rownames(ScoringMatrix))==substring(V,i,i)),which(colnames(ScoringMatrix)==substring(W,j,j))]
    S[i+1,j+1]<-max(S[i,j+1]-5,S[i,j]+temp,S[i+1,j]-5)
  }
}
V_rev<-paste(rev(strsplit(V,NULL)[[1]]),collapse="")
W_rev<-paste(rev(strsplit(W,NULL)[[1]]),collapse="")
S_rev<-matrix(0,nrow=n+1,ncol=m-middle+1)
for (i in 2:(n+1))
{
  S_rev[i,1]<--5*(i-1)
}
for (j in 2:(m-middle+1))
{
  S_rev[1,j]<--5*(j-1)
}
for (i in 1:n)
{
  for (j in 1:(m-middle))
  {
    temp<-ScoringMatrix[which((rownames(ScoringMatrix))==substring(V_rev,i,i)),which(colnames(ScoringMatrix)==substring(W_rev,j,j))]
    S_rev[i+1,j+1]<-max(S_rev[i,j+1]-5,S_rev[i,j]+temp,S_rev[i+1,j]-5)
  }
}
middle_length<-S[,middle+1]+rev(S_rev[,m-middle+1])
middle_node_x<-which.max(middle_length)
middle_node_y<-middle+1
middle_node<-paste("(",middle_node_x,", ",middle_node_y,")",sep="")
FindMiddleEdge<-function(middle_node_x,middle_node_y,S)
{
  temp<-ScoringMatrix[which((rownames(ScoringMatrix))==substring(V,middle_node_x-1,middle_node_x-1)),which(colnames(ScoringMatrix)==substring(W,middle_node_y-1,middle_node_y-1))]
  if (S[middle_node_x-1,middle_node_y-1]==S[middle_node_x,middle_node_y]-temp)
  {
    middle_edge<-paste("(",middle_node_x-1,", ",middle_node_y-1,") ",middle_node,sep="")
  }
  else if (S[middle_node_x,middle_node_y-1]==S[middle_node_x,middle_node_y]+5)
  {
    middle_edge<-paste("(",middle_node_x,", ",middle_node_y-1,") ",middle_node,sep="")
  }
  else if (S[middle_node_x-1,middle_node_y]==S[middle_node_x,middle_node_y]+5)
  {
    middle_edge<-paste("(",middle_node_x-1,", ",middle_node_y,") ",middle_node,sep="")
  }
  return(middle_edge)
}
FindMiddleEdge(middle_node_x,middle_node_y,S)
###############GreedySorting
setwd("C:/Users/Nan Sun/Desktop")
P<-"-16 -20 +11 +12 -14 -13 -15 -6 -8 -19 -18 -17 -10 +4 -5 -2 +7 -3 +1 -9"
P<-strsplit(P,split=" ")[[1]]
p<-length(P)
for (k in 1:p)
{
  if (substring(P[k],2,nchar(P[k]))!=k)
  {
    temp<-which(substring(P,2,nchar(P[]))==k)
    string_temp<-strsplit(toString(as.numeric(rev(P[k:temp]))*(-1)),split=", ")[[1]]
    string_temp[which(substring(string_temp,1,1)!="-")]<-paste("+",string_temp[which(substring(string_temp,1,1)!="-")],sep="")
    P[k:temp]<-string_temp
    cat(P,file="outfile.txt",sep=' ',append=TRUE)
    cat("\n",file="outfile.txt",append=TRUE)
  }
  if (substring(P[k],1,1)=="-")
  {
    substring(P[k],1,1)<-"+"
    cat(P,file="outfile.txt",sep=' ',append=TRUE)
    cat("\n",file="outfile.txt",append=TRUE)
  }
}
result<-readLines("outfile.txt")
result<-paste("(",result,")",sep="")
write.table(result,file="outfile1.txt",row.names=FALSE,quote=FALSE)
###############
setwd("C:/Users/Nan Sun/Desktop")
P<-readLines("data.txt")
P<-strsplit(P,split=" ")[[1]]
p<-length(P)
P<-as.numeric(P)
P_E<-seq(0,p+1,1)
P_E[2:(p+1)]<-P
no_breakpoints<-0
for (i in 1:(length(P_E)-1))
{
  if ((P_E[i+1]-P_E[i])!=1)
  {
    no_breakpoints<-no_breakpoints+1
  }
}
#######################################
#######################################
#######################################
#####replace "(" with "( "#############
#####then repalce ")" with ") "########
#####finally replace ")" with " )"#####
setwd("C:/Users/Nan Sun/Desktop")
PQChromosomes<-read.table("PQChromosomes.txt",fill=TRUE)
#######################################
PChromosomes<-data.frame(PQChromosomes[1,])
QChromosomes<-data.frame(PQChromosomes[2,])
P_Begin_Index<-which(PChromosomes=="(")
Q_Begin_Index<-which(QChromosomes=="(")
P_End_Index<-which(PChromosomes==")")
Q_End_Index<-which(QChromosomes==")")
P_Number<-length(which(PChromosomes==")"))
Q_Number<-length(which(QChromosomes==")"))
#####Use Number Largest################
Block_Number<-length(which((t(QChromosomes)=="(")==FALSE & (t(QChromosomes)==")")==FALSE))
#######################################
ChromosomeToCycle<-function(Chromosome)
{
  Nodes<-rep(NA,2*length(Chromosome))
  for (j in 1:length(Chromosome))
  {
    i<-Chromosome[j]
    if (i>0)
    {
      Nodes[2*j-1]<-2*i-1
      Nodes[2*j]<-2*i
    }
    else
    {
      Nodes[2*j-1]<--2*i
      Nodes[2*j]<--2*i-1
    }
  }
  return(Nodes)
}
#######################################
#####Nodes<-ChromosomeToCycle(c(+1,-2,-3,+4))
#######################################
CycleToChromosome<-function(Nodes)
{
  Chromosome<-rep(NA,length(Nodes)/2)
  for (j in 1:(length(Nodes)/2))
  {
    if (Nodes[2*j-1]<Nodes[2*j])
    {
      Chromosome[j]<-Nodes[2*j]/2
    }
    else
    {
      Chromosome[j]<--Nodes[2*j-1]/2
    }
  }
  return(Chromosome)
}
########################################
#####CycleToChromosome(Nodes)
########################################
########################################
ColoredEdges<-function(Number,Begin_Index,End_Index,Chromosomes)
{
  Edges<-NULL
  for (i in 1:Number)
  {
    Chromosome<-as.numeric(unlist(t(Chromosomes[(Begin_Index[i]+1):(End_Index[i]-1)])))
    Nodes<-ChromosomeToCycle(Chromosome)
    Nodes<-append(Nodes,Nodes[1])
    for (j in 1:length(Chromosome))
    {
      Edges<-rbind(Edges,c(Nodes[2*j],Nodes[2*j+1]))
    }
  }
  return(Edges)
}
########################################
BlackEdges<-function(Number,Begin_Index,End_Index,Chromosomes)
{
  Edges<-NULL
  for (i in 1:Number)
  {
    Chromosome<-as.numeric(unlist(t(Chromosomes[(Begin_Index[i]+1):(End_Index[i]-1)])))
    Nodes<-ChromosomeToCycle(Chromosome)
    for (j in 1:length(Chromosome))
    {
      Edges<-rbind(Edges,c(Nodes[2*j-1],Nodes[2*j]))
    }
  }
  return(Edges)
}
#########################################
P_ColoredEdges<-ColoredEdges(P_Number,P_Begin_Index,P_End_Index,PChromosomes)
Q_ColoredEdges<-ColoredEdges(Q_Number,Q_Begin_Index,Q_End_Index,QChromosomes)
P_BlackEdges<-BlackEdges(P_Number,P_Begin_Index,P_End_Index,PChromosomes)
Q_BlackEdges<-BlackEdges(Q_Number,Q_Begin_Index,Q_End_Index,QChromosomes)
#########################################
library(plyr)
#########################################
FindFinalCycles<-function(P_ColoredEdges,Q_ColoredEdges)
{
  Cycles<-NULL
  for(kk in 1:nrow(P_ColoredEdges))
  {
    VisitedNumbers<-NULL
    VisitedNumbers<-append(VisitedNumbers,P_ColoredEdges[kk,1])
    CurrentNumber<-P_ColoredEdges[kk,2]
    VisitedNumbers<-append(VisitedNumbers,CurrentNumber)
    while (length(VisitedNumbers)<max(P_ColoredEdges))
    {
      Qtemp<-Q_ColoredEdges[apply(Q_ColoredEdges,1,function(x) CurrentNumber%in%x),]
      CurrentNumber<-Qtemp[which(Qtemp!=CurrentNumber)]
      VisitedNumbers<-append(VisitedNumbers,CurrentNumber)
      Ptemp<-P_ColoredEdges[apply(P_ColoredEdges,1,function(x) CurrentNumber%in%x),]
      CurrentNumber<-Ptemp[which(Ptemp!=CurrentNumber)]
      VisitedNumbers<-append(VisitedNumbers,CurrentNumber)
    }
    Cycles<-rbind(Cycles,VisitedNumbers)
  }
  Cycles<-t(apply(Cycles,1,sort))
  Cycles<-unique.data.frame(Cycles)
  FinalCycles<-data.frame(NULL)
  for (i in 1:nrow(Cycles))
  {
    FinalCycles<-rbind.fill(FinalCycles,data.frame(t(unique(Cycles[i,]))))
  }
  FinalCycles<- unique.data.frame(FinalCycles)
  return(FinalCycles)
}
############################################
#####FindFinalCycles(P_ColoredEdges,Q_ColoredEdges)
############################################
#####GraphToGenome######used#for#P##or##Q###
GraphToGenome<-function(P_BlackEdges,P_ColoredEdges)
{
  Chromosomes<-NULL
  Cycles<-FindFinalCycles(P_BlackEdges,P_ColoredEdges)
  for (i in 1:nrow(Cycles))
  {
    NodesTemp<-Cycles[i,][!is.na(Cycles[i,])]
    VisitedNumbers<-NULL
    InitialTemp<-P_BlackEdges[which(apply(P_BlackEdges,1,function(x) x[1]%in%NodesTemp)==TRUE),]
    InitialTemp<-t(as.matrix(InitialTemp))
    ttttemp<-ifelse(nrow(InitialTemp)==1,min(InitialTemp[2]),min(InitialTemp[,1]))
    VisitedNumbers<-append(VisitedNumbers,ttttemp)
    FirstTemp<-P_BlackEdges[apply(P_BlackEdges,1,function(x) ttttemp%in%x),]
    CurrentNumber<-FirstTemp[which(FirstTemp!=ttttemp)]
    VisitedNumbers<-append(VisitedNumbers,CurrentNumber)
    while (length(VisitedNumbers)<length(NodesTemp))
    {
      Qtemp<-P_ColoredEdges[apply(P_ColoredEdges,1,function(x) CurrentNumber%in%x),]
      CurrentNumber<-Qtemp[which(Qtemp!=CurrentNumber)]
      VisitedNumbers<-append(VisitedNumbers,CurrentNumber)
      Ptemp<-P_BlackEdges[apply(P_BlackEdges,1,function(x) CurrentNumber%in%x),]
      CurrentNumber<-Ptemp[which(Ptemp!=CurrentNumber)]
      VisitedNumbers<-append(VisitedNumbers,CurrentNumber)
    }
    ChromosomeTemp<-CycleToChromosome(VisitedNumbers)
    ChromosomeTemp<-paste("(",paste(ChromosomeTemp, collapse=" "),") ")
    Chromosomes<-append(Chromosomes,ChromosomeTemp)
  }
  return(Chromosomes)
}
###########################################
#####GraphToGenome(P_BlackEdges,P_ColoredEdges)
###########################################
Two_BreakGenomeGraph<-function(P_ColoredEdges,ice,iceice,jack,jackjack)
{
  P_ColoredEdges<-P_ColoredEdges[!apply(P_ColoredEdges,1,function(x) identical(sort(x),sort(c(ice,iceice)))),]
  P_ColoredEdges<-P_ColoredEdges[!apply(P_ColoredEdges,1,function(x) identical(sort(x),sort(c(jack,jackjack)))),]
  P_ColoredEdges<-rbind(P_ColoredEdges,c(ice,jack),c(iceice,jackjack))
  return(P_ColoredEdges)
}
###########################################
#####Two_BreakGenomeGraph(P_ColoredEdges,1,6,3,8)
###########################################
Two_BreakOnGenome<-function(PChromosomes,ice,iceice,jack,jackjack)
{
  PChromosomes<-unlist(strsplit(PChromosomes," "))
  P_Begin_Index<-which(PChromosomes=="(")
  P_End_Index<-which(PChromosomes==")")
  P_Number<-length(which(PChromosomes==")"))
  P_ColoredEdges<-ColoredEdges(P_Number,P_Begin_Index,P_End_Index,PChromosomes)
  P_BlackEdges<-BlackEdges(P_Number,P_Begin_Index,P_End_Index,PChromosomes)
  P_ColoredEdges<-Two_BreakGenomeGraph(P_ColoredEdges,ice,iceice,jack,jackjack)
  NewPChromosomes<-GraphToGenome(P_BlackEdges,P_ColoredEdges)
  return(NewPChromosomes)
}
##############################
##############################
ShortestRearrangementScenario<-function(P_ColoredEdges,P_BlackEdges,Q_ColoredEdges)
{
  RedEdges<-P_ColoredEdges
  BlueEdges<-Q_ColoredEdges
  BreakPointCycles<-FindFinalCycles(RedEdges,BlueEdges)
  while(ncol(BreakPointCycles)>2)
  {
    NonTrivialCycles<-BreakPointCycles[which(is.na(BreakPointCycles[,3])==FALSE),]
    BlueEdges_NTC<-BlueEdges[which(apply(BlueEdges,1,function(x) all(x%in%NonTrivialCycles[1,]))==TRUE),]
    ArbitraryBlueEdge<-BlueEdges_NTC[sample(1:nrow(BlueEdges_NTC),1),]
    jack<-ArbitraryBlueEdge[1]
    iceice<-ArbitraryBlueEdge[2]
    jacktemp<-RedEdges[which(apply(RedEdges,1,function(x) jack%in%x)==TRUE),]
    iceicetemp<-RedEdges[which(apply(RedEdges,1,function(x) iceice%in%x)==TRUE),]
    ice<-jacktemp[which(jacktemp!=jack)]
    jackjack<-iceicetemp[which(iceicetemp!=iceice)]
    RedEdges<-RedEdges[!apply(RedEdges,1,function(x) identical(sort(x),sort(c(ice,jack)))),]
    RedEdges<-RedEdges[!apply(RedEdges,1,function(x) identical(sort(x),sort(c(iceice,jackjack)))),]
    RedEdges<-rbind(RedEdges,c(jack,iceice),c(jackjack,ice))
    BreakPointCycles<-FindFinalCycles(RedEdges,BlueEdges)
    tempresult<-GraphToGenome(P_BlackEdges,RedEdges)
    cat(tempresult,file="outfile.txt",sep=' ',append=TRUE)
    cat("\n",file="outfile.txt",append=TRUE)
  }
}
ShortestRearrangementScenario(P_ColoredEdges,P_BlackEdges,Q_ColoredEdges)
##############################
##############################
setwd("C:/Users/Nan Sun/Desktop")
AdjacencyList<-read.table("al.txt",fill=TRUE)
library(stringr)
AdjacencyList<-str_split_fixed(AdjacencyList$V1,"->",2)
AdjacencyList<-data.frame(AdjacencyList)
AdjacencyList<-data.frame(AdjacencyList$X1,str_split_fixed(AdjacencyList$X2,":",2))
names(AdjacencyList)<-c("from","to","weight")
allnumbers<-as.numeric(as.matrix(unique(AdjacencyList$from)))
distancematrix<-matrix(Inf,length(allnumbers),length(allnumbers))
for (i in 1:length(allnumbers))
{
  for (j in 1:length(allnumbers))
  {
    temp<-as.numeric(as.matrix(AdjacencyList[which(AdjacencyList$from==(i-1) & AdjacencyList$to==(j-1)),]$weight))
    if (length(temp)!=0)
    {
      distancematrix[i,j]<-temp
    }
  }
}
diag(distancematrix)<-0
for (k in 1:length(allnumbers))
{
  for (i in 1:length(allnumbers))
  {
    for (j in 1:length(allnumbers))
    {
      if (distancematrix[i,j]>distancematrix[i,k]+distancematrix[k,j])
      {
        distancematrix[i,j]<-distancematrix[i,k]+distancematrix[k,j]
      }
    }
  }
}
nn<-32
result<-distancematrix[1:nn,1:nn]
for (i in 1:nn)
{
  cat(result[i,],file="outfile.txt",sep=' ',append=TRUE)
  cat("\n",file="outfile.txt",append=TRUE)
}
##############################
##############################
setwd("C:/Users/Nan Sun/Desktop")
distancematrix<-read.table("al.txt",fill=TRUE)
jj<-11
Limb<-function(distancematrix,jj)
{
  limblength_jj<-Inf
  for (i in 1:dim(distancematrix)[1])
  {
    for (k in 1:dim(distancematrix)[1])
    {
      temp<-(distancematrix[i,jj]+distancematrix[jj,k]-distancematrix[i,k])/2
      if (temp<limblength_jj & temp>0)
      {
        limblength_jj<-temp
      }
    }
  }
  return(limblength_jj)
}
Limb(distancematrix,jj)
####################
####################
setwd("C:/Users/Nan Sun/Desktop")
distancematrix<-read.table("al.txt",fill=TRUE)
originaldistancematrix<-read.table("al.txt",fill=TRUE)
colnames(distancematrix)<-seq(0,dim(distancematrix)[2]-1,1)
rownames(distancematrix)<-seq(0,dim(distancematrix)[2]-1,1)
colnames(originaldistancematrix)<-seq(0,dim(originaldistancematrix)[2]-1,1)
rownames(originaldistancematrix)<-seq(0,dim(originaldistancematrix)[2]-1,1)
for (l in 1:(dim(originaldistancematrix)[1]-2))
{
  n<-dim(originaldistancematrix)[1]
  limblength<-Limb(originaldistancematrix,n)
  for (j in 1:(n-1))
  {
    originaldistancematrix[j,n]<-originaldistancematrix[j,n]-limblength
    originaldistancematrix[n,j]<-originaldistancematrix[j,n]
  }
  for (ii in 1:n)
  {
    for (kk in 1:n)
    {
      if (originaldistancematrix[ii,kk]==originaldistancematrix[ii,n]+originaldistancematrix[n,kk] & originaldistancematrix[ii,n]!=0 & originaldistancematrix[n,kk]!=0)
      {
        i<-ii
        k<-kk
      } 
    }
  }
  x<-originaldistancematrix[i,n]
  distancematrix<-rbind(distancematrix,Inf)
  distancematrix<-cbind(distancematrix,Inf)
  coltemp<-colnames(distancematrix)
  coltemp[length(coltemp)]<-length(coltemp)-1
  rowtemp<-rownames(distancematrix)
  rowtemp[length(rowtemp)]<-length(rowtemp)-1
  colnames(distancematrix)<-coltemp
  rownames(distancematrix)<-rowtemp
  diag(distancematrix)<-0
  distancematrix[i,dim(distancematrix)[1]]<-x
  distancematrix[dim(distancematrix)[1],i]<-x
  originaldistancematrix<-originaldistancematrix[-n,-n]
}
originaldistancematrix<-read.table("al.txt",fill=TRUE)
colnames(originaldistancematrix)<-seq(0,dim(originaldistancematrix)[2]-1,1)
rownames(originaldistancematrix)<-seq(0,dim(originaldistancematrix)[2]-1,1)
for (i in (dim(originaldistancematrix)[1]+1):dim(distancematrix)[1])
{
  temp<-as.numeric(as.matrix(sort(distancematrix[i,])[2]))
  tempindex<-which(as.vector(distancematrix[i,])==temp)
  for (j in 1:(i-1))
  {
    if (j!=tempindex & distancematrix[i,j]==Inf)
   {
     distancematrix[i,j]<-as.numeric(as.matrix(distancematrix[tempindex,j]-distancematrix[i,tempindex]))
     distancematrix[j,i]<-distancematrix[i,j]
   }
  }
}
##############################
##############################
########UPGMA#################
##############################
##############################
library(phangorn)  
setwd("C:/Users/Nan Sun/Desktop")
distancematrix<-read.table("al.txt",fill=TRUE)
colnames(distancematrix)<-seq(0,dim(distancematrix)[2]-1,1)
rownames(distancematrix)<-seq(0,dim(distancematrix)[2]-1,1)
result<-upgma(distancematrix,method="average")
UPGMA_edge<-result$edge-1
UPGMA_edgelength<-result$edge.length
UPGMA<-data.frame(UPGMA_edge,round(UPGMA_edgelength,digits=2))
colnames(UPGMA)<-c("to","from","weight")
N<-dim(UPGMA)[1]
for (i in 1:N)
{
  UPGMA<-rbind(UPGMA,c(UPGMA[i,]$from,UPGMA[i,]$to,UPGMA[i,]$weight))
}
result<-UPGMA[order(UPGMA$to),]
for (i in 1:dim(result)[1])
{
  cat(paste(result[i,][1],"->",result[i,][2],":",result[i,][3],sep=""),file="outfile.txt",append=TRUE)
  cat("\n",file="outfile.txt",append=TRUE)
}
#################################
#################################
##########Neighbor Joining#######
#################################
#################################
library(phangorn)  
setwd("C:/Users/Nan Sun/Desktop")
distancematrix<-read.table("al.txt",fill=TRUE)
colnames(distancematrix)<-seq(0,dim(distancematrix)[2]-1,1)
rownames(distancematrix)<-seq(0,dim(distancematrix)[2]-1,1)
result<-NJ(as.matrix(distancematrix))
NJ_edge<-result$edge-1
NJ_edgelength<-result$edge.length
NeiborJoining<-data.frame(NJ_edge,round(NJ_edgelength,digits=2))
colnames(NeiborJoining)<-c("to","from","weight")
N<-dim(NeiborJoining)[1]
for (i in 1:N)
{
  NeiborJoining<-rbind(NeiborJoining,c(NeiborJoining[i,]$from,NeiborJoining[i,]$to,NeiborJoining[i,]$weight))
}
result<-NeiborJoining[order(NeiborJoining$to),]
for (i in 1:dim(result)[1])
{
  cat(paste(result[i,][1],"->",result[i,][2],":",result[i,][3],sep=""),file="outfile.txt",append=TRUE)
  cat("\n",file="outfile.txt",append=TRUE)
}
###################################
###################################

