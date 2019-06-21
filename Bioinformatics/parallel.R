source("execute.R")
library(snow)
cl<-makeCluster(3,type='SOCK') 
clusterEvalQ(cl,source("execute.R")) 
clusterEvalQ(cl,f<-function(i){RandomizedMotifSearch(Dna,k)})
N<-780
result<-parLapply(cl,1:N,f)
stopCluster(cl) 
temp<-rep(NA,N)
for (i in 1:N)
{
  temp[i]<-sum(ScoreConstruction(result[[i]]))
}
temp<-result[[which.min(temp)]]
sum(ScoreConstruction(temp))
print(cat(temp,sep='\n'))