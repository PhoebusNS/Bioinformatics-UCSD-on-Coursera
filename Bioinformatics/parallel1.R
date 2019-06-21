source("execute1.R")
library(snow)
cl<-makeCluster(3,type='SOCK') 
clusterEvalQ(cl,source("execute1.R")) 
clusterEvalQ(cl,ff<-function(i){GibbsSampler(Dna,k,N)})
M<-42
result<-parLapply(cl,1:M,ff)
stopCluster(cl) 
temp<-rep(NA,M)
for (i in 1:M)
{
  temp[i]<-sum(ScoreConstruction(result[[i]]))
}
temp<-result[[which.min(temp)]]
sum(ScoreConstruction(temp))
print(cat(temp,sep='\n'))