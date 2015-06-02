# last HPV linear array
library(data.table)
# make a mix of HPV 6 and 11
#mixes <- sort(sample(1:9, size=4, replace=FALSE))
#save(mixes, file="HPVmixes.RData")
data.table(HPV6=mixes/10, HPV11=1-mixes/10)
