# last HPV linear array
library(data.table)

# make a mix of HPV 6 and 11
#mixes <- sort(sample(1:9, size=4, replace=FALSE))
#save(mixes, file="HPVmixes.RData")
load(file = "HPVmixes.RData")
data.table(HPV6=mixes/10, HPV11=1-mixes/10)

grep("hpv", names(affected), ignore.case = TRUE, value = TRUE)
summary(affected$hpv)


# Strategy

# all laryngeal specimens that have not been typed yet
# random sample of 6 to verify (will not have enough linear array wells)
# random sample of 11 to verify (will not have enough linear array wells)
# all that were neither 6 or 11 but that we have a laryngeal specimen on
# specimens where we got a different result on sequence clustering versus allele specific PCR (will not have enough linear array wells)
# highest priority would be those that were subjected to whole exome or HLA typing

# hpv.type is union list of all HPV types


# neither 6 or 11

hpv.retest <- hpv.type[!hpvtype %in% c(6,11),] #hpv.retest laryngeal specimens that have been done but not clear cut results

# untyped laryngeal specimens
# compare specimens table to hpv.type and select all those laryngeal specimens that have not been typed.

hpv.never.done <- specimens[TissCode=="LRNXP"&!nuc.ac.nr %in% hpv.type$nucleic.acid.nr, list(nuc.ac.nr, SubjLastName, RecdDate)] #hpv.never.done have specimens but have never been typed

# rbindlist of the two above lists

# prioritize by finding out which ones have been HLA typed and whole exome typed.



