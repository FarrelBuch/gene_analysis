# Select samples for Whole Exome Sequencing
# multiplex families
# affecteds
# most severe clinical course
# earliest onset
# good DNA (ie whole blood, good ABS ratio, non-extreme concentration)
# have HumanOmni1-Quad data
# part of trios favored over duos
# HPV typing available favored over not



library(ggplot2)


source(file="Multiplex family germline.r")# brings in nuc.ac.multiplex which is nuc.ac.nr of all the those in the multiplex family

mendel.error.fam <- c("BJW18006", "DET10003", "FJB01117", "FKK24002", "GHP04013", "JEM04008", "JWT08002", "PXC08007")#families that were deleted from the TDTae alogrithm "delete 8 families with highest Mendel errors from the dataset", https://mail.google.com/mail/u/0/#apps/subject%3A(Some+additional+details+on+the+plink+run+)/n25/1358bc48a0e710bb

express.rank <- affected[,list(PtCode, dxage, rankage = rankage <- frank(dxage), aggr.max.freq, rank.max.freq = rank.max.freq <- rank(-aggr.max.freq), avg.annual.frq, rank.annual.freq = rank.annual.freq <- rank(-avg.annual.frq), aggr.distal, rank.distal = rank.distal <-  rank(aggr.distal), aggr.tracheostomy, rank.tracheost = rank.tracheost <-  rank(aggr.tracheostomy),aggr.surgcount, rank.surgcount = rank.surgcount <- rank(-aggr.surgcount),hpv, rank.express=(aggr.max.freq*(1+avg.annual.frq)*(1+aggr.surgcount)*rank.tracheost*rank.distal)/(dxage))]#rank.express provides a rank (higher is more expressive) of the expressiviity of any susceptibility. They are people who were diagnosed the youngest and had the most aggressive clinical course. Added a 1 to avg.annual.frq and aggr.surgcount when calculating rank.express since some zeros messed with the multiplication.
express.rank[PtCode=="RYS16028PT", rank.express:=NA]
express.rank[order(-rank.express)][seq(from=1, to=0.5*.N, by = 10 )]
express.rank[is.na(rank.express), rank.express:=signif(quantile(express.rank$rank.express, na.rm = TRUE, probs = 0.1), digits = 2)]# some ranks are NA because there was just one piece of missing data. I therefore assigned all of them to have the 10th percentile rank.exome score
express.rank[rank.express==0]





DNA.rank.WB <- all.subj.byspec[TissCode=="WB",list(PtCode, relationship, nuc.ac.nr,dnaconcngpmicl, absratio, absratiolow, rank.absratio=rank.absratio <- rank(absratiolow), dnayieldmicrog, yieldlow =  yieldlow <- ifelse(-scale(log(dnayieldmicrog))<1,1,-scale(log(dnayieldmicrog))), dna.rank=1/(yieldlow*(0.1+absratiolow)), TissCode)] #DNA.rank provides a rank (higher number is more favorable) of our DNA stock solutions. We chose the dna stock yield that was within 1 sd of the log dnayield or above.
all.subj.byspec[TissCode=="WB",list(dnayieldmicrog, ifelse(-scale(log(dnayieldmicrog))<1,1,-scale(log(dnayieldmicrog))))]
DNA.rank.WB[order(-dna.rank)][seq(from=1, to=.N, by = 15 )]

DNA.rank.MW <- all.subj.byspec[TissCode=="MW",list(PtCode, relationship, nuc.ac.nr,dnaconcngpmicl, absratio, absratiolow, rank.absratio=rank.absratio <- rank(absratiolow), dnayieldmicrog, dna.rank=dnayieldmicrog/(0.1+absratiolow), TissCode)] #DNA.rank provides a rank (higher number is more favorable) of our DNA stock solutions. For MW it is calculated differently. Here we do not favor the average yield but rather the hihgest yield.
DNA.rank.MW[dnayieldmicrog<=0, dna.rank:=0.0001]
DNA.rank.MW[order(-dna.rank)][seq(from=1, to=.N, by = 15 )]

#did dnayield from mw improveover time
ggplot(all.subj.byspec[TissCode=="MW"], aes(x=DNAextract.date, y=log2(dnayieldmicrog))) + geom_point()  + stat_smooth(method = "lm")
summary(lm(log2(dnayieldmicrog)~DNAextract.date, data =all.subj.byspec[TissCode=="MW"]))
# library(BayesianFirstAid)
# fit2 <- bayes.t.test(log2(dnayieldmicrog)~DNAextract.date>as.IDate("2009-01-01"), data=all.subj.byspec[TissCode=="MW"])
# summary(fit2)
# plot(fit2)

t.test(log2(dnayieldmicrog)~DNAextract.date>as.IDate("2009-01-01"), data=all.subj.byspec[TissCode=="MW"])
ggplot(data=all.subj.byspec[TissCode=="MW"], aes(x=DNAextract.date>as.IDate("2009-01-01"), y=log2(dnayieldmicrog))) + geom_boxplot()

mw.stock <- subset(all.subj.byspec, TissCode=="MW", c(dnayieldmicrog,DNAextract.date)  )
mw.stock[, extract.era:=ifelse(DNAextract.date>as.IDate("2009-01-01"), "since2009", "before2009")]
ggplot(data=mw.stock, aes(x=log2(dnayieldmicrog))) + geom_histogram(binwidth = 0.8) + facet_grid(extract.era~.) + aes(y = ..density..)





l.all.germ <- list(DNA.rank.WB, DNA.rank.MW)
DNA.rank <- rbindlist(l.all.germ,use.names=TRUE, fill=TRUE)
DNA.rank[, dna.rank:=dna.rank*(rank(TissCode)^2)]#dna.rank therefore accounts for best ratio and good yield and the type of tissue that it was extracted from where WB is better than mouthwash
DNA.rank[order(-dna.rank)][seq(from=1, to=.N, by = 25 )]

ggplot(all.subj.byspec[TissCode!="BUC"], aes(x=log10(dnayieldmicrog))) + geom_histogram(binwidth=0.25) + facet_grid(TissCode~.) + aes(y = ..density..) + ggtitle("Yield of DNA from each tissue type")

illumina.omni <- rpinfinwrk.dt[day2blue=="Y",nuc.ac.nr]# nuc.ac.nr of all specimens that we have illumina.omni data on

setkey(DNA.rank, PtCode)
setkey(express.rank, PtCode)


ordered.exome <- merge(express.rank, DNA.rank[relationship=="PT",.SD, .SDcols=-relationship])#do not want the DNA specimens from relations getting in  here therefore filter for affected patients only
ordered.exome[,c("rankage", "rank.max.freq", "rank.annual.freq", "rank.distal", "rank.tracheost", "rank.surgcount", "absratiolow", "rank.absratio"):= NULL]

# get a field in for illumina genotyping
ordered.exome[, illumina:="noillumina"]
ordered.exome[nuc.ac.nr %in% illumina.omni, illumina:="haveillumina"]
ordered.exome[, rank.illumina:=rank(illumina)]#lower number is better so it should go in the denominator

# get a field in for trio vs duo vs solo
fam.dt <- all.subj.byspec[,list(father=any(relationship=="FA"), mother=any(relationship=="MO")),by=family]# fam.dt is a data.table to put the family status in one table
fam.dt[, PtCode:=paste0(family,"PT")]
fam.dt[, trio.duo.sol:=factor(1+father+mother, labels = c("solo", "duo", "trio"))]
setkey(fam.dt,PtCode)
ordered.exome[fam.dt, trio.duo.sol:= i.trio.duo.sol]# read about this at data.table join then add columns to existing data.frame without re-copy at http://stackoverflow.com/questions/19553005/data-table-join-then-add-columns-to-existing-data-frame-without-re-copy
# add a score for HPV type to multiply to numerator to find best affecteds to send for whole exome sequencing.
ordered.exome[,rank.hpv:=1]
ordered.exome[!is.na(hpv), rank.hpv:=2]
ordered.exome[hpv==6|hpv==11, rank.hpv:=4]


ordered.exome[,rank.exome := dna.rank^2*rank.express*as.numeric(trio.duo.sol)*rank.hpv/rank.illumina^2]


ordered.exome[is.na(rank.exome), rank.exome:=signif(quantile(ordered.exome[,rank.exome], na.rm = TRUE, probs = 0.1), digits = 2)]# some ranks are NA because there was just one piece of missing data. I therefore assigned all of them to have the 10th percentile rank.exome score

# Generate list of suitable samples-------------------------
exome.list <- ordered.exome[!nuc.ac.nr %in% nuc.ac.multiplex][sample(.N,22, replace=FALSE, prob=rank.exome)][order(-rank.exome)]#gives us random sample weighted for most suitable in descending order of suitability
nuc.ac.exome.list <- exome.list[,nuc.ac.nr]

#how about some adult onset since it may be a different disease
exome.ao <- ordered.exome[dxage>18][sample(.N,8, replace=FALSE, prob=(rank.exome))][order(-rank.exome)]
nuc.ac.exome.ao <- exome.ao[,nuc.ac.nr]


#how about some sets of parents from trios but not with mendelian errors
# consulted geneticist will not do trio parents at this time Wednesday, 15 Apr 2015 15:01
exome.trio <- exome.list[!PtCode %chin% paste0(mendel.error.fam,"PT")&trio.duo.sol=="trio"&TissCode=="WB"][sample(.N, 0, replace=FALSE, prob=rank.express),PtCode]
exome.trio <- str_sub(exome.trio,end=8L)#to convert the names from PtCode to trio name
exome.parents <- DNA.rank.WB[str_sub(PtCode, end=8L) %chin% exome.trio&relationship %chin% c("FA", "MO")]

DNA.rank.WB[relationship %chin% c("FA", "MO")]

nuc.ac.exome.parents <- exome.parents[,nuc.ac.nr]

#how about some children with mild disease
# instead of weighting by rank.exome, we will weight by rank.exome divided by square of expressivitiy
# set a limit that their age of diagnosis must be <12
exome.indolent <- ordered.exome[!nuc.ac.nr %in% nuc.ac.multiplex&dxage<12&aggr.max.freq<4&aggr.surgcount<10&aggr.distal=="cleardist"&aggr.tracheostomy=="Never"&!(is.na(aggr.surgcount)|is.na(aggr.max.freq)|is.na(avg.annual.frq)|is.na(aggr.distal)|is.na(aggr.tracheostomy)|is.na(dxage))][sample(.N,8, replace=FALSE, prob=rank.exome/rank.express^2)][order(-rank.exome)]#gives us random sample weighted for most suitable in descending order of suitability, Eliminated rank
summary(ordered.exome$rank.express)
nuc.ac.exome.indolent <- exome.indolent[,nuc.ac.nr]

nuc.ac.for.exome <- unique(c(nuc.ac.multiplex, nuc.ac.exome.list, nuc.ac.exome.ao, nuc.ac.exome.parents, nuc.ac.exome.indolent))
length(nuc.ac.for.exome)

#save(nuc.ac.for.exome, file = "nucleic acid whole exome.RData")


multip.n <- length(nuc.ac.multiplex)
high.penetrance.n <- length(nuc.ac.exome.list)
ao.n <- length(nuc.ac.exome.ao)
parents <- length(nuc.ac.exome.parents)
indolent <- length(nuc.ac.exome.indolent)
sum(multip.n, high.penetrance.n, ao.n, parents, indolent)

# find the samples for evan----------------
load(file="plating data minus 106.RData")
setkey(plating,nuc.ac.nr)
gofind <- plating[nuc.ac.nr %in% nuc.ac.for.exome,list(nuc.ac.nr,plate,well)][order(plate,well)]#samples in the big agena send out

gofind.106 <- setdiff(nuc.ac.for.exome, plating$nuc.ac.nr)#samples that were sent out in the initial 106
go.find.all.exome <- rbindlist(list(gofind, data.table(nuc.ac.nr=gofind.106)), fill=TRUE, use.names=TRUE)
write.csv(go.find.all.exome, file = "go find all exome.csv", row.names = FALSE)

# problem specimens notified 2015-04-22--------------
# http://gsl.hudsonalpha.org/projects/haib15FJB3120/view_project

problem.sampl <- data.table(hdslph=c(2, 1, 30, 17, 8, 6), nuc.ac.nr=c(760, 736, 600, 793, 759, 384), problem=c("low quantity", "low quantity", "degraded", "degraded", "low quantity", "degraded"), key="nuc.ac.nr")
problem.but.multiplex <- intersect(problem.sampl$nuc.ac.nr, nuc.ac.multiplex)
problem.sampl[J(problem.but.multiplex), resolve:="continue because valuable multiplex"]

setkey(ordered.exome,nuc.ac.nr)
ordered.exome[problem.sampl]# appears as if nuc.ac.nr 600 was selected because it came from a an affected person with a high expressivity of the disease and their absorbance ratio and yield of DNA was good notwithstanding that it was from a mouthwash specimen. Let us replace it

problem.sampl[J(c(600, 384)),resolve:="replace with another sample from high expressivity"]
load("nucleic acid whole exome.RData")# do this so we get the actual original specimens sent

exome.list <- ordered.exome[!nuc.ac.nr %in% c(nuc.ac.multiplex, nuc.ac.for.exome) ][sample(.N,1, replace=FALSE, prob=rank.exome)][order(-rank.exome)]#gives us random sample weighted for most suitable in descending order of suitability
#nuc.ac.exome.substitute <- c(exome.list[,nuc.ac.nr], nuc.ac.exome.substitute)
#save(nuc.ac.exome.substitute, file = "nucleic acid whole exome substitute.RData")



