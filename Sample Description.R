# This script is used to provide descriptive statistics of the sample that goes into a paper

# Run Munge RRP.r first

source("Munge RRP.R")
# eg. we will use affected to describe the disease and all.subj to describe the genders and the ethnicities and the races and the family compositions.

# need all.subj and affected

# first we select based on which sample we are describing

require(ggplot2)
require(reshape2)

# select individuals genotyped for 296 illumina----------------
# they are in PLINK.pedigree
# some subjects were excluded from the analysis

subj.296 <- PLINK.pedigree[,PtCode]
subj.296 <- subj.296[-grep(pattern="^JLJ09006|^BJW18006|^DET10003|^FJB01117|^FKK24002|^GHP04013|^JEM04008|^JWT08002|^PXC08007", x=subj.296)] #Delete two members with pedigree ID JLJ09006 as there is no child in this family, delete 

affec.296 <- affected[PtCode %chin% subj.296,]
all.sub.296 <- all.subj[PtCode %chin% subj.296,]



# decide which sample will be analyzed--------------
# affec.samp will be only the patients (affected individuals) who we are analyzing
# all.sub.samp will be the all the subjects in the group to be analyzed
# affec.samp <- affec.296
# all.sub.samp <- all.sub.296
affec.samp <- affected
all.sub.samp <- all.subj

# Numbers------------------------

# how many subjects
n <- affec.samp[,.N]

# how many families
n.fams <- length(unique(all.sub.samp[,family]))

# family composition
fam.composition <- all.sub.samp[,list(.N, relationships=str_c(relationship,collapse=", ")),by=family]
fam.composition[,list(n=nr <- .N, perc=round(nr/n.fams,2)) , by=relationships]
all.sub.samp[,.N,by=relationship]



# gender of affected, race of affected, ethnicity of affected
gend <- affec.samp[,.N,by=Sex]
gend[,list(Sex,N,perc=round(N/sum(N),2))]
race <- affec.samp[,.N,by=Race]
race[,list(Race,N,perc=round(N/sum(N),2))]
eth <- affec.samp[,.N,by=Ethnicity]
eth[,list(Ethnicity,N,perc=round(N/sum(N),2))][order(-N)]

# what age of diagnosis, what age now
# age at enrollment
summary(affec.samp[,list(enrolage, dxage, follup.yr)])
ggplot(affec.samp, aes(x=enrolage)) + geom_histogram(binwidth=5)
ggplot(affec.samp, aes(x=enrolage)) + stat_ecdf()
ecdf.enrolage <- ecdf(affec.samp[,enrolage])
ages <- sort(c(seq(from=5,to=30,by=5),3, 18,50))
round(ecdf.enrolage(ages),2)

# age at diagnosis
ggplot(affec.samp, aes(x=dxage)) + geom_histogram(binwidth=5)
ggplot(affec.samp, aes(x=dxage)) + stat_ecdf()
ecdf.dxage <- ecdf(affec.samp[,dxage])
round(ecdf.dxage(ages),2)

# duration of followup (diagnosis to last contact for which we have data reporting)
ggplot(affec.samp, aes(x=follup.yr)) + geom_histogram(binwidth=5)
ggplot(affec.samp, aes(x=follup.yr)) + stat_ecdf()
quantile(affec.samp[,follup.yr],c(0.20, 0.25, 0.5, 0.75, 0.80), na.rm=TRUE)
ecdf.follup <- ecdf(affec.samp[,follup.yr])
round(ecdf.follup(c(0.9, 1:10)),2)

  # makes a table of ages and what cumulative percentage 
age.percentile <- data.table(ages, perc.dx=round(ecdf.dxage(ages),2), perc.enrl=round(ecdf.enrolage(ages),2))
age.percentile
write.csv(age.percentile, file="age percentile sample.csv", row.names=FALSE)

# Geography
country <- affec.samp[,.N, by=country]
country[,list(country, N, perc = round(N/sum(N),2))]


# HPV type
hpv <- affec.samp[,.N,by=hpv][order(hpv)]
hpv[,list(hpv, N, perc = round(N/sum(N),3))]
hpv[hpv==6|hpv==11,list(hpv, N, perc = round(N/sum(N),3))]# limiting the distribution of hpv type to only those that are a known 6 or a known 11, excluding those that are unknown or those that are 6 and 11

affec.samp[,.N,by=hpv][,prop := N/sum(N)]

hpv.onset <- affec.samp[,.N,by=list(hpv,onset)][order(onset,hpv)]
hpv.onset.wide <- dcast(hpv.onset, formula=onset~hpv, value.var = "N")
hpv.onset
hpv.onset.wide


# Agggressive or not, component of aggressive

summary(affec.samp[, list(trach=factor(aggr.tracheostomy), distal=factor(aggr.distal))])
summary(affec.samp[, list(aggr.surgcount, surg.count.10 = factor(ifelse(aggr.surgcount>=10,"aggr","indolent")), aggr.max.freq, surg.freq = factor(ifelse(aggr.max.freq>=4,"aggr","indolent")), SurgIntervDay) ])

ggplot(affec.samp, aes(aggr.tracheostomy)) + geom_bar(width=0.5) + scale_y_continuous("number of affected subjects", limits=c(0,90))
ggplot(affec.samp, aes(aggr.distal)) + geom_bar(width=0.5) + scale_y_continuous("number of affected subjects", limits=c(0,90))

theme_set(theme_bw())
affec.samp[,.N, by=list(aggr.tracheostomy,aggr.distal )]
aggressive.matrix <- affec.samp[, list(PtCode, aggr.tracheostomy = aggr.tracheostomy=="trached", aggr.distal = aggr.distal=="involvdist", aggr.surgcount = aggr.surgcount>=10, aggr.max.freq = aggr.max.freq>=4, aggr.composite)]
aggr.true  <-  aggressive.matrix[,list(aggr.tracheostomy = sum(aggr.tracheostomy, na.rm=TRUE), aggr.distal = sum(aggr.distal, na.rm=TRUE), aggr.surgcount = sum(aggr.surgcount, na.rm=TRUE), aggr.max.freq = sum(aggr.max.freq, na.rm=TRUE), aggr.composite = sum(aggr.composite, na.rm=TRUE))]
aggr.prop  <-  aggr.true/aggressive.matrix[,.N]
rbindlist(list(round(aggr.prop,2), aggr.true))
setnames(aggressive.matrix, old=c("aggr.tracheostomy", "aggr.distal", "aggr.surgcount", "aggr.max.freq", "aggr.composite"), new = c("trach", "distal", "surg", "freq", "composite"))

aggr.compostion.count <- aggressive.matrix[,.N,by=list(trach, distal, surg, freq)][order(trach, distal, surg, freq)]
aggr.compostion.count

#barplot raw metrics of aggressiveness
melted.aggr<-melt(data=aggressive.matrix,id.vars="PtCode", value.name="aggressive") 

ggplot(melted.aggr, aes(x=variable,fill=aggressive, order = -as.numeric(aggressive))) + geom_bar()  + labs(title="Aggressiveness of clinical course",x = "clinical characteristic", y = "count of subjects") + theme_gray(22)

ggsave(filename="aggresivness of clinical course.svg", width=88, height=88, units="mm", scale=2)

# any malignant transformation------------------
names(affec.samp)
ldt <- affec.samp[, lapply(.SD, function(x) grepl("malig|carci|transf|cancer|dysplas", x))] # searches for text string in any columng and marks true or false which is the same as 1 or 0
affec.samp[rowSums(ldt)>0, list(PtCode, SitePulmonay, SiteBronchus, SiteTrachea, CoComments, MoComments, FaComments, CareComments, Comments)] # lists the rows in which there was anything true

#ftable of aggr criteria and composite

