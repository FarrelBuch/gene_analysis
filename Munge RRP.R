# This script relies on "read genomics in.r" to read from GENOMICS and Google. At that stage it is either ready for use by this script or it is ready for encryption.
source(file = "read genomics in.r")
#load(file="V:/GRRPsnapshot.RData")#if I need to load from a encrypted file
#load(file="D:/Users/Farrel/My Documents/RRPGenSucep/GRRPsnapshot.RData")#if I need to load office computer
require(data.table)
require(stringr)
require(zoo)
require(reshape2)

#todo Something is missing from the first row or a row of specimens that has NA in it and it is messing up the IDateTime transformation


# create data tables with names and keys------------------

rrplog <- data.table(rrplog)
clin.upd <- data.table(clin.upd)#read in clinical update
psga <- data.table(psga)
collab <- data.table(collab)
hpv.early <- data.table(hpv.early)
hpv.provtous <- data.table(hpv.provtous)
specimens <- data.table(specimens)
rrp.related <- data.table(rrp.related)
rpdna <- data.table(rpdna)
setnames(x=rpdna,old="rpdnanr",new="nuc.ac.nr")
PLINK.pedigree <- data.table(PLINK.pedigree)
setnames(x=PLINK.pedigree,old=c("V1","SubjLastName"),new=c("nuc.ac.nr","PtCode"))
rpinfinwrk.dt <- data.table(rpinfinwrk)
remove(rpinfinwrk)
setnames(x=rpinfinwrk.dt, old=c("Sample ID", "rpdnanr"), new=c("sample.id", "nuc.ac.nr"))
rpinfinwrk.dt[, ':='(sample.id = paste0(chipnr, "_", region),
										 assaystarted=as.IDate(assaystarted, origin="1899-12-30"))]

sapply(mget(ls()[sapply(mget(ls(), .GlobalEnv),is.data.table)], .GlobalEnv),dim)


setkey(rrplog,PtCode)
setkey(clin.upd,PtCode)
setkey(psga,FJBcode)
setkey(collab,CollabCode)
setkey(specimens,SubjLastName)
setkey(rrp.related,PtCode)
setkey(hpv.early, pt.code)
setkey(hpv.provtous, pt.code)
setkey(rpdna,nuc.ac.nr)
setkey(PLINK.pedigree,PtCode)
setkey(rpinfinwrk.dt, nuc.ac.nr)
setkey(approved.sat,InstCode)


# There are some horrid naming conventions HS# for human subject specimen number and NA# for nucleic acid number are serriously messing up R script
setnames(specimens,old=c("HS#","NA#"),new=c("hum.sub.nr","nuc.ac.nr"))
hpv.early[,Kbddate := as.IDate(Kbddate, origin="1899-12-30")]# have done this date conversion earlier than in other tables becaues may need it earlier.


# correct typographical error--------------------

# HS#37292 appears to be incorrectly entred as RYS16002PT yet in the paper "bible" it is RYS16022PT, typo that must be fixed
specimens[hum.sub.nr==37292, SubjLastName:="RYS16022PT"]

# The FJBCode in psga whould all have the suffix of "PT"
psga[,FJBcode := paste0(FJBcode,"PT")]
setkey(psga, FJBcode)

# RP1392 and RP1393 are both whole blood from RYS16015PT received 2010-11-15 but suspect 1393 is LRNXP #TODO confirm this is true
specimens[nuc.ac.nr==1393, TissCode:="LRNXP"]


# clean up collaborator's table to only include those relevant to us---------------
collab.codes <- sort(unique(c(rrplog[,CollabCode],rrp.related[,CollabCode])))
collab <- collab[collab.codes]
remove(collab.codes)
collab[is.na(Country),]
# Vanderbilt got put in twice as an insitution,
# deal with it in the collaborators table
collab["CTW",InstCode:="VUMC"]

# deal with in the rrplog
rrplog[InstCode=="CHV", ':='(InstCode="VUMC",
                             InstName="Vanderbilt University Medical Center")]
unique(rrplog[grep("V", InstCode), list(InstCode, InstName)])


# add city and country to each rrplog row------------------------

setkey(rrplog,CollabCode)
setkey(collab,Country)
collab["United States of America",Country:="USA"]
setkey(collab,CollabCode)
rrplog[collab, country:=Country] # adds a column to the rrplog with the hpv type, lookup by assignment
rrplog[collab, c("city", "zip") := list(City, ZipCode) ] #
rrplog[CollabCode=="FB3",list(PtCode,country,zip)] #multiple assignment in one line
setkey(rrplog, PtCode)
# FJB collab were acquired through patient support group and we must look to psga for location
psga[grep(pattern="United States of Ame", x=PTCntry, ignore.case=TRUE), PTCntry := "USA"]
rrplog[psga, c("country", "city", "zip") := list(PTCntry, PTCity, PTZip) ] # multiple lookup and assignment in one line


# merge HPV type data-------------------------

# TODO merge the different sources of HPV typing data
# as of Saturday, 26 Oct 2013 14:41 there is
# 1) Joseph Donfack hpv.early, done in the mid 00's
# 2) hpv.provtous (hpv typing that was provided to us)
# 3) calls.by.ABI.cluster
# Start by merging SubjLastName into calls.by.ABI.cluster by using specimens table, join on RP nucleic acid number
# Create a union list of all HPV typing
# Look to see if we have any results done twice, then normalize
# merge the HPV type into the rrplog where it can then flow through the data
# 2013-10-26 14:46

ABI.cluster.HPV <- merge(x=calls.by.ABI.cluster,y=specimens[,list(nucleic.acid.nr=nuc.ac.nr, SubjLastName)], all.x=TRUE, all.y=FALSE)
names(ABI.cluster.HPV)
setnames(ABI.cluster.HPV, old="SubjLastName", new="pt.code")
ABI.cluster.HPV[,source := "sequence clust"]
hpv.type <- rbind(ABI.cluster.HPV,hpv.early[,list(nucleic.acid.nr,hpvtype,pt.code, source= "aspcr rflp")])
hpv.type <- rbind(hpv.type,hpv.provtous[,list(nucleic.acid.nr=NA , hpvtype, pt.code, source = "provided to us")])# hpv.type is union list of all HPV types
setnames(hpv.type, old="pt.code", new="PtCode")
setkey(hpv.type, PtCode)
# For one or other reason we have some duplicated typing (good for checking)
hpv.retyped <- hpv.type[duplicated(hpv.type[,PtCode])|duplicated(hpv.type[,PtCode], fromLast=TRUE),][order(nucleic.acid.nr)]#patients who have had hpv genotyping twice
save(list=c("ABI.cluster.HPV", "hpv.early","hpv.provtous", "hpv.type", "hpv.retyped" ), file="hpv typing tables.Rdata")
# lets evaluate specimens that have been typed a second time
hpv.retyped.wide  <-  dcast(data=hpv.retyped, formula=PtCode + nucleic.acid.nr ~ source, value.var="hpvtype")
write.csv(hpv.retyped.wide, file = "hpv retypted.csv", row.names = FALSE)
# remove duplicates from hpv.type to normalize
hpv.type[,source := factor(source, levels=c("aspcr rflp", "sequence clust", "provided to us"))]
hpv.unique <- unique(hpv.type)
rrplog[hpv.unique, hpv:=hpvtype]# adds a column to the rrplog with the hpv type, lookup by assignment
rrplog[,hpv:=factor(hpv, levels=c( "6", "11", "16", "611", "retest"))]
write.csv(hpv.unique[,list(PtCode,hpvtype,source)], file = "hpv type share.csv", row.names = FALSE)


# get rid of missing values-------------
rrplog[MoEduPast=="Not Recorded",MoEduPast:=NA_character_]
rrplog[MoEduBef=="Not Recorded",MoEduBef:=NA_character_]
rrplog[AOPtEdu=="Not Recorded",AOPtEdu:=NA_character_]
rrplog[,MoEduPast:=factor(MoEduPast,levels=c("Grade9","Grade12","HighSchool","CollegeSome","AssociateDeg","Bachelor","Graduate"))]
rrplog[,MoEduBef:=factor(MoEduBef,levels=c("Grade9","Grade12","HighSchool","CollegeSome","AssociateDeg","Bachelor","Graduate"))]
rrplog[,AOPtEdu:=factor(AOPtEdu,levels=c("Grade9","Grade12","HighSchool","CollegeSome","AssociateDeg","Bachelor","Graduate"))]
rrplog[,MaxEdu:=factor(pmax(as.numeric(AOPtEdu),as.numeric(MoEduBef),as.numeric(MoEduPast),na.rm=TRUE),labels=c("Grade9","Grade12","HighSchool","CollegeSome","AssociateDeg","Bachelor","Graduate"))]
rrplog[grepl("999",IncomeYBef),IncomeYBef:=NA]
rrplog[grepl("999",IncomeYPast),IncomeYPast:=NA]
rrplog[grepl("999",AOPtIncomeY),AOPtIncomeY:=NA]
rrplog[grepl("9999",YearDx),YearDx:=NA]
rrplog[,max.incom:=pmax(IncomeYBef, IncomeYPast, AOPtIncomeY,na.rm=TRUE)]
#  TODO At one stage it appears as if any income of less than $1000 in the USA
#  was declared NA. At this stage I do not have countries in the rrplog
#  take care of this later if it does indeed need to be taken care of
#  rrplog[max.incom<1000&Country=="USA", max.incom:=NA]
rrplog[grepl("999",PeopleNrBef),PeopleNrBef:=NA]
rrplog[grepl("999",PeopleNrNow),PeopleNrNow:=NA]
rrplog[grepl("99",MonthDx),MonthDx:=NA]
rrplog[grepl("99",MoDOBYear),MoDOBYear:=NA]
rrplog[grepl("99",MoDOBMon),MoDOBMon:=NA]
rrplog[grepl("99",StageScore),StageScore:=NA]
rrplog[StageScore==0,StageScore:=-1L]#since zero will kill the transformation
rrplog[,SurgIntervType:=factor(SurgIntervType, levels=c("pre-determined intervals","hybrid intervals","prn intervals"))]
rrplog[grepl("99",SurgNumLife),SurgNumLife:=NA]
rrplog[grepl("99",SurgIntervDay),SurgIntervDay:=NA]
rrplog[SurgIntervDay==0,SurgIntervDay:=NA] # some affecteds were enrolled at the time of the first surgery and therefore there was no previous surgery. They may have been makred as 0 but in fact they should be missing
rrplog[grepl("99",SurgNum12Mon),SurgNum12Mon:=NA]
rrplog[grepl("99",MoSexAge),MoSexAge:=NA]
clin.upd[grepl("99",NrPrecTwelM),NrPrecTwelM:=NA]
clin.upd[grepl("99",UpdSurgNumLife),UpdSurgNumLife:=NA]
rrplog[grepl("99",SurgMaxNum12Mon),SurgMaxNum12Mon:=NA]



# transform all date fields in to IDate-----------------
# some months of diagnosis are missing. Make it 6, this may lead to trouble but lets see
rrplog[is.na(MonthDx),MonthDx:=6L]
rrplog[,dxdate:=paste(YearDx, MonthDx, "15",sep="-")]
rrplog[grepl("NA",dxdate),dxdate:=as.character(NA)]
rrplog[,dxdate:=as.IDate(dxdate)]

rrplog[is.na(MoDOBMon),MoDOBMon:=6L]
rrplog[,moth.dob:=paste(MoDOBYear, MoDOBMon, "15",sep="-")]
rrplog[grepl("NA",moth.dob),moth.dob:=as.character(NA)]
rrplog[,moth.dob:=as.IDate(moth.dob)]
rrplog[,moth.dob]

date.cols.rrp <- c("SpecDate", "DOB", "PostEnrollment12MonthUpdateFormMailedOn")
for(date.col in date.cols.rrp) {rrplog[,date.col:=as.IDate(rrplog[[date.col]]),with=FALSE]}# tip: this is how to perform a function on multiple columns of a data.table at the same time. It is quicker than lapply #tip
clin.upd[,ClinUpdDate:=as.IDate(ClinUpdDate)]
clin.upd[,c("kbddateclin.upd","kbdtimeclin.upd"):=IDateTime(KbdDate)]
specimens[,c("kbddatespecim","kbdtimespecim"):=IDateTime(KbdDate)]
rpdna[,specdate:=as.IDate(specdate, origin="1899-12-30")]
rpdna[,extractdate:=as.IDate(extractdate, origin="1899-12-30")]
setnames(x=rpdna,old=c("specdate","extractdate"),new=c("spectophotom.date","DNAextract.date"))

#now that we have transformed date fields, correct typographical errors
rrplog[PtCode=="RYS16028PT", DOB:=as.IDate("1974-08-28")]#error corrected by email from Riaz Seedat subject:Please check DOB and year of diagnosis on RYS16028PT https://mail.google.com/mail/u/0/#inbox/14cbd9b96bd2f34d


# generate some useful fields------------------
# quality DNA has a A260/A280 of 1.8, low, but not high, ratio indicative of a problem
rpdna[absratio>1.8, absratiolow := 0]
rpdna[absratio<=1.8, absratiolow := 1.8-absratio]
# need to create a column in specimens that specifies the yield of DNA in µg
rpdna[,dnayieldmicrog:=dnaconcngpmicl*volumemicl*0.001]


affected <- merge(clin.upd[,.SD, .SDcol=-kbdtimeclin.upd],rrplog, all.y = TRUE )
clin.upd[,.SD, .SDcol=-kbdtimeclin.upd]
key(clin.upd)
key(rrplog)

# # remove illegitimates-----------------
####
# If analysis took place before removal then this must be commented out
# otherwise we end up with unexplaned missing values
####

# Some specimens are illegitimate. For intance JED08002 consent was not done correctly and therefore should be expelled. All the enrollees from Brazil failed to get CONEP approval
# All the rows from  specimens table must be removed
# All the rows from affected table must be removed
# JED08002 to be deleted
affected <- affected[!"JED08002PT"]
specimens <- specimens[grep(pattern="^JED08002",x=SubjLastName,invert=TRUE),]
# Brazil CONEP rejected, Patrick Froelich from EHH sent papilloma but not enrolled through study.
reject.instit <- c("HSP","FUSP","EHH")
setkey(affected,InstCode)#index by InstCode
setkey(specimens,InstCode)
affected <- affected[!reject.instit,]
specimens <- specimens[!reject.instit,]
# some specimens are listed in the specimens table but not present in the corresponding nucleic acid table
# I think that they were entered in the HSlog which wrote them to specimens and then when removed from HSlog in Genomics they remained in the specimens table


setkey(affected,PtCode)# put the index back as expecting
setkey(specimens,SubjLastName)

# GJC03001 is duplicate of GJC03005 and therefore GJC03001 must be deleted post hoc
# see email subject:GJC03001 and GJC03005 with Graeme Copley
specimens <- specimens[!SubjLastName %chin% c("GJC03001PT", "GJC03001MO")]
affected <- affected[PtCode!="GJC03001PT"]
specimens[nuc.ac.nr==257]

# neaten the data-------------
# InstCountry has some values as United States of America and also has USA, combine has USA
setkey(specimens,InstCountry)
specimens["United States of America",InstCountry:="USA"]
setkey(specimens,SubjLastName)
# need a column that specifies family
specimens[,family:=str_sub(SubjLastName,end=8L)]
# need a column that specifies type of relationship
specimens[,relationship:=str_sub(SubjLastName,start=9L)]
specimens[grep("^S",relationship),relationship:="SI"]#each sibling had a different suffix, make them all SI
# do the relationship thing in specs.post.crash as well
specs.post.crash[,relationship:=str_sub(SubjLastName,start=9L)]#there were no siblings after 2013-11-21






# Need to get one union list of affected plus relatives with all the affected information therein------------
#cat(shQuote(names(specimens)), sep=", ")# gives me a list of all the column names so that I can use it in statement
same.col <- intersect(names(specimens),names(affected))
col.in.rrplog <- setdiff(names(affected),names(specimens))
# start by removing junk columns from specimens table
specimens[,c("CoInstCode", "StudyCode", "SpeciesCode","TissDesc", "TissRef", "CollectBy", "CollectLab", "PtName",  "MolDiag", "ResNum", "MolOncNum", "PresumDiag", "Freezer", "Box", "Depleted", "SpecQty", "CultQty", "SerQty", "SerLoc", "QLQty", "QLLoc", "Verified", "OldInfo", "File#", "SNP", "SubjFirstName", "SubjMiddleName", "StudyShort", "Slot", "Rack", "DeleteMsg", "BSNum", "subjSSN", "Stu#", "SubjDegree", "SubjPrefix", "SubjSuffix", "PtCode", "CoSuffix", "CoPrefix","CoDegree"):=NULL,with=FALSE]# PtCode in the Specimen Table is useless, instead we rely on SubjLastName
# Somehow the date on which we received FB301004PT was never entered. We believe it is 2010-01-14
specimens["FB301004PT",RecdDate:=as.Date("2010-01-14")]

# GENOMICS Database crashed 2013-11-21, get later specimens in------------------
l = list(specs.post.crash, specimens)
l = list(specimens[,c(1:26, 28, 29), with=FALSE], specs.post.crash)# kbdtimespecim (27th column was messing with results)
specimens <- rbindlist(l, fill = TRUE)







# need a union table of subject to feed numbers to regulators--------------
affec.demograph <- affected[,list(PtCode,Sex, Ethnicity, Race, DOB, SpecDate)]
setkey(affec.demograph,PtCode)
setkey(specimens,SubjLastName)
all.subj.byspec <- affec.demograph[specimens]
setkey(all.subj.byspec,nuc.ac.nr)
all.subj.byspec <- merge(x=all.subj.byspec,y=rpdna,all.x=TRUE)#we need all the rows from all.subj.byspec but we do not need specimens that have been excluded such as Patrick Froelich PAF and Brazil
setkey(all.subj.byspec,PtCode)
all.subj.byspec[,enroll.date:=pmin(SpecDate, CollectDate, RecdDate,kbddatespecim, na.rm =TRUE)]
setkey(all.subj.byspec,PtCode,enroll.date)# we have now sorted by PtCode and enroll.date which could be different from specimen to specimen because of how enroll.date was created, which was by specimen
setkey(all.subj.byspec,PtCode)#now they should remain sorted but the only key is Ptcode, so when we issue unique() command it will take the first row

all.subj <- unique(all.subj.byspec)#we had duplicate entries since 2 specimens from the same subject generated two rows, one for each specimen
all.subj[,c("TissCode", "comments",  "slotnr",  "boxnr",  "freezernr",  "spectophotom.date",  "extractby",  "DNAextract.date",  "volumemicl",  "dilfactor",  "absratio",  "abs280",  "abs260",  "dnaconcngpmicl",  "hsnr",  "hum.sub.nr",  "nuc.ac.nr", "absratiolow"):=NULL]#do not want specification of tissue here in listing of subjects
# use the RRPlog DOB as most accurate, where there is no RRPlog DOB use the SubjDOB from the specimens log, once used, get rid of SubjDOB from Specimens
all.subj[is.na(DOB)&!is.na(SubjDOB),DOB:=SubjDOB]
all.subj[,SubjDOB:=NULL]
all.subj[str_sub(PtCode, start = 9L, end = 9L)=="S"&is.na(DOB)&grepl("[[:digit:]]{6}$",PtCode),DOB:=as.Date(x=paste(str_sub(PtCode,start=10L, end=15L),"15",sep=""),format="%Y%m%d")]#extracts siblings PtCode to create a date of birth






all.subj[,age:=round(as.numeric(pmin(SpecDate, RecdDate,CollectDate,na.rm =TRUE)-DOB)/365.25,1)]
all.subj[,child.adult:=cut(age,breaks=c(-Inf,18.0,+Inf),labels=c("child","adult"))]#classify into adult and child
all.subj[str_sub(PtCode,start=9L) %chin% c("MO", "FA", "MGM", "MGF"), child.adult:="adult"]# many grand parents or parents do not have dob but we know they are adult
all.subj[str_sub(PtCode,end=8L)=="GHP04015",child.adult:="adult"]# complicated family laid out in email GHP04015 missing parent's questionnaire dated 2009-02-17
# for siblings where we do not know the year of birth, we will assume that the sibling is the same category as the affected individual
no.age.classif <- all.subj[is.na(child.adult),PtCode]
all.subj[paste(str_sub(no.age.classif,end=8L),"PT", sep = ""),child.adult]
all.subj[str_sub(PtCode, start = 9L, end = 9L)=="S"&is.na(child.adult),age.hook:=paste(str_sub(PtCode,end=8L),"PT", sep = "")]
setkey(all.subj,PtCode)
all.subj[str_sub(PtCode, start = 9L, end = 9L)=="S"&is.na(child.adult),child.adult:=all.subj[as.character(PtCode) %chin% age.hook,child.adult]]
# for relatives we will assume that the race and ethnicity is the same as the affected individual
# use zoo na.locf to carry known values forwards and backwards in a family
setkey(all.subj,PtCode)
all.subj[,`:=`(Ethnicity=na.locf(Ethnicity,na.rm=FALSE),Race=na.locf(Race,na.rm=FALSE)),by=str_sub(PtCode,end=8L)]
all.subj[,`:=`(Ethnicity=na.locf(Ethnicity,na.rm=FALSE,fromLast=TRUE),Race=na.locf(Race,na.rm=FALSE,fromLast=TRUE)),by=str_sub(PtCode,end=8L)]
# Gender
# SubjSex came from specimen and Sex came from rrplog
# Sex has had greater integrity checking so let that be the default

# Where Sex is NA,copy SubjSex
setkey(all.subj,Sex,SubjSex)
all.subj[is.na(Sex)&SubjSex=="F",Sex:="Female"]
all.subj[is.na(Sex)&SubjSex=="M",Sex:="Male"]
all.subj[,`:=`(SubjSex=NULL,age.hook=NULL)]#no longer need the variables
setkey(all.subj,PtCode)


# Work on parsing affected table rrplog + clin.upd---------------
# previously done in #Amalgamate the Follow Up data with the enrollment data in RRPDataIntegrity
affected[,LastFU:=pmax(ClinUpdDate,SpecDate,na.rm=TRUE)] # LastFU for Last Follow up date

#Tracheotomy

# fix problem of NA in TrachStat. It causes NA to permeate all fields.
# if anniversary update not available then will assign value that was present at enrollment
affected[is.na(TrachStat),TrachStat:= Tracheostomy]
# This problem all seems to stem from GHP04005PT having status as "Used to Have" instead of "Used to have" just a capitalization issue
affected[Tracheostomy=="Never"|TrachStat=="Never",aggr.tracheostomy:="Never"]#Tracheostomy is status at initial enrollment, Trachstat is status at anniversary follow up, aggr. tracheostomy is a variable summing up over the two other variables.
affected[TrachStat=="Currently present"|TrachStat=="Used to have"|Tracheostomy=="Currently present"|Tracheostomy=="Used to have",aggr.tracheostomy:="trached"]
# Total number of surgeries in life so far
affected[,aggr.surgcount:=pmax(UpdSurgNumLife,SurgNumLife,na.rm=T)]#SurgNumLife is the total number of surgeries at the time of enrollment, UpdSurgNumLife is the same number at the anniversary update
# Distal invovlement
affected[,site.distal.enrol:=ifelse(pmax(SitePulmonay, SiteBronchus, SiteTrachea,na.rm=T)==0,"cleardist","involvdist")]
affected[,EverDist:=ifelse(EverDist==0,"cleardist","involvdist")]#todo might need to be fixed
affected[PapBelowVoiceBox=="Not Recorded"|PapBelowVoiceBox=="Unsure",PapBelowVoiceBox:=NA_character_]
affected[PapBelowVoiceBox=="No",PapBelowVoiceBox:="cleardist"]
affected[PapBelowVoiceBox %like% "Yes",PapBelowVoiceBox:="involvdist"]
ftable(affected[,list(site.distal.enrol,
                      EverDist,
                      PapBelowVoiceBox)],exclude = NaN)
affected[site.distal.enrol=="cleardist"|EverDist=="cleardist"|PapBelowVoiceBox=="cleardist",aggr.distal:="cleardist"]
affected[site.distal.enrol=="involvdist"|EverDist=="involvdist"|PapBelowVoiceBox=="involvdist",aggr.distal:="involvdist"]

#Maximum Frequency in any one year period year could come from
#number in the year preceding enrolement - satellite
#number in the 12 months since enrollment - satelliete update
#max number in 12 months is PSGA questionnaire
#calculated average  (as long as follow up has been at least 10 months
affected[,follup.yr:=round(as.numeric(LastFU-dxdate)/365.25,digits=2)]
affected[follup.yr < 0 & follup.yr > -0.044, follup.yr := 0 ] # affected enrolled at the time of diagnosis may have appeared to be last seen (some date in the early part of the month) before their diagnosis because the diagnosis only recrods the month and takes the middle of the month as the date of diagnosis. 0.044 y is equal to 16 days.
affected[,avg.annual.frq:=round(aggr.surgcount/follup.yr,digits=1)]
affected[follup.yr<(9/12),avg.annual.frq:=NA]
affected[,aggr.max.freq:=round(pmax(SurgNum12Mon, NrPrecTwelM, SurgMaxNum12Mon, avg.annual.frq,na.rm=T))]
#minimum of maximum is always going to be 1
affected[aggr.max.freq<1,aggr.max.freq:=1]
affected[,aggr.composite:=aggr.tracheostomy=="trached"|  aggr.distal=="involvdist"|  aggr.surgcount>=10|  aggr.max.freq>=4]

# date and duration calculations--------------
# for example, need to determine age of diagnosis, age of enrollment, age of mother at birth, sexual debut of mother to birth of child
# the date is in all.subj enroll.date
# use http://stackoverflow.com/questions/15005666/seeking-an-better-way-to-add-columns-in-data-table-from-lookup-table
all.enroll <- all.subj[,PtCode,enroll.date]
setkey(all.enroll,PtCode)
affected <- all.enroll[affected]
duration.years <- function(later,earlier) {round(as.numeric(later-earlier)/365.25,2)}
affected[,maternal.debut:=moth.dob+(MoSexAge+0.5)*365.25]#approximate date on which mother lost virginity. They only told us age so I added half year to number to get middle of the year in which they were that age
affected[,dxage := duration.years(dxdate, DOB)]
affected[,enrolage := duration.years(enroll.date, DOB)]
affected[,motherage := duration.years(DOB,moth.dob)]
affected[,virginbirth := duration.years(DOB,maternal.debut)]
affected[,onset:= ifelse(dxage>=18,"adult","child")]
affected[virginbirth< -2|virginbirth>40,list(PtCode,virginbirth)]
class(affected[,virginbirth])

# internal work, find everywhere where TissCode is
# all.data.frames <- ls()[sapply(mget(ls(), .GlobalEnv),is.data.frame)]#this is how we get to manipulate all the data.frames
# for(df in all.data.frames) {print(paste(which(grepl(pattern="TissCode",names(get(df)))),df))}
# grepl(pattern="TissCode",names(get(df))


# Final end products
# affected data table provides listing of affected individuals with all their characterization
# all.subj data.table shows all human subjects in the study with their ethnicity









# Merge hpv type to come later

