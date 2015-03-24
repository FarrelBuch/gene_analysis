# After the data has been read in
# we need to clean it
# Before we can clean it, we need to detect violations of integrity

require(stringr)

source("Munge RRP.R")

# Some specimens from affected subjects have been logged in the HSlog (and thus are in the specimen table) but are not in the RRPlog-----------------

not.in.log <- grep("PT", setdiff(y=rrplog[,PtCode],x=specimens[,SubjLastName]),value=TRUE)
specs.not.in.log <- specimens[not.in.log,list(subject=SubjLastName,date.received=RecdDate,institution=str_sub(InstName,1,25)),mult="first"][order(date.received)]
specs.not.in.log #eventual report table
# it is possible to enter a unique subject code in the clinical update table twice----------------

#list rows with duplicated PtCode
dups.clin.upd <- clin.upd[as.character(clin.upd[which(duplicated(clin.upd)),PtCode]),]
dups.clin.upd #eventual report table

# # find what the second RYS16001PT belongs to
# # it must belong to someone
# # enrolled on or before 2011-08-03
# # enrolled on or after 2009-08-03
# rrplog[,SpecDate:=as.IDate(SpecDate)]
# setkey(rrplog,SpecDate)
# #find out which subjects it could be
# potentials <- rrplog[SpecDate<=as.IDate("2011-08-03")&SpecDate>=as.IDate("2009-08-03"),PtCode]
# #exclude those who already have an entry in the clinical update
# potentials <- grep("RYS",setdiff(potentials,clin.upd[,PtCode]),value=TRUE)
# setkey(rrplog,PtCode)
# rrplog[as.character(potentials),][SurgNumLife<=2]
# remove(potentials)

# Specimens have no date of origin
# How is it possible that we have a spcemin with no specimen date, no received date and no collection date?
nodate <- all.subj.byspec[is.na(SpecDate)&is.na(RecdDate)&is.na(CollectDate),list(PtCode, family, relationship, Sex, Ethnicity, Race, DOB, hum.sub.nr, CollabCode, InstCode, TissCode, CoFirstName, CoLastName, InstName, SubjSex, SubjDOB, SpecDate, RecdDate, CollectDate, kbddatespecim, KbdDate, InKbd, UpdateUser, UpdateDate)]
nodate #eventual report table specimens that seem to arive out of thin air with no date whatsoever

# Need check that every RRP type subjectlastname in specimen table is coded as being in the rp or rrprel study---------------
# start by bringin in the full length of the specimen table but only the columns of study code and subjectslastname
study.check <- data.table(sqlQuery(channellab,"SELECT SubjLastName,StudyCode from CGS.Specimen", believeNRows = FALSE))
setkey(study.check,SubjLastName)
study.check[,rp.pattern:=grepl("^[[:alpha:]]{3}[[:digit:]]{5}[[:alpha:]]|^FB301[[:digit:]]{3}[[:alpha:]]",SubjLastName)]#"FB3" is an exception to the subject code structure and was used for patients I enrolled here at AGH
table(study.check$rp.pattern,study.check$StudyCode)
spec.maybe.rp <- study.check[rp.pattern==TRUE&(StudyCode!="RP"&StudyCode!="RPREL"),]#these are specimens in specimens table that are not in RP or RPREL but maybe should be.
spec.in.rp.wrong.code <- study.check[rp.pattern!=TRUE&(StudyCode=="RP"|StudyCode=="RPREL"),]#these are specimens in specimens table that are in RP or RPREL but do not seem to have the correct code structure
spec.maybe.rp #eventual report table
spec.in.rp.wrong.code #eventual report table

# Check that every mother is female and every father is male--------------------
# List male mothers and female fathers
# for use during script development
# table(specimens$SubjSex,str_sub(specimens$SubjLastName, start = 9L, end = 10L))
trangender.parents <- specimens[(str_sub(SubjLastName, start = 9L, end = 10L)=="FA"&SubjSex=="F")|(str_sub(SubjLastName, start = 9L, end = 10L)=="MO"&SubjSex=="M"),list(SubjLastName,SubjSex,CollabCode,InstCode,InstCity,CollectDate,RecdDate,KbdDate,InKbd,UpdateDate,UpdateUser,Comments)]
trangender.parents #eventual report table

# Check that dates of birth and dates of collection make sense--------------------
# anyone enrolled before they were born
setkey(all.subj,PtCode, age, DOB, CollectDate, RecdDate)
in.utero <- all.subj[DOB>CollectDate,]
in.utero #eventual report table
impos.receive <-all.subj[RecdDate<CollectDate,]
impos.receive #eventual report table, impossible to receive a specimen before it was collected.

# Some Surgeons did not provide their style--------------
# Find missing surgeon's style
# Exclude those who came in through the support group
# Exclude those who were enrolled soon after diagnosis (meaning that at the time of enrollment they had had very few surgeries(n=0 or 1))
no.style <- affected[is.na(SurgIntervType)&InstCode!="ASRI"&SurgNumLife>2,list(CoLastName,SpecDate,PtCode,InstName)][order(CoLastName,SpecDate)]
no.style #eventual report table of affected where the surgeon did not tell us of their style, despite having more than one surgery

# missing the total number of surgeries at the time of enrollment--------------
# how is it possible
no.nr.of.surg <- affected[is.na(SurgNumLife),list(PtCode,StudyCode,CollabCode,InstCode)]
no.nr.of.surg# #eventual report table, find out why these subjets have no number of surgeries at the time of enrollment.

#some rows in the specimen table seem to be duplicates--------------

more.than.two <- all.subj.byspec[,.N,by=PtCode][N>2]
more.than.two#  eventual report table, should never have more than two specimens for a subject
# in one instance John deliberately took another specimen, told me about this in an email Fri, Dec 12, 2008 at 3:44 PM subject: Re: The curious duplication of JMS12004
setkey(all.subj.byspec,PtCode)
all.subj.byspec[more.than.two[,PtCode],list(SubjDOB, enroll.date, TissCode, hum.sub.nr, nuc.ac.nr, CoLastName, SpecDate,  CollectDate,  RecdDate,  kbddatespecim, InKbd)]

more.than.one <- all.subj.byspec[relationship!="PT",.N, by=PtCode][N>1]
more.than.one#  eventual report table, should never have more than one specimen for a subject who is a relative

# some nuc.ac.nr were entered into the database incorrectly resulting in duplicates---------------
setkey(specimens, nuc.ac.nr)
nuc.ac.nr.duplicated <- specimens[nuc.ac.nr %in% specimens[duplicated(specimens),nuc.ac.nr], list(nuc.ac.nr, SubjLastName, TissCode, RecdDate, SubjSex, InstCity)]
nuc.ac.nr.duplicated


# Key Operator miskeyed the Subject name into the HSLog----------------
# Search Specimens table to find all these potential problems
# TODO Look for PTCode where the first two digits and the three letters are differnt
# TODO look for two whole bloods or two laryngeal specimens for a single subject that would be a red flag
# TODO Look for larynx and whole blood that come on a diffent day: some indeed did but they need to be manually checked
# TODO Look for family with two mothers or two fathers
# TODO Least useful and may or may not be used would be specimens from family who are not all on the same day
# TODO Every hs# in the specimen table should be in the rpdna table unless it has not yet been processed. This would detect spurious HS# items in the specimen table that were deleted from the electronic HSlog but did not get deleted from the specimen table.
# TODO The specimens table has the nucleic acid number. Find nucleic acid numbers that have two HS#s.
# TODO The HS# should increase with CollectDate and RecdDate. Find where this is not the case since those rows in the specimen data.table need to be checked against the paper HSLog bible for possible errors






# Abberant time lapses between dates----------------------
# create a collective enrollment date which would be the minimum of SpecDate (from the rrplog) and the RecdDate and CollectDate (from specimens table). You would think SpecDate and CollectDate would always be the same but that will be a separate investigation
all.subj.byspec[,list(PtCode,TissCode,SpecDate,CollectDate,RecdDate,KbdDate,UpdateDate,difflog=CollectDate-SpecDate,shipping=RecdDate-CollectDate)] #TODO HEREIAM

# Specimens from the same subject may have different "enroll.date". Therefore for any subject the enroll.date is actually their  minimum enroll date.
specimens.on.diff.days <- all.subj.byspec[,list(engap=max(enroll.date)-min(enroll.date),min=min(enroll.date), max=max(enroll.date)),by=PtCode][engap>0]
specimens.on.diff.days#  eventual report table, must make sure that the enrollment date of the subject is the lesser of two "enrollment dates" since each specimen could have had its own enrollment date.

# There seem to be more families than affected subjects.
# Need to withdraw parents if no sibling or at least be aware
names(affec.samp)
affected.fam.names  <- affec.samp[,str_sub(PtCode,end=8L)]
anyDuplicated(affected.fam.names)
#HEREIAM
# TODO compare family names from affec.samp to family names of all.sub.samp
