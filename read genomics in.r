# Read in all RRP data from wherever it is needed then there is the option to save it encrypted



# read from GENOMICS and Google---------------
# RRP log, RRP Clin Update, RRP
# PSGA, , RRP Site, HPV Type, Specimen I did not bring in Follow-up After 12
# Months since I could not find its CGS name. I probably do not need it. todo:
# still need to fetch data on google such as slected specimens for illumina
# genome wide typing

require(RGoogleDocs)
require(RODBC)
require(data.table)
require(reshape2)
require(tidyr)

# # read from GENOMICS-------------
# # on 2013-11-21 the Sybase Adaptive Server Anywhere went down
# # Untill it comes back up, I am going to comment this all out.
# # We will get the data from the saved .Rdata file
load(file="D:/Users/Farrel/My Documents/RRPGenSucep/GRRPsnapshot.RData")
#load(file="V:/GRRPsnapshot.RData")
# channellab <- odbcConnect("labdatafarrel")#connects to genomics database where most RRP data is stored, I do not know how but I got the labdata connection to work when I created it as labdatafarrel sticking all the parameters into SQL Anwywhere ("UID=fjb;PWD=***;ENG=cgs9;DBN=labdata;LINKS=TCPIP(IP=10.254.13.25)")
#
# #to find name of table in CGS database, open the table and click on "Query" and the new window will show the name as the database knows it.
#
# rrplog <- sqlQuery(channellab, "select * from CGS.RRP",stringsAsFactors=FALSE)
# clin.upd <- sqlQuery(channellab, "select * from CGS.RrpClinUpd",stringsAsFactors=FALSE)
# psga <- sqlQuery(channellab, "select * from CGS.RRP_PSGA_CFU",stringsAsFactors=FALSE)
# collab <-sqlQuery(channellab, "select * from CGS.COLLABORATOR",stringsAsFactors=FALSE)
# # hpv.type.jd <- sqlQuery(channellab, "select * from CGS.HPVTYPE") #there is no need to have this read since a copy of it is on Google Drive
# # purge corrupt rows from specimens table
# # CAR12021PT HS#4053 LRNXP received 2012-09-20 could not have been true, but instead it was HS#40353. We know that since HS# as low as 4053 took place in 1998 before RRP started. I think this occurred when it was miskeyed into HSlog but then was deleted. Deleting from HSlog apparently does not delete from Specimens table. Therefore the row with HS 4053 must be deleted
# # RYS16057PT has the exact same problem as CAR12021PT. There are two laryngeal specimens. One is with HS 40338 and the other with HS 4038, HS 4038 was before we started collecting, therefore we must delete the HS 4038 row from our specimens data.table
# # we have two whole blood specimens for LJS09014PT. In the specimens table we have both received by us on 2013-04-16. One is HS# 41270 and the other is HS#41244. His laryngeal papilloma specimen was HS# 41271. HS 41270 is the correct one, 41244 is a completely different sample from the METRC study. Therefore the row with HS 41244 must be deleted. The study code was changed before 2013-10-28 17:48
# specimens <- sqlQuery(channellab,"SELECT * from CGS.Specimen WHERE StudyCode= 'RP'AND (HS#<>4053 AND HS#<>4038)", believeNRows = FALSE,stringsAsFactors=FALSE,na.strings="") #the HS# part of the statement is to avoid reading in corrupt rows
# rrp.related <- sqlQuery(channellab,"SELECT * from CGS.Specimen WHERE StudyCode= 'RPREL'", believeNRows = FALSE,stringsAsFactors=FALSE)

# read from Google--------------------
# Monday, 1 Jun 2015 21:40 RGoogleDocs not working, will transition to googlesheets


# options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
# #WARNING: this would prevent curl from detecting a 'man in the middle' attack
# sheets.con = getGoogleDocsConnection(getGoogleAuth(login = getOption("GoogleDocsPassword"),service ="wise"))#the username and password is coming from the profile
# hpv.google=getWorksheets("hpv type",sheets.con)#hpv.google is the HPV typing data that was stored in Google
# hpv.early <- sheetAsMatrix(hpv.google$"Sheet 1",header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE) #Get one sheet that is what was typed in early years by Joseph Donfack using allele specific PCR and restriction fragrment length polymorphism
# hpv.provtous <- sheetAsMatrix(hpv.google$"provided to us",header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE) #Get one sheet that is what other offices or labs provided to us for HPV typing
#
# rpdna.google=getWorksheets("rpdna",sheets.con)#rpdna is the extracted DNA details and location from each is stored in Google
# rpdna<-sheetAsMatrix(rpdna.google$dnastock,header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE)
# mergeddata3=getWorksheets("mergeddata3",sheets.con)#spreadsheet supplied to Jurg Ott to run PLINK on the 296 samples that first underwent genotyping
# PLINK.pedigree<-sheetAsMatrix(mergeddata3$"Main Sheet",header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE)
# rpinfinwrk.google <- getWorksheets("rpinfinwrk",sheets.con)#rpinfinwrk is the workflow that we used to do the infinium assay in preparation for the HumanOmni1-Quad DNA analysis
# rpinfinwrk<-sheetAsMatrix(rpinfinwrk.google$"Main Sheet",header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE)
#
# specs.post.crash.google <- getWorksheets("specimens post database crash",sheets.con) # "specimens post database crash" are the specimens received after GENOMICS went down and now stored in Google
# specs.post.crash <- as.data.table(sheetAsMatrix(specs.post.crash.google$Sheet1, header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE))
# specs.post.crash[,`:=`(nuc.ac.nr = as.integer(nuc.ac.nr),
#                        CollectDate = as.Date(CollectDate,origin="1899-12-30"),
#                        RecdDate = as.Date(RecdDate,origin="1899-12-30"))]
# class(specs.post.crash$CollectDate)
# setkey(specs.post.crash, SubjLastName)
#
# # Read Satellite Hospitals
# satellites.google=getWorksheets("Satellite Hospitals",sheets.con) #satellites.google is a listing of each satellite to link PI to institution and date of IRB approval. Spreadsheet held in  Google
# satellites <- as.data.table(sheetAsMatrix(satellites.google$Hospitals,header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE)) #satellites is a data.table of each satellite to link PI to institution and date of IRB approval. Spreadsheet held in  Google
# satellites[,IRBApproval := as.IDate(IRBApproval,origin="1899-12-30")]
# setnames(satellites, old="Hospital Name", new="hospital")
# approved.sat  <- satellites[!is.na(IRBApproval),list(PI, PI.first,  InstCode, hospital, IRBApproval, City, Country, zip, FWA, ArmstongPap)] #approved.sat is Institutional Review Board (IRB) approved satellites
# rm(satellites.google)
#
# # Read Collaborators curated from Farrel
# coinvest.google=getWorksheets("coinvestigator contact details",sheets.con) #coinvest.google is a listing of each prinicipal investigator no matter how much enrollment. It was generated from Google Contacts
# coinvest <- as.data.table(sheetAsMatrix(coinvest.google$coinvest,header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE)) #coinvest is a data.table of each pi to link PI to emial and held in  Google

# googlesheets way #HEREIAM annotate this all from above
library(googlesheets)
hpv.type <- gs_title("hpv type")
hpv.early <- data.table(gs_read(ss = hpv.type , ws = "Sheet 1"))
hpv.provtous <- data.table(gs_read(ss = hpv.type , ws = "provided to us"))
rpdna <- data.table(gs_read(ss = gs_title("rpdna") , ws = "dnastock"))
PLINK.pedigree <- data.table(gs_read(ss = gs_title("mergeddata3") , ws = "Main Sheet"))
rpinfinwrk <- data.table(gs_read(ss = gs_title("rpinfinwrk") , ws = "Main Sheet"))
specs.post.crash <- data.table(gs_read(ss = gs_title("specimens post database crash") , ws = "Sheet1"))
specs.post.crash[,`:=`(nuc.ac.nr = as.integer(nuc.ac.nr),
											 CollectDate = as.Date(CollectDate,origin="1899-12-30"),
											 RecdDate = as.Date(RecdDate,origin="1899-12-30"))]
setkey(specs.post.crash, SubjLastName)
satellites <- data.table(gs_read(ss = gs_title("Satellite Hospitals") , ws = "Hospitals"))
satellites[,IRBApproval := as.IDate(IRBApproval,origin="1899-12-30")]
setnames(satellites, old="Hospital.Name", new="hospital")
approved.sat  <- satellites[!is.na(IRBApproval),list(PI, PI.first,  InstCode, hospital, IRBApproval, City, Country, zip, FWA, ArmstongPap)] #approved.sat is Institutional Review Board (IRB) approved satellites
coinvest <-  data.table(gs_read(ss = gs_title("coinvestigator contact details") , ws = "coinvest"))
setnames(coinvest,make.names(names(coinvest)))
# get rid of all columns besides name and email
coinvest[,c(5:27,36:88) := NULL]
coinvest.long  <- reshape(data=coinvest, direction="long" , varying = list(c( "E.mail.1...Type", "E.mail.2...Type", "E.mail.3...Type", "E.mail.4...Type"), c( "E.mail.1...Value", "E.mail.2...Value", "E.mail.3...Value", "E.mail.4...Value")))
setkey(coinvest.long,Name)
coinvest.long[,c("time", "id") := NULL]
setnames(coinvest.long, old=c("E.mail.1...Type", "E.mail.1...Value"), new=c("email.type","email"))
coinvest.long  <-  coinvest.long[!((is.na(email)|grepl("old",x=email.type, ignore.case=TRUE ))),]
coinvest.wide.mail <- coinvest.long[,list(emails=paste(email, collapse = ",")), by=list(PI=Family.Name,PI.first=Given.Name)]
coinvest.long  <- coinvest.long[,list(PI=Family.Name, PI.first=Given.Name,email )]

# get ready to attach emails to approved.sat
setkey(approved.sat, PI, PI.first )
setkey(coinvest.wide.mail,  PI, PI.first)
setkey(coinvest.long,  PI, PI.first)




# create key for consent codes
# the consent codes were  entered into the rrplog for the patients
# would need to go to paper charts to find consent codes for relatives
consent <- data.table(ConsentCode=0:3,consent.descr=c("The blood and the wart specimen may only be used for this research project",
                                                      "The blood and the wart specimen may be used for any unspecified research project, including genetic reseach. There is no need to contact me",
                                                      "The blood and the wart specimen may be used for any unspecified research project, but not for genetic research. There is no need to contact me",
                                                      "Please contact me to obtain consent for the blood and the wart specimen to be used in a specific future project"),key="ConsentCode")

# read from RData on my hard drive----------------
load(file="C:/Users/Farrel/Google Drive/RRPGenetics/calls by ABI cluster.Rdata")

# cleaning up electronic data entry typographical issues-----------
# I detected that there were a few RP nuc.ac.nr that were duplicates. Evan reconciled against paper bible on Monday, 23 Mar 2015.
# https://mail.google.com/mail/ca/u/0/#inbox/14c0050a5c948790
specimens$`NA#`[specimens$SubjLastName=="AJD18001PT"&specimens$TissCode=="LRNXP"] <- 1612
specimens$`NA#`[specimens$`NA#`==1611&specimens$SubjLastName=="FGD01004PT"&specimens$TissCode=="LRNXP"] <- 1613
specimens <- specimens[-which(specimens$`NA#`==1751)[2],] # this deals with the second instance of RP1751 which is a duplicate except for the date and time of entry into the database.
specimens$`NA#`[specimens$`NA#`==1882&specimens$SubjLastName=="AMW26009PT"&specimens$TissCode=="LRNXP"] <- 1885
specimens$`NA#`[specimens$`NA#`==1952&specimens$SubjLastName=="VDU29043PT"&specimens$TissCode=="LRNXP"] <- 1954
specimens$`NA#`[specimens$`NA#`==1612&specimens$SubjLastName=="FGD01004PT"&specimens$TissCode=="WB"] <- 1614



# some final issues and notification of being done-------------


all.data.frames <- ls()[sapply(mget(ls(), .GlobalEnv),is.data.frame)]#this is how we get to manipulate all the data.frames

#save(list=all.data.frames,file="V:/GRRPsnapshot.RData")
#save.image("D:/Users/Farrel/My Documents/RRPGenSucep/GRRPsnapshot.RData")

shell.exec("http://upload.wikimedia.org/wikipedia/commons/e/e1/Activated_fire_alarm_%28sound%29.ogg")#play a file when done
