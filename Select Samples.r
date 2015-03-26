##THIS WILL SAMPLE THE SPECIMEN COLLECTION FOR THE FIRST 100 probands and their relatives
##prefer whole blood over mouth wash
#
#Prerequisites
#Run PedigreeKD.r to get the data frame ?PedigreeTable?
#
#Goals
#Family name by trio status
#Family name by whole blood samples, laryngeal bx
#Randomly order
# 
#Output rp nucleic acid numbers so that they can be picked
#must have whole family together.

#Family name by holding score
require(RODBC)

channellab <- odbcConnect("Labdata")# connect to Genomics now "labdata" via ODCB 
holdings <-subset(sqlQuery(channellab, "SELECT HS#, SubjLastName, TissCode, StudyCode from CGS.Specimen WHERE STUDYCODE='RP'",believeNRows = FALSE),select=-StudyCode)
require(gdata)
holdings <-drop.levels(holdings)
head(holdings)
require(car)

holdings <-transform(holdings, FA=grepl("FA", substr(holdings$SubjLastName,9,10)),MO=grepl("MO", substr(holdings$SubjLastName,9,10)), tiss.score=recode(TissCode,"'WB'=10;'MW'=4;'BUC'=3;'LRNXP'=5;else=0", as.factor.result=FALSE))

require(plyr)
familyname <-function(x) substr(x,1,8)#function to get family name out of subject number
sum(holdings$MO,holdings$FA)

fam.holdings <-ddply(holdings,.(substr(holdings$SubjLastName,1,8)), function(df) c(tiss.score=sum(df$tiss.score), trio=1+sum(df$FA,df$MO), indivs=length(unique(df$SubjLastName))))#adds all the tissue scores to get a sum for the whole family; add another function to find out how many unique people, if father, if mother, trio, duo, solo, how many non mothers and fathers)
names(fam.holdings)[1] <-"familyname"
fam.holdings$trio <-factor(fam.holdings$trio, labels=c("solo","duo","trio"))

addmargins(table(fam.holdings$trio))
prop.table(table(fam.holdings$trio))
fam.holdings <-transform(fam.holdings, value=((as.numeric(trio)-1)^3*tiss.score*tiss.score/indivs))
sample <-sample(1:nrow(fam.holdings),120*46/56,prob=fam.holdings$value)
table(fam.holdings[sample,"trio"])

#Now that we have the families we need the individual specimens.
#It must be everyone of a family but not only the whole blood and the mouthwash
illumina.specimens <-subset(holdings,familyname(SubjLastName) %in% fam.holdings$familyname[sample]& (TissCode=="MW" | TissCode=="WB"),select=c(HS.,SubjLastName)) 
ord <-order(illumina.specimens$SubjLastName)
illumina.specimens <-illumina.specimens[ord,]


#Hold on since we have already done specimens, 
#first.run <-c(1078, 1081, 1080, 543, 545, 672, 546, 673, 674, 934, 933, 935, 68, 712, 312, 836)
#We now need to exclude everything that was already done




#Need to get nucleic acid numbers for all members of selected families
#import used to be from Genomics but since nanodrop to genomics macro broke will rather take it from google docs

library(RGoogleDocs)

options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
sheets.con <- getGoogleDocsConnection(getGoogleAuth(login = getOption("GoogleDocsPassword"),service ="wise")) 
ts2=getWorksheets("rpdna",sheets.con)
dna<-sheetAsMatrix(ts2$dnastock,header=TRUE, as.data.frame=TRUE, trim=TRUE)
#now merge dna with illumina.speicmens
rp.infin.wrk <-merge(illumina.specimens, dna, by.x="HS.", by.y="hsnr", all.x=TRUE)




#Pulling out specific families for genotyping
multiplexes <- c("CMM12003PT","CMM12003FA","CMM12003MO","RTC13002PT","RTC13002MO","FJB01123PT","FJB01123FA","FJB01120PT","FJB01123MO","FJB01120MO")#individuals from a multiplex family
# specimens (and holdings) table has subjectlastname and HS# and the specimen type (blood vs mouthwash vs larynx), nanodrop rpdna (here already imported as dna) has the link of hs# to dna#
#strategy
#get specimen table (already have in holdings)
#get rpdna    (already have in dna)
#merge them on HS number and keep HS, nucleic acid number, specimen tissue
hsdnanrs <-merge(holdings,dna,by.x="HS.",by.y="hsnr")
hsdnanrs.notlr <-subset(hsdnanrs,TissCode!="LRNXP",select=c(HS.:rpdnanr))#merged data but without laryngeal specimens
multiplexes.dnanr <-subset(hsdnanrs,TissCode!="LRNXP"&SubjLastName %in% multiplexes,select=(names(rp.infin.wrk)))
#make sure that multiplexes are to be genotyped even if they have not been selected by random sampling
rp.infin.wrk <-rbind(multiplexes.dnanr,rp.infin.wrk)



#Hold on since we have already done specimens, 
#already done
ts1=getWorksheets("rpinfinwrk",sheets.con)[[1]]# the [[1]] simply gets the first worksheet
#getExtent(ts1)[2,1] #getExtent enables us to determine the number of cells with data in it
already.run <-as.vector(as.matrix(ts1[2:getExtent(ts1)[2,1],1])[,1])#get the vector of rpdna nucleic acid number that have been used already

rp.infin.wrk <-subset(rp.infin.wrk, !is.element(rpdnanr,already.run),select=rpdnanr)
write.csv(rp.infin.wrk,"rpdnanr.csv")