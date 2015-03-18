# Trying to figure out what SNPS to use in Sequenom

library(data.table)
require(RGoogleDocs)
require(reshape2)
require(stringi)



snps.dt <- fread(input = "C:/Users/Farrel/Downloads/296autosmendelsort.csv")
snps.dt[,.N, by=SNP %like% "SNP"]
snps.dt[SNP %like% "SNP", list(SNP)][sample(.N,20)]

PLINK.snps <- data.table(read.table(file = "C:/Users/Farrel/Google Drive/RRPGenetics/296autosmendel.txt",  header = TRUE))# had to change single name from HPA#_SNP7 to HPAWTFSNP7 so that it could be read in, also had to change HPA#_SNP1 to HPAWTFSNP1, change HPA#_SNP8 to HPAWTFSNP8 .    

manifest <- fread(input = "C:/Users/Farrel/Google Drive/RRPGenetics/HumanOmni1-Quad_v1-0_H manifest.csv", sep=",",  nrows=1134514,select=c("Name", "Chr")) # nrows=1134514 is to exclude all the details about the controls that are not in the standard format. Got this from Illumina website.http://support.illumina.com/downloads/humanomni1-quadv1_product_files.html
manifest[sample(.N,60)]
sort(table(stri_extract_first(str = manifest[,Name], regex = "[^0-9]+")), decreasing=TRUE)# table of all the types of features named in the manifest.

seq.jxd  <- fread("D:/Users/Farrel/My Documents/RRPGenSucep/PedSeqJurgNov2007.txt")   #sequenome genoytping yielded by Joseph Donfack in 2007. We also merged human subject code and pedigree data into the genotyping.
seq.jxd.rs <- grep(pattern = "^rs", x = names(seq.jxd), value = TRUE) #sequenome genoytping SNP numbers yielded by Joseph Donfack in 2007. We had others which were discarded

#Get Jessica LaRusch suggestions from Google Spreadsheets-------------
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
#WARNING: this would prevent curl from detecting a 'man in the middle' attack
sheets.con = getGoogleDocsConnection(getGoogleAuth(login = getOption("GoogleDocsPassword"),service ="wise"))#the username and password is coming from the profile 
jessica=getWorksheets("2014SequenomSNPs",sheets.con)#jessica is the workbook that contains her working for determing SNPS to get for free in Sequenome
jessica.final.dt <- data.table(sheetAsMatrix(jessica$"Final group",header=TRUE, as.data.frame=TRUE, trim=TRUE, stringsAsFactors=FALSE)) #Get one sheet that is what Jessica recommended with proxies also stated
setnames(x = jessica.final.dt, make.names(names(jessica.final.dt), unique=TRUE))

tdtae.snps <- c("rs3819692", "rs12868666", "rs662041")
ever12.snps <- "rs11656744"
foxp3.snps <- c("rs3761549", "rs3761547", "rs4824747")
prior7.snps <- c(tdtae.snps, ever12.snps, foxp3.snps)

intersect(prior7.snps,jessica.final.dt[,SNP.rs.])
intersect(prior7.snps, seq.jxd.rs)
intersect(prior7.snps, manifest[,Name])
intersect(jessica.final.dt[,SNP.rs.], manifest[,Name])


# Melt all the proxies so that they can be searched for too
jessica.long <- melt.data.table(data = jessica.final.dt, id.vars = "SNP.rs.", measure.vars = 9:65, value.name = "proxy")[,variable:=NULL]
setkey(jessica.long, SNP.rs.)
names(jessica.final.dt)
length(manifest[,Name])
length(unique(jessica.long[,SNP.rs.]))

jes.prox.illu <- intersect(jessica.long[,proxy], manifest[,Name])
jes.snp.illu <- intersect(jessica.long[,SNP.rs.], manifest[,Name])
jes.final.snp <- jessica.final.dt[,SNP.rs.]


proxy.interest <- jessica.long[proxy %chin% jes.prox.illu& !SNP.rs. %chin% jes.snp.illu, list(any.one.proxy = paste0(proxy, collapse=",")), by=SNP.rs.]
proxy.interest

# Finally go for SNP.rs. for that were not on chip
jes.snp.noillu <- setdiff(jessica.final.dt$SNP.rs.,c(jes.snp.illu, proxy.interest[,SNP.rs.]))
jes.snp.noillu.butproxy <- setdiff(jessica.final.dt$SNP.rs.,c(jes.snp.illu, jes.snp.noillu))

jes.prox.noillu <- setdiff(jessica.long[,proxy], c(manifest[,Name],NA))

proxy.interest.noillu <- jessica.long[proxy %chin% jes.prox.noillu, list(any.one.proxy = paste0(proxy, collapse=",")), by=SNP.rs.]
proxy.interest.noillu

