# Analyze Agena Quality and Quantity Control assay on 106 specimens.

library(data.table)
library(readxl)
library(zoo)


# Read their report
agena.106 <- as.data.table(read_excel("Allegheny_106-samples.xlsx", sheet = "Allegheny_106-samples", na="NA"))
setnames(agena.106, make.names(names(agena.106), unique=TRUE))
setnames(agena.106, old="Gender", new="qc.gender")
agena.106[, ':='(position = paste0(Plate,Well, sep=""),
								 Plate = NULL,
								 Well = NULL)] #position is concatenation of plate and well so that it can be merged between two tables
agena.106[, qc.est.conc := Amplifiable.Copies..Avg. * (1/320) / 2] #concenteration in nanograms per microliter
# Robin Everts said "We always take about 2 uL to start PCR with so the data in column C is based on that, Assuming 10 ng is 3200 copies"

setkey(agena.106, position)


# Read what what Evan Sent them
evan.sent <- as.data.table(read_excel("C:/Users/Farrel/Google Drive/RRPGenetics/gene_analysis/Agena_106samples Evan Sent.xlsx"))
setnames(evan.sent, make.names(names(evan.sent), unique=TRUE))
setnames(evan.sent, old=1, new="plate")
evan.sent$plate <- na.locf(evan.sent$plate)
evan.sent$plate <- substr(evan.sent$plate, start=7, stop=7)
evan.sent[, ':='(position = paste0(plate,well, sep=""),
								 plate = NULL,
								 well = NULL)] #position is concatenation of plate and well so that it can be merged between two tables
evan.sent[, concentration..ng.ul. := as.numeric(concentration..ng.ul.)]
setkey(evan.sent, position)


agena.evan <- merge(agena.106, evan.sent)
setnames(agena.evan, old=c("Amplifiable.Copies..Avg.", "Esimated.DNA.amount..ng.", "Amplifiable.Copies..SD.", "concentration..ng.ul.", "total.ng.DNA", "total.volume..ul.", "sample.name"), new=c("qc.amplif.copies", "qc.est.dna.ng", "qc.amplif.copies.sd", "ahn.sent.conc.ngperul", "dna.sent.ng", "vol.sent.ul", "nuc.ac.nr"))


# compare AHN belief and Agena quality control
agena.evan[,":="(rankAHN = rank(ahn.sent.conc.ngperul),
								 rankqc = rank(qc.amplif.copies))]
agena.evan[, rank.diff := rankAHN-rankqc]
summary(agena.evan$rank.diff)
agena.evan[, nuc.ac.nr := as.integer(nuc.ac.nr)]
setkey(agena.evan, nuc.ac.nr)

all.subj.byspec[, nuc.ac.nr := as.integer(nuc.ac.nr)]
setkey(all.subj.byspec, nuc.ac.nr)

# merge AHN belief and Agena quality control with our specimens
agena.ahn.106 <- merge(x = agena.evan, y = all.subj.byspec, all.x = TRUE, all.y = FALSE)



# plot the AHN vs Agena
library(ggplot2)

# lets looks at some data
agena.ahn.106[ ,list(position, TissCode, ahn.sent.conc.ngperul, absratio, qc.est.conc, qc.amplif.copies)][order(ahn.sent.conc.ngperul)][1:25] #sorted by spectrophotometry


ggplot(data = agena.ahn.106,  aes(y =  log10(qc.est.conc), x=log10(ahn.sent.conc.ngperul) )) +
	geom_point(aes(colour=TissCode), size=3, alpha=0.6) +
	ggtitle("Comparing concentraton of DNA sent out based on spectrophotometry\n vs estimated concentration based on Sequenom")



#describe data set wrt pass or fail qc

feats <- quote(c(N=.N,
								 as.list(range(qc.amplif.copies))))
agena.ahn.106[,eval(feats), by=TissCode]
agena.ahn.106[,list(N=.N,
										min.sent.ngperul = min(ahn.sent.conc.ngperul),
										min.avg.qc.conc = min(qc.est.conc),
										perc10.sent.ngperul = quantile(ahn.sent.conc.ngperul, probs=0.1),
										perc10.avg.qc.conc = quantile(qc.est.conc, probs=0.1),
										perc50.avg.qc.conc = median(qc.est.conc),
										max.avg.qc.conc = max(qc.est.conc)), keyby=list(Passed.QC,TissCode)]

# describe data wrt snp calls
# NA means no result
agena.ahn.106[, missed.snps:=rowSums(is.na(.SD)), .SDcols=SIDv1_SNP01:SIDv1_SNP44] #tip found a way to add a column that counts the number of NA's across a row.

ggplot(data=agena.ahn.106, aes(x=missed.snps)) +
	geom_histogram(aes(y=..density..)) +
	facet_wrap(~TissCode) +
	ggtitle("Density of missing SNP calls out of 44\n1 buccal, 22 mouthwash, 83 whole blood")

# how many missed calls
print(agena.ahn.106[,list(.N), keyby=missed.snps][,cumulative:=round(cumsum(N)/agena.ahn.106[,.N],3)])
# how many missed calls amongst those with low stock concentrations
print(agena.ahn.106[dnaconcngpmicl<36,list(.N), keyby=missed.snps][,cumulative:=round(cumsum(N)/agena.ahn.106[,.N],3)])# we discovered in this sample of 106 that 75% of stock solutions with concentations of <36 ng/µl have 0, 1 or 2 missed SNPs out of 44.
# todo figure out what proportion of all stocks are <36 ng/µl

12/106


# lets show the characteristics of the samples with missing SNPS in descending order
agena.ahn.106[,list(position, TissCode, ahn.sent.conc.ngperul, absratio, qc.est.conc, qc.amplif.copies, missed.snps)][order(-missed.snps)][1:20]

# does number of missed calls related to concentration in our stock solution or perhaps our absorbanc ratio

ggplot(data = agena.ahn.106,  aes(y =  missed.snps, x=log10(dnaconcngpmicl) )) +
	geom_point(aes(colour=TissCode), size=3, alpha=0.6) +
	ggtitle("Comparing concentraton of DNA of our stock based on sectrophotometry\n vs number of missed SNP calls out of 44 tested")


ggplot(data = agena.ahn.106,  aes(y =  missed.snps, x=absratio )) +
	geom_point(aes(colour=TissCode), size=3, alpha=0.5) +
	geom_vline(xintercept = 1.8)
	ggtitle("Comparing absorbance ratio of our stocks \n vs number of missed SNP calls, line at 1.8")


summary(lm(missed.snps~absratiolow, data = agena.ahn.106))
agena.ahn.106[,list(missed.snps, absratio)][order(-missed.snps)]
agena.ahn.106[absratio>2.5, list(qc.amplif.copies, agena.est.conc=round(qc.est.conc,2), missed.snps, dnaconcngpmicl, abs260,  abs280,  absratio,  DNAextract.date,  nuc.ac.nr)]

summary(lm(missed.snps~dnaconcngpmicl, data = agena.ahn.106))
summary(lm(missed.snps~qc.est.conc, data = agena.ahn.106))
summary(lm(missed.snps~dnaconcngpmicl+TissCode, data = agena.ahn.106[TissCode!="BUC"]))
summary(lm(missed.snps~dnaconcngpmicl*TissCode, data = agena.ahn.106[TissCode!="BUC"]))
summary(lm(missed.snps~dnaconcngpmicl, data = agena.ahn.106[TissCode!="BUC"]))

summary(aov(formula = qc.est.conc ~ cut(missed.snps, breaks = c(-1, 0.5, 2.5, 44)), data = agena.ahn.106))
ahn.conc.miss <- aov(formula = ahn.sent.conc.ngperul ~ cut(missed.snps, breaks = c(-1, 0.5, 2.5, 44)), data = agena.ahn.106)
summary(ahn.conc.miss)
TukeyHSD(ahn.conc.miss)
library(gplots)
plotmeans(agena.ahn.106$ahn.sent.conc.ngperul ~ cut(agena.ahn.106$missed.snps, breaks = c(-1, 0.5, 2.5, 44), labels = c("none", "one or two", "several or lots")), xlab="number of missing SNPS",
					ylab="concentration of DNA sent out by AHN (ng/µl)", main="The more SNPS that were missed the lower the average concentration")

summary(aov(formula = qc.est.conc ~ cut(missed.snps, breaks = c(-1, 0.5, 2.5, 44)), data = agena.ahn.106))
ahn.conc.miss <- aov(formula = dnaconcngpmicl ~ cut(missed.snps, breaks = c(-1, 0.5, 2.5, 44)), data = agena.ahn.106)
summary(ahn.conc.miss)
TukeyHSD(ahn.conc.miss)

library(gplots)
plotmeans(agena.ahn.106$dnaconcngpmicl ~ cut(agena.ahn.106$missed.snps, breaks = c(-1, 0.5, 2.5, 44), labels = c("none", "one or two", "several or lots")), xlab="number of missing SNPS",
					ylab="concentration of DNA stock held by AHN (ng/µl)", main="The more SNPS that were missed the lower the average concentration")

# check gender
agena.ahn.106[,list(nuc.ac.nr, qc.gender, SubjSex, gender.same=qc.gender==SubjSex)][,.N,by=gender.same]
