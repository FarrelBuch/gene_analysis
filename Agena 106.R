# Analyze Agena Quality and Quantity Control assay on 106 specimens.

library(data.table)
library(readxl)
library(zoo)


# Read their report
agena.106 <- as.data.table(read_excel("Allegheny_106-samples.xlsx", sheet = "Allegheny_106-samples"))
setnames(agena.106, make.names(names(agena.106), unique=TRUE))
names(agena.106)
agena.106[, ':='(position = paste0(Plate,Well, sep=""),
								 Plate = NULL,
								 Well = NULL)] #position is concatenation of plate and well so that it can be merged between two tables
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
setnames(agena.evan, old=c("Amplifiable.Copies..Avg.", "Esimated.DNA.amount..ng.", "Amplifiable.Copies..SD.", "concentration..ng.ul.", "total.ng.DNA", "total.volume..ul.", "sample.name"), new=c("qc.amplif.copies", "qc.est.dna.ng", "qc.amplif.copies.sd", "ahn.stock.conc.ngperul", "dna.sent.ng", "vol.sent.ul", "nuc.ac.nr"))


# compare AHN belief and Agena quality control
agena.evan[,":="(rankAHN = rank(ahn.stock.conc.ngperul),
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

agena.ahn.106$ahn.stock.conc.ngperul

ggplot(data = agena.ahn.106,  aes(y =  log(qc.amplif.copies), x=log(ahn.stock.conc.ngperul) )) + geom_point(aes(colour=TissCode), size=3, alpha=0.6)
class(agena.ahn.106$ahn.stock.conc.ngperul)
class(agena.ahn.106$qc.amplif.copies)


#describe data set wrt pass or fail qc

feats <- quote(c(N=.N,
								 as.list(range(qc.amplif.copies))))
agena.ahn.106[,eval(feats), by=TissCode]
agena.ahn.106[,list(N=.N,
										min.stock.ngperul = min(ahn.stock.conc.ngperul),
										min.avg.copy = min(qc.amplif.copies),
										perc10.stock.ngperul = quantile(ahn.stock.conc.ngperul, probs=0.1),
										perc10.avg.copy = quantile(qc.amplif.copies, probs=0.1),
										perc50.avg.copy = median(qc.amplif.copies),
										max.avg.copy = max(qc.amplif.copies)), keyby=list(Passed.QC,TissCode)]
