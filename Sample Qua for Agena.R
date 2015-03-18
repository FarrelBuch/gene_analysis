# Determine a random sample of 106 to send to Sequenom for quality and quantity determination

summary(specimens)
summary(all.subj.byspec[TissCode!="LRNXP",list(nuc.ac.nr, InstCode = factor(InstCode), RecdDate, TissCode = factor(TissCode), InstName = factor(InstName), InstCountry=factor(InstCountry), spectophotom.date, dnaconcngpmicl, abs260, abs280, absratio, dilfactor, volumemicl, DNAextract.date, extractby = factor(extractby), enroll.date)])

sample.106.dna <- all.subj.byspec[TissCode!="LRNXP",list(nuc.ac.nr)][sample(.N,106)]

summary(all.subj.byspec[nuc.ac.nr %in% sample.106.dna,list(nuc.ac.nr, InstCode = factor(InstCode), RecdDate, TissCode = factor(TissCode), InstName = factor(InstName),  InstCountry=factor(InstCountry), spectophotom.date, dnaconcngpmicl, abs260, abs280, absratio, dilfactor, volumemicl, DNAextract.date, extractby = factor(extractby), enroll.date)])

save(sample.106.dna, file="nucleic acid RP qua control.RData" )

#create output for Evan to work from---------------
table.106 <- all.subj.byspec[nuc.ac.nr %in% sample.106.dna,list(nuc.ac.nr, dnaconcngpmicl, absratio, TissCode = factor(TissCode), DNAextract.date, spectophotom.date, extractby = factor(extractby))][order(nuc.ac.nr)]
write.csv(table.106, file = "sample 106 germline DNA.csv", row.names = FALSE)
                