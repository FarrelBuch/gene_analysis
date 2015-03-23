# Send all germline (bucal, mouthwash, whole blood) DNA extacts to Agena Bioscience

summary(all.subj.byspec[TissCode!="LRNXP",list(nuc.ac.nr, InstCode = factor(InstCode), RecdDate, TissCode = factor(TissCode), InstName = factor(InstName), InstCountry=factor(InstCountry), spectophotom.date, dnaconcngpmicl, abs260, abs280, absratio, dilfactor, volumemicl, DNAextract.date, extractby = factor(extractby), enroll.date)])

all.subj.byspec[TissCode!="LRNXP",list(nuc.ac.nr, TissCode = factor(TissCode), dnaconcngpmicl, abs260, abs280, absratio, dilfactor, volumemicl, extractby = factor(extractby), RecdDate, DNAextract.date, spectophotom.date  )][order(dnaconcngpmicl)][1:7]
