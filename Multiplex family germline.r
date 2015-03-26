# Multiplex samples


multiplexes <- c("CMM12003PT","CMM12003FA","CMM12003MO","RTC13002PT","RTC13002MO","FJB01123PT","FJB01123FA","FJB01120PT","FJB01123MO","FJB01120MO")# these are the people who come from a multiplex family
all.subj.byspec[TissCode!="LRNXP"& PtCode %chin% multiplexes,list(PtCode, nuc.ac.nr, TissCode = factor(TissCode), dnaconcngpmicl, abs260, abs280, absratio, DNAextract.date)][order(-TissCode)]

