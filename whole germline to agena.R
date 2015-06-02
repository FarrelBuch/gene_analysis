# Send all germline (bucal, mouthwash, whole blood) DNA extacts to Agena Bioscience

library(data.table)

summary(all.subj.byspec[TissCode!="LRNXP",list(nuc.ac.nr, InstCode = factor(InstCode), RecdDate, TissCode = factor(TissCode), InstName = factor(InstName), InstCountry=factor(InstCountry), spectophotom.date, dnaconcngpmicl, abs260, abs280, absratio, dilfactor, volumemicl, DNAextract.date, extractby = factor(extractby), enroll.date)])

all.subj.byspec[TissCode!="LRNXP",list(nuc.ac.nr, TissCode = factor(TissCode), dnaconcngpmicl, abs260, abs280, absratio, dilfactor, volumemicl, extractby = factor(extractby), RecdDate, DNAextract.date, spectophotom.date  )][order(dnaconcngpmicl)][1:7]

plates <- ceiling(nrow(all.subj.byspec[TissCode!="LRNXP"])/95)

plating <- all.subj.byspec[TissCode!="LRNXP",list(nuc.ac.nr,
																						 dnaconcngpmicl,
																						 absratio,
																						 TissCode = factor(TissCode),
																						 extractby = factor(extractby))][order(dnaconcngpmicl)][, ':=' (plate=rep(1:plates, each=95),
																						 																															 position=rep(1:95, times=plates))]


# Generate random spots for water------------
water.holes <- sample(1:96, size = plates, replace = TRUE)

plating[, water.hole:=water.holes[plate], by=plate]
plating[, real.position:=position]#set it up
plating[position>=water.hole, real.position := position+1]

# now add water rows
plating <- rbind(plating,data.table(nuc.ac.nr = rep(0, length.out=plates), plate=seq(from=1,to=plates), water.hole=water.holes, real.position=water.holes), use.names=TRUE, fill = TRUE)
setkey(plating,plate,real.position)

# well nomenclature
plating[,well:=paste0(rep(LETTERS[1:8],times=12),rep(1:12, each=8))]

plating[,':='(position=NULL, water.hole=NULL)]
setnames(plating, old="real.position", new="position")
setcolorder(plating, c("plate",  "well",  "position",  "nuc.ac.nr",  "dnaconcngpmicl",  "absratio",  "TissCode",  "extractby"))
plating[, dnaconcngpmicl := round(dnaconcngpmicl)]

write.csv(plating, file = "ginormous germline to Agena.csv", row.names = FALSE)


# plate concentration characteristics
conc.by.plate <- plating[,list(min.conc=min(dnaconcngpmicl, na.rm = TRUE), median=median(dnaconcngpmicl, na.rm = TRUE), max.conc=max(dnaconcngpmicl, na.rm = TRUE)), by=plate]
write.csv(conc.by.plate, "concentration by plate.csv", row.names = FALSE)



#save(plates, plating, water.holes, file = "plating data.RData")



