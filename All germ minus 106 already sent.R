# Send all germline (bucal, mouthwash, whole blood) DNA extacts but exclude the 106 that were already sent to Agena Bioscience

library(data.table)

load(file="nucleic acid RP qua control.RData")

all.subj.byspec[TissCode!="LRNXP"&!nuc.ac.nr %in% sample.106.dna$nuc.ac.nr ,list(nuc.ac.nr, TissCode = factor(TissCode), dnaconcngpmicl, abs260, abs280, absratio, dilfactor, volumemicl, extractby = factor(extractby), RecdDate, DNAextract.date, spectophotom.date  )][order(dnaconcngpmicl)]

plates <- ceiling(nrow(all.subj.byspec[TissCode!="LRNXP"&!nuc.ac.nr %in% sample.106.dna$nuc.ac.nr])/95)

plating <- all.subj.byspec[TissCode!="LRNXP"&!nuc.ac.nr %in% sample.106.dna$nuc.ac.nr,list(nuc.ac.nr,
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




# plate concentration characteristics
conc.by.plate <- plating[,list(min.conc=min(dnaconcngpmicl, na.rm = TRUE), median=median(dnaconcngpmicl, na.rm = TRUE), max.conc=max(dnaconcngpmicl, na.rm = TRUE)), by=plate]
write.csv(conc.by.plate, "concentration by plate minus 106.csv", row.names = FALSE)





# DNA to add--------------
conc.by.plate[median<=20, dnatoadd:=6]#volume of stock to add
conc.by.plate[median>20, dnatoadd:=5]
conc.by.plate[,watertoadd:=round(((median*dnatoadd)/20)-dnatoadd)]
conc.by.plate[watertoadd<0,watertoadd:=0]
# merge median plate details into plating
setkey(plating, plate)
setkey(conc.by.plate, plate)
conc.by.plate[plating]
plating <- plating[conc.by.plate[,list(plate,dnatoadd, watertoadd)]]
plating[,totalngdna := dnaconcngpmicl*dnatoadd]
plating[,totalvol := dnatoadd+watertoadd]
plating[,dilconcngpµl := totalngdna/totalvol]
setcolorder(plating, c("plate",  "well",  "position",  "nuc.ac.nr",  "dnaconcngpmicl", "dnatoadd", "watertoadd", "dilconcngpµl", "totalngdna", "totalvol", "absratio",  "TissCode",  "extractby"))
write.csv(plating, file = "ginormous germline minus 106 to Agena.csv", row.names = FALSE)
#save(plates, plating, water.holes, file = "plating data minus 106.RData")

