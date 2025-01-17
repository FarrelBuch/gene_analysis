#Sample for HLA typing 200 specimens

#Criteria

## Must not have been done already
## Must have good DNA (exclude stock specimens < 30 ng/µl)
## Must be part of trio
## enrich for clear HPV type
## enrich for illumina
## indiferent if whole exome sequencing
## must not be in mendel.error.fam

source(file = "exome sequence sample selection.r")

# read in HLA typing already done----------
library(readxl)
first.hla <- as.data.table(read_excel("20131023TCA results checked.xls"))
setkey(first.hla, PID)
first.hla <- first.hla[!is.na(PID)]
first.hla[,nuc.ac.nr:=as.numeric(PID)]
setkey(first.hla, nuc.ac.nr)
setkey(all.subj.byspec, nuc.ac.nr)
first.hla[all.subj.byspec, PtCode:= i.PtCode]# read about this at data.table join then add columns to existing data.frame without re-copy at http://stackoverflow.com/questions/19553005/data-table-join-then-add-columns-to-existing-data-frame-without-re-copy
first.hla[all.subj.byspec, family:= i.family]


# use exome squence selection but alter formula for HLA criteria
rank.hla.affected <- ordered.exome[,rank.hla:=rank.exome/rank.express]
rank.hla.affected[order(-rank.hla)]

# bind rows of affected and unaffected
rank.hla.fams <- rbindlist(list(rank.hla.affected, DNA.rank[relationship!="PT"]), use.names = TRUE, fill=TRUE)
rank.hla.fams[,family:=str_sub(PtCode,end=8L)]




# Exclude some families--------------------

# ensure no mendelian error families
rank.hla.fams <- rank.hla.fams[!family %chin% mendel.error.fam]

# remove HLA already done
rank.hla.fams <- rank.hla.fams[!nuc.ac.nr %in% first.hla$nuc.ac.nr]

# remove  "SI"  "MGM" "MGF" for now
rank.hla.fams <- rank.hla.fams[!relationship %chin% c("SI" , "MGM", "MGF")]

# remove duo and solo
trio.families <- rank.hla.fams[trio.duo.sol=="trio", list(family)]
rank.hla.fams <- rank.hla.fams[family %chin% trio.families$family]

# remove low concentration stock specimens
rank.hla.fams <- rank.hla.fams[dnaconcngpmicl>=30]


#rank order desirability
#rank.hla.fam multiplys rank.hla of affected by the dna.rank of each family member.
rank.hla.fams[,rank.hla.fam:=(prod(dna.rank)*max(rank.hla, na.rm=TRUE)), by=family]#rank.hla.fam multiplys rank.hla of affected by the dna.rank of each family member.
rank.hla.fams[order(-rank.hla.fam,family)]

# Generate list of suitable samples-------------------------
hla.200.fams <- unique(rank.hla.fams[rank.hla.fam>0,list(family, rank.hla.fam)])[sample(.N,67, replace=FALSE, prob=rank.hla.fam)]
nuc.ac.hla.200.list <- rank.hla.fams[family %chin% hla.200.fams$family, nuc.ac.nr]

# find the samples for evan----------------
load(file="plating data minus 106.RData")
setkey(plating,nuc.ac.nr)
gofind.hla <- plating[nuc.ac.nr %in% nuc.ac.hla.200.list,list(nuc.ac.nr,plate,well)][order(plate,well)]#samples in the big agena send out

gofind.hla.106 <- setdiff(nuc.ac.hla.200.list, plating$nuc.ac.nr)#samples that were sent out in the initial 106
go.find.all.hla <- rbindlist(list(gofind.hla, data.table(nuc.ac.nr=gofind.hla.106)), fill=TRUE, use.names=TRUE)
setkey(go.find.all.hla,nuc.ac.nr)

# merge with data requested by Carrington
master.hla.200 <- all.subj.byspec[,.SD, .SDcols=c(nuc.ac.nr,dnaconcngpmicl:volumemicl, dnayieldmicrog, family, relationship, PtCode)][go.find.all.hla]
setkey(master.hla.200,PtCode)

master.hla.200 <- all.subj[,list(PtCode, Sex, Ethnicity, Race)][master.hla.200]

# pipette workflow------------------
# never pipette <5 µl of stock
# never send < 20 µl of final volume

# for samples with stock concentration of 30 ng/µl to 50 ng/µl we will provide not just 500 ng but 1000 ng
# for the rest of the samples we will give 750 ng at concentration of 50 - 100 ng/µl
master.hla.200[dnaconcngpmicl<50, pip.stock:=signif(1000/dnaconcngpmicl,2)]
master.hla.200[dnaconcngpmicl>=50, pip.stock:=signif(750/dnaconcngpmicl,2)]
master.hla.200[pip.stock<5, pip.stock:=5]

master.hla.200[dnaconcngpmicl<50, pip.diluent:=0]
master.hla.200[dnaconcngpmicl>=50, pip.diluent:=signif(((pip.stock*dnaconcngpmicl)/75)-pip.stock, digits=2)]
master.hla.200[(pip.diluent+pip.stock)<20, pip.diluent:=signif(((pip.stock*dnaconcngpmicl)/50)-pip.stock, digits=2)]
master.hla.200[pip.stock+pip.diluent<20,list(max(dnaconcngpmicl))]

write.csv(master.hla.200, file = "master hla 200.csv", row.names = FALSE)
#save(master.hla.200,nuc.ac.hla.200.list, file = "hla 200.RData")#contains two objects: a vector of nuc.ac.nr of the random sample enriched for most suitable to undergo HLA typing, the entire master.hla.200 that includes PtCode, demographics etc

master.hla.200[,.N,by=family]
