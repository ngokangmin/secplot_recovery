
# using new functions from CTFSRPackage
library(date)
# Making the data into the right format # abundFitBT.rdata

splist <- read.delim('/Users/nkmstar/Dropbox/Front Royal 2013/Hokkaido/splist.txt')
splist <- read.delim('/Users/KangMin/Dropbox/Front Royal 2013/Hokkaido/splist.txt')

# secplot checks # don't deal with multiple stems, too complicated
sec1.ori <- biomass.CTFSdb(bukittimahsecondary.stem1, bukittimahsecondary.full1, plot='Bukit Timah Secondar')
sec2.ori <- biomass.CTFSdb(bukittimahsecondary.stem2, bukittimahsecondary.full2, plot='Bukit Timah Secondar')
sec3.ori <- biomass.CTFSdb(bukittimahsecondary.stem3, bukittimahsecondary.full3, plot='Bukit Timah Secondar')
# all have same number of rows
sec3.ori$dbh[sec3.ori$dbh=="0"] <- NA # sec3agb DBH has 0mm 
# dead trees have no dbh # good
# 4 trees recorded dead in census 2 came back to life in census 3 # probably couldn't be found
# species are the same in all 3 censuses # good
sec1.ori[is.na(sec1.ori$date),]$ExactDate <- "2014-08-27"
sec1.ori[is.na(sec1.ori$date),]$date <- 16310 # one date in census 1 is NA
# 1 missing tree in census 2, carried over to census 3, 8 more missing trees in census 3

# remove quadrats K1, K2, K3, N1, O1, P1, L2, M2, N2 # quadrats with >20% primary forest
exclude <- c("K1", "K2", "K3", "N1", "O1", "P1", "L2", "M2", "N2")
sec1 <- subset(sec1.ori, !(quadrat %in% exclude))
sec2 <- subset(sec2.ori, !(quadrat %in% exclude))
sec3 <- subset(sec3.ori, !(quadrat %in% exclude))


# remove unidentified species
sec1.idsp <- subset(sec1, sp != "AAAAAA")
sec2.idsp <- subset(sec2, sp != "AAAAAA")
sec3.idsp <- subset(sec3, sp != "AAAAAA")

# modification from Ryan
# mindbh = 1cm
# function uses status=="A"
bt.modelR12 = model.littleR.Gibbs(cns1=sec1.idsp, cns2=sec2.idsp, modeltype='asymexp', mindbh=10, sptable=splist, start.param=c(-3,.8,.01,.5), bad.modelparam=bad.asymexp.param, steps=1200, burn=200, showstep=50)
bt.modelR23 = model.littleR.Gibbs(cns1=sec2.idsp, cns2=sec3.idsp, modeltype='asymexp', mindbh=10, sptable=splist, start.param=c(-3,.8,.01,.5), bad.modelparam=bad.asymexp.param, steps=1200, burn=200, showstep=50)
bt.modelR13 = model.littleR.Gibbs(cns1=sec1.idsp, cns2=sec3.idsp, modeltype='asymexp', mindbh=10, sptable=splist, start.param=c(-3,.8,.01,.5), bad.modelparam=bad.asymexp.param, steps=1200, burn=200, showstep=50)
# mindbh = 10cm
bt.modelR12.10 = model.littleR.Gibbs(cns1=sec1.idsp, cns2=sec2.idsp, modeltype='asymexp', mindbh=100, sptable=splist, start.param=c(-3,.8,.01,.5), bad.modelparam=bad.asymexp.param, steps=1200, burn=200, showstep=50)
bt.modelR23.10 = model.littleR.Gibbs(cns1=sec2.idsp, cns2=sec3.idsp, modeltype='asymexp', mindbh=100, sptable=splist, start.param=c(-3,.8,.01,.5), bad.modelparam=bad.asymexp.param, steps=1200, burn=200, showstep=50)
bt.modelR13.10 = model.littleR.Gibbs(cns1=sec1.idsp, cns2=sec3.idsp, modeltype='asymexp', mindbh=100, sptable=splist, start.param=c(-3,.8,.01,.5), bad.modelparam=bad.asymexp.param, steps=1200, burn=200, showstep=50)

# export files, so no need to re-run (takes time)
bt.sec.model1cm <- list(census1.2=bt.modelR12, census2.3=bt.modelR23, census1.3=bt.modelR13)
bt.sec.model10cm <- list(census1.2=bt.modelR12.10, census2.3=bt.modelR23.10, census1.3=bt.modelR13.10)
bt.sec.model <- list(dbh10=bt.sec.model1cm, dbh100=bt.sec.model10cm)
save(bt.sec.model, file="/Users/nkmstar/Dropbox/secplot_recovery/abundFitBT_sec_23Nov2016.rdata")

# make new fits file
hyper.means <- as.data.frame(rbind(bt.modelR12$means, bt.modelR23$means, bt.modelR13$means))
CTFSabundmodel_new <- as.data.frame(cbind(hyper.means[,c(1:3)], hyperSDlow=hyper.means$hyperSDlow, hyperSDup=hyper.means$hyperSDup), row.names = c("BT_sec_1.2", "BT_sec_2.3", "BT_sec_1.3"))

hyper.means.10 <- as.data.frame(rbind(bt.modelR12.10$means, bt.modelR23.10$means, bt.modelR13.10$means))
CTFSabundmodel_new.10 <- as.data.frame(cbind(hyper.means.10[,c(1:3)], hyperSDlow=hyper.means.10$hyperSDlow, hyperSDup=hyper.means.10$hyperSDup), row.names = c("BT_sec_1.2", "BT_sec_2.3", "BT_sec_1.3"))

write.table(CTFSabundmodel_new, "/Users/nkmstar/Dropbox/secplot_recovery/CTFS_new_fits_BT.txt", quote=F, sep="\t")
write.table(CTFSabundmodel_new.10, "/Users/nkmstar/Dropbox/secplot_recovery/CTFS_new_fits_BT.10.txt", quote=F, sep="\t")









