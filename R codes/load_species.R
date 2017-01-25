
# all trees alive, dead and recruited, take from census 3
# campau_all <- subset(sec3, sp=="CAMPAU")
# dillsu_all <- subset(sec3, sp=="DILLSU")
# rhodci_all <- subset(sec3, sp=="RHODCI")
# adindu_all <- subset(sec3, sp=="ADINDU")
# streel_all <- subset(sec3, sp=="STREEL")
# prunpo_all <- subset(sec3, sp=="PRUNPO")
# calowa_all <- subset(sec3, sp=="CALOWA")
# elaems_all <- subset(sec3, sp=="ELAEMS")
# fagrfr_all <- subset(sec3, sp=="FAGRFR")
# macaba_all <- subset(sec3, sp=="MACABA")
# artoda_all <- subset(sec3, sp=="ARTODA")

# includes only trees alive in each census

campau.1_alive <- subset(sec1.alive, sp == "CAMPAU")
campau.2_alive <- subset(sec2.alive, sp == "CAMPAU") # includes recruits
campau.3_alive <- subset(sec3.alive, sp == "CAMPAU") # includes recruits

dillsu.1_alive <- subset(sec1.alive, sp == "DILLSU")
dillsu.2_alive <- subset(sec2.alive, sp == "DILLSU")
dillsu.3_alive <- subset(sec3.alive, sp == "DILLSU")

rhodci.1_alive <- subset(sec1.alive, sp == "RHODCI")
rhodci.2_alive <- subset(sec2.alive, sp == "RHODCI")
rhodci.3_alive <- subset(sec3.alive, sp == "RHODCI")

adindu.1_alive <- subset(sec1.alive, sp == "ADINDU")
adindu.2_alive <- subset(sec2.alive, sp == "ADINDU")
adindu.3_alive <- subset(sec3.alive, sp == "ADINDU")

streel.1_alive <- subset(sec1.alive, sp == "STREEL")
streel.2_alive <- subset(sec2.alive, sp == "STREEL")
streel.3_alive <- subset(sec3.alive, sp == "STREEL")

prunpo.1_alive <- subset(sec1.alive, sp == "PRUNPO")
prunpo.2_alive <- subset(sec2.alive, sp == "PRUNPO")
prunpo.3_alive <- subset(sec3.alive, sp == "PRUNPO")

calowa.1_alive <- subset(sec1.alive, sp == "CALOWA")
calowa.2_alive <- subset(sec2.alive, sp == "CALOWA")
calowa.3_alive <- subset(sec3.alive, sp == "CALOWA")

elaems.1_alive <- subset(sec1.alive, sp == "ELAEMS")
elaems.2_alive <- subset(sec2.alive, sp == "ELAEMS")
elaems.3_alive <- subset(sec3.alive, sp == "ELAEMS")

fagrfr.1_alive <- subset(sec1.alive, sp == "FAGRFR")
fagrfr.2_alive <- subset(sec2.alive, sp == "FAGRFR")
fagrfr.3_alive <- subset(sec3.alive, sp == "FAGRFR")

macaba.1_alive <- subset(sec1.alive, sp == "MACABA")
macaba.2_alive <- subset(sec2.alive, sp == "MACABA")
macaba.3_alive <- subset(sec3.alive, sp == "MACABA")

artoda.1_alive <- subset(sec1.alive, sp == "ARTODA")
artoda.2_alive <- subset(sec2.alive, sp == "ARTODA")
artoda.3_alive <- subset(sec3.alive, sp == "ARTODA")


# create objects for dead trees in census 3
campau.3_dead <- subset(sec3, status == "D" & sp == "CAMPAU")
campau.3_surv <- subset(sec3.surv, sp == "CAMPAU")

dillsu.3_dead <- subset(sec3, status == "D" & sp == "DILLSU")
dillsu.3_surv <- subset(sec3.surv, sp == "DILLSU")

rhodci.3_dead <- subset(sec3, status == "D" & sp == "RHODCI")
rhodci.3_surv <- subset(sec3.surv, sp == "RHODCI")

adindu.3_dead <- subset(sec3, status == "D" & sp == "ADINDU")
adindu.3_surv <- subset(sec3.surv, sp == "ADINDU")

# recruits in census 3

streel.3_recr <- subset(sec1, status == "P" & sp == "STREEL")

prunpo.3_recr <- subset(sec1, status == "P" & sp == "PRUNPO")

calowa.3_recr <- subset(sec1, status == "P" & sp == "CALOWA")

elaems.3_recr <- subset(sec1, status == "P" & sp == "ELAEMS")

fagrfr.3_recr <- subset(sec1, status == "P" & sp == "FAGRFR")

macaba.3_recr <- subset(sec1, status == "P" & sp == "MACABA")

artoda.3_recr <- subset(sec1, status == "P" & sp == "ARTODA")


# create point pattern objects using ppp function in spatstat
library(spatstat)

# declining species    # I need to account for existing distributions of trees
campau.1.ppp <- ppp(campau.1_alive$gx, campau.1_alive$gy, c(0,100), c(0,200), marks = campau.1_alive$dbh)
campau.3_dead.ppp <- ppp(campau.3_dead$gx, campau.3_dead$gy, c(0,100), c(0,200))
campau.3_surv.ppp <- ppp(campau.3_surv$gx, campau.3_surv$gy, c(0,100), c(0,200))

dillsu.1.ppp <- ppp(dillsu.1_alive$gx, dillsu.1_alive$gy, c(0,100), c(0,200), marks = dillsu.1_alive$dbh)
dillsu.3_dead.ppp <- ppp(dillsu.3_dead$gx, dillsu.3_dead$gy, c(0,100), c(0,200))
dillsu.3_surv.ppp <- ppp(dillsu.3_surv$gx, dillsu.3_surv$gy, c(0,100), c(0,200))

rhodci.1.ppp <- ppp(rhodci.1_alive$gx, rhodci.1_alive$gy, c(0,100), c(0,200), marks = rhodci.1_alive$dbh)
rhodci.3_dead.ppp <- ppp(rhodci.3_dead$gx, rhodci.3_dead$gy, c(0,100), c(0,200))
rhodci.3_surv.ppp <- ppp(rhodci.3_surv$gx, rhodci.3_surv$gy, c(0,100), c(0,200))

adindu.1.ppp <- ppp(adindu.1_alive$gx, adindu.1_alive$gy, c(0,100), c(0,200), marks = adindu.1_alive$dbh)
adindu.3_dead.ppp <- ppp(adindu.3_dead$gx, adindu.3_dead$gy, c(0,100), c(0,200))
adindu.3_surv.ppp <- ppp(adindu.3_surv$gx, adindu.3_surv$gy, c(0,100), c(0,200))

# recruiting species
streel.1.ppp <- ppp(streel.1_alive$gx, streel.1_alive$gy, c(0,100), c(0,200), marks = streel.1_alive$dbh)
streel.3_recr.ppp <- ppp(streel.3_recr$gx, streel.3_recr$gy, c(0,100), c(0,200))

calowa.1.ppp <- ppp(calowa.1_alive$gx, calowa.1_alive$gy, c(0,100), c(0,200), marks = calowa.1_alive$dbh)
calowa.3_recr.ppp <- ppp(calowa.3_recr$gx, calowa.3_recr$gy, c(0,100), c(0,200))


# duplicate points checks
# duplicated(rhodci.3_surv.ppp)
# duplicated(adindu.3_surv.ppp)
# duplicated(streel.3_recr.ppp)
# duplicate points checks end


# generate 99% envelopes for random mortality hypothesis 

# all surviving trees at census 3
alive_table <- list(campau.3_alive, dillsu.3_alive, rhodci.3_alive, adindu.3_alive)
surv_table <- list(campau.3_surv, dillsu.3_surv, rhodci.3_surv, adindu.3_surv)
pattern_list.2 <- list(campau=NULL, dillsu=NULL, rhodci=NULL, adindu=NULL)
for (i in 1:length(species_table)){
  pattern_list.1 <- list()
  for (j in 1:99){
    temp <- alive_table[[i]][sample(nrow(alive_table[[i]]), nrow(surv_table[[i]]), replace=F),] # random simulation of tree coordinates that survive in census 3
    pattern_list.1[[j]] <- ppp(x=temp$gx, y=temp$gy, c(0,100), c(0,200))
  }
  pattern_list.2[[i]] <- pattern_list.1
}

campau.3_surv.env <- envelope(campau.3_surv.ppp, fun=Kest, simulate=pattern_list.2[[1]])
dillsu.3_surv.env <- envelope(dillsu.3_surv.ppp, fun=Kest, simulate=pattern_list.2[[2]])
rhodci.3_surv.env <- envelope(rhodci.3_surv.ppp, fun=Kest, simulate=pattern_list.2[[3]])
adindu.3_surv.env <- envelope(adindu.3_surv.ppp, fun=Kest, simulate=pattern_list.2[[4]])

# all dead trees at census 3
alive_table <- list(campau.3_alive, dillsu.3_alive, rhodci.3_alive, adindu.3_alive)
dead_table <- list(campau.3_dead, dillsu.3_dead, rhodci.3_dead, adindu.3_dead)
pattern_list.2 <- list(campau=NULL, dillsu=NULL, rhodci=NULL, adindu=NULL)
for (i in 1:length(species_table)){
  pattern_list.1 <- list()
  for (j in 1:99){
    temp <- alive_table[[i]][sample(nrow(alive_table[[i]]), nrow(dead_table[[i]]), replace=F),] # random simulation of tree coordinates that died by census 3
    pattern_list.1[[j]] <- ppp(x=temp$gx, y=temp$gy, c(0,100), c(0,200))
  }
  pattern_list.2[[i]] <- pattern_list.1
}

campau.3_dead.env <- envelope(campau.3_dead.ppp, fun=Kest, simulate=pattern_list.2[[1]])
dillsu.3_dead.env <- envelope(dillsu.3_dead.ppp, fun=Kest, simulate=pattern_list.2[[2]])
rhodci.3_dead.env <- envelope(rhodci.3_dead.ppp, fun=Kest, simulate=pattern_list.2[[3]])
adindu.3_dead.env <- envelope(adindu.3_dead.ppp, fun=Kest, simulate=pattern_list.2[[4]])








