
# Siew Chin's paper used 2004 census

# load both pri and sec plot data first # codes from growth manuscript 

# note that primary forest portions are not included

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

sec1[duplicated(sec1$tag) == T] # no duplicated tags
sec2[duplicated(sec2$tag) == T]
sec3[duplicated(sec3$tag) == T]

# remove quadrats K1, K2, K3, N1, O1, P1, L2, M2, N2, O2 # same as Siew Chin's manuscript
exclude <- c("K1", "K2", "K3", "N1", "O1", "P1", "L2", "M2", "N2", "O2")
sec1 <- subset(sec1.ori, !(quadrat %in% exclude))
sec2 <- subset(sec2.ori, !(quadrat %in% exclude))
sec3 <- subset(sec3.ori, !(quadrat %in% exclude))

# corrections #
sec1[sec1$tag == "K4-01570",]$dbh <- 39.8 # clearly typo error
# move the coordinates of duplicate points by a little, because they cause problems for functions in spatstat
# coordinates of smaller tree is changed by gx-0.1
sec1$gx[sec1$tag=="K5-02015"] <- sec1$gx[sec1$tag=="K5-02015"] - 0.1
sec2$gx[sec2$tag=="K5-02015"] <- sec2$gx[sec2$tag=="K5-02015"] - 0.1
sec3$gx[sec3$tag=="K5-02015"] <- sec3$gx[sec3$tag=="K5-02015"] - 0.1
sec1$gx[sec1$tag=="N3-08642"] <- sec1$gx[sec1$tag=="N3-08642"] - 0.1
sec2$gx[sec2$tag=="N3-08642"] <- sec2$gx[sec2$tag=="N3-08642"] - 0.1
sec3$gx[sec3$tag=="N3-08642"] <- sec3$gx[sec3$tag=="N3-08642"] - 0.1
sec1$gy[sec1$tag=="N3-08645"] <- sec1$gy[sec1$tag=="N3-08645"] - 0.1 # gy changed
sec2$gy[sec2$tag=="N3-08645"] <- sec2$gy[sec2$tag=="N3-08645"] - 0.1 # gy changed
sec3$gy[sec3$tag=="N3-08645"] <- sec3$gy[sec3$tag=="N3-08645"] - 0.1 # gy changed
sec1$gx[sec1$tag=="O3-11011"] <- sec1$gx[sec1$tag=="O3-11011"] - 0.1
sec2$gx[sec2$tag=="O3-11011"] <- sec2$gx[sec2$tag=="O3-11011"] - 0.1
sec3$gx[sec3$tag=="O3-11011"] <- sec3$gx[sec3$tag=="O3-11011"] - 0.1
sec1$gy[sec1$tag=="O3-11010"] <- sec1$gy[sec1$tag=="O3-11010"] - 0.1 # gy changed
sec2$gy[sec2$tag=="O3-11010"] <- sec2$gy[sec2$tag=="O3-11010"] - 0.1 # gy changed
sec3$gy[sec3$tag=="O3-11010"] <- sec3$gy[sec3$tag=="O3-11010"] - 0.1 # gy changed
sec1$gx[sec1$tag=="P2-13119"] <- sec1$gx[sec1$tag=="P2-13119"] - 0.1
sec2$gx[sec2$tag=="P2-13119"] <- sec2$gx[sec2$tag=="P2-13119"] - 0.1
sec3$gx[sec3$tag=="P2-13119"] <- sec3$gx[sec3$tag=="P2-13119"] - 0.1
# LITSMA changed to LITSAC because they are the same
sec1[sec1$sp=="LITSMA",]$sp <- "LITSAC"
sec2[sec2$sp=="LITSMA",]$sp <- "LITSAC"
sec3[sec3$sp=="LITSMA",]$sp <- "LITSAC"
# ARCHCL merge with ARCHMI, all individuals to take code ARCHMI. ARCHCL has incorrect spelling of species name. 
sec1[sec1$sp=="ARCHCL",]$sp <- "ARCHMI"
sec2[sec2$sp=="ARCHCL",]$sp <- "ARCHMI"
sec3[sec3$sp=="ARCHCL",]$sp <- "ARCHMI"
# corrections end

# add basal area column in m^2
sec1$ba <- pi*(sec1$dbh/ 1000 / 2)^2
sec2$ba <- pi*(sec2$dbh/ 1000 / 2)^2
sec3$ba <- pi*(sec3$dbh/ 1000 / 2)^2
# add species classification, sourced from Tree Flora of Malaya, Wayside Trees of Malaya, Flora Malesiana, Plants of Southeast Asia website
sp_class <- read.csv("/Users/nkmstar/Dropbox/secplot_recovery/R codes/sp_abundance_table.csv")
sec1$sp_class <- sp_class[match(sec1$sp, sp_class$sp),]$class
sec2$sp_class <- sp_class[match(sec2$sp, sp_class$sp),]$class
sec3$sp_class <- sp_class[match(sec3$sp, sp_class$sp),]$class
# trees alive
sec1.alive <- subset(sec1, status == "A")
sec2.alive <- subset(sec2, status == "A")
sec3.alive <- subset(sec3, status == "A")
sec1.alive$size <- cut(sec1.alive$dbh, c(10,50,100,10000), include.lowest = T, right = F)
sec2.alive$size <- cut(sec2.alive$dbh, c(10,50,100,10000), include.lowest = T, right = F)
sec3.alive$size <- cut(sec3.alive$dbh, c(10,50,100,10000), include.lowest = T, right = F)
sec2.surv <- subset(sec2, status=="A" & tag %in% subset(sec1, status=="A")$tag) # survivors from census 1
sec3.surv <- subset(sec3, status=="A" & tag %in% subset(sec1, status=="A")$tag) # survivors from census 1

# create a combined file of all censuses with unique quadrat names corresponding to each census
# so that it's easy to create an abundance table for each species and each census
sec.alive <- rbind(sec1.alive, sec2.alive, sec3.alive)
sec.alive$quadrat.census <- c(paste(sec1.alive$quadrat, 1, sep = "."), paste(sec2.alive$quadrat, 2, sep = "."), paste(sec3.alive$quadrat, 3, sep = "."))



# check Siew Chin's files
sc_sec <- read.delim("/Users/nkmstar/Dropbox/secplot_recovery/R codes/s1_1.txt")
sc_sec <- read.delim("/Users/KangMin/Dropbox/secplot_recovery/R codes/s1_1.txt")
sc_sec <- droplevels(subset(sc_sec, !(quadrat %in% exclude)))
nrow(subset(sc_sec, status=="A"))
nrow(subset(sc_sec, status=="A"))/40*25 # this is close to what I got
nrow(sec1.alive)/40*25
nrow(sc_sec)/40*25
table(subset(sc_sec, status=="A")$quadrat)
table(droplevels(subset(sc_sec, status=="A"))$sp) # only MACAHY is missing, a few other species have abundance discrepancies of 1 individual
table(sec1.alive$sp)
sc_sec$ba <- pi*(sc_sec$dbh/ 1000 / 2)^2
sum(sc_sec$ba, na.rm=T) / 40 * 25
sum(sec1.alive$ba, na.rm=T) / 40 * 25 # discrepancy
temp2 <- sec1.alive[FALSE,]
for (i in 1:nrow(sc_sec)){
  temp1 <- sec1.alive[match(sc_sec$tag[i], sec1.alive$tag),]
  if (!is.na(temp1$dbh)){
    if (sc_sec$dbh[i] == temp1$dbh) temp2 <- temp2
    if (sc_sec$dbh[i] != temp1$dbh) temp2 <- rbind(temp2, temp1)
  }
}
subset(sc_sec, tag=="K4-01570") # this is the problem, probably the problem disappeared as data from the later two censuses were used for checking
sort(tapply(sc_sec$ba, sc_sec$sp, sum, na.rm=T)) # using her file I'm getting exactly the same results as I did with my copy, so I'm not sure how she got the results in the supplementary materials she published in the 2013 paper





# look at degree of similarity between censuses
# use hierarchical clustering to group all quadrats from all censuses

library(vegan)

sec.spec <- tapply(sec.alive$sp, list(sec.alive$quadrat.census, sec.alive$sp), length)
sec.spec[is.na(sec.spec)] <- 0
sec.dist <- vegdist(sec.spec, method = "bray")
sec.wardD <- hclust(sec.dist, "ward.D")
sec.wardD2 <- hclust(sec.dist, "ward.D2")

plot(sec.wardD)
plot(sec.wardD2) # both methods show little plot-wide change in species composition, but census 3 is constantly more different than both census 1 and 2. 

# use NMDS to look at distances between quadrats from all censuses
library(MASS)
isoMDS(sec.dist) # error

# ANOSIM # may not be reliable
census <- rep(c("census1", "census2", "census3"), times=41)
sec.ano <- anosim(sec.dist, census)
summary(sec.ano)
plot(sec.ano)

# ADONIS # assumes that sites are interchangeable under a true null hypothesis, but mine are temporally correlated
adonis(sec.spec ~ census, permutations=100) # no significant differences between censuses
adonis(sec.dist ~ census, permutations=100)

# look at degree of similarity between censuses END


# species abundances
sort(table(sec1.alive$sp)); length(table(sec1.alive$sp))-2 # 125 species; minus AAAAAA and PANDDP
sort(table(sec2.alive$sp)); length(table(sec2.alive$sp))-2 # 139 species (+14)
sort(table(sec3.alive$sp)); length(table(sec3.alive$sp))-2 # 148 species (+11)
# number of species has been increasing at each census
# number of species in each species class
table(as.data.frame(cbind(sp=unique(sec1.alive$sp), sp_class=as.character(sp_class[match(unique(sec1.alive$sp), sp_class$sp),"class"])))$sp_class)
table(as.data.frame(cbind(sp=unique(sec2.alive$sp), sp_class=as.character(sp_class[match(unique(sec2.alive$sp), sp_class$sp),"class"])))$sp_class)
table(as.data.frame(cbind(sp=unique(sec3.alive$sp), sp_class=as.character(sp_class[match(unique(sec3.alive$sp), sp_class$sp),"class"])))$sp_class)
# mean species richness in each quadrat
no.sp.quad1 <- numeric()
for (i in 1:40){
  quad <- levels(factor(sec1.alive$quadrat))[i]
  temp1 <- subset(sec1.alive, quadrat == quad & sp !="AAAAAA")
  no.sp.quad1[i] <- length(table(temp1$sp))
}
# poster figure: no. species in each species class in 2004, 2008 and 2012 
sp_class_sp <- cbind( table(as.data.frame(cbind(sp=unique(sec1.alive$sp), sp_class=as.character(sp_class[match(unique(sec1.alive$sp), sp_class$sp),"class"])))$sp_class), table(as.data.frame(cbind(sp=unique(sec2.alive$sp), sp_class=as.character(sp_class[match(unique(sec2.alive$sp), sp_class$sp),"class"])))$sp_class), table(as.data.frame(cbind(sp=unique(sec3.alive$sp), sp_class=as.character(sp_class[match(unique(sec3.alive$sp), sp_class$sp),"class"])))$sp_class)  )
dimnames(sp_class_sp)[[2]] <- c("2004", "2008", "2012")
png("/Users/nkmstar/Dropbox/secplot_recovery/ESJ_2017_poster/sp_class_sp.png", width=650, height=600)
par(mar=c(5.1, 4.5, 4.1, 2.1))
barplot(sp_class_sp, ylim=c(0,200), ylab="No. species", main="No. species", legend=c("exotic species", "primary forest species", "primary + secondary forest species", "secondary forest species"), args.legend = list(x=3, y=200, bty="n", y.intersp = 0.8, cex=2), cex.axis = 2, cex.names = 2, cex.lab=2, cex.main=2)
dev.off()
# combine these into one table # refer to sp_abundance_table.xlsx


# stem density
# whole plot
nrow(sec1.alive)
nrow(sec2.alive)
nrow(sec3.alive) # no major changes in stem density, which means that losses in species are balanced by gains in other species
# per ha
nrow(sec1.alive)/40*25
nrow(sec2.alive)/40*25
nrow(sec3.alive)/40*25
# 99% CI by bootstrapping quadrats
ind1.quad <- tapply(sec1.alive$sp, sec1.alive$quadrat, length) * 25
ind2.quad <- tapply(sec2.alive$sp, sec2.alive$quadrat, length) * 25
ind3.quad <- tapply(sec3.alive$sp, sec3.alive$quadrat, length) * 25
ind1.sim <- numeric()
ind2.sim <- numeric()
ind3.sim <- numeric()
for (i in 1:1000){
  temp <- sample(1:40, 40, replace=T)
  ind1.sim[i] <- mean(ind1.quad[temp])
  ind2.sim[i] <- mean(ind2.quad[temp])
  ind3.sim[i] <- mean(ind3.quad[temp])
}
quantile(ind1.sim, c(0.005, 0.995)); mean(ind1.sim)
quantile(ind2.sim, c(0.005, 0.995)); mean(ind2.sim)
quantile(ind3.sim, c(0.005, 0.995)); mean(ind3.sim)
# by size class per ha
table(sec1.alive$size)/40*25
table(sec2.alive$size)/40*25
table(sec3.alive$size)/40*25
# by species class per ha
table(sec1.alive$sp_class)/40*25
table(sec2.alive$sp_class)/40*25
table(sec3.alive$sp_class)/40*25
# poster figure: no. individuals in each species class in 2004, 2008 and 2012 
sp_class_ind <- cbind(table(sec1.alive$sp_class)/40*25, table(sec2.alive$sp_class)/40*25, table(sec3.alive$sp_class)/40*25 )
dimnames(sp_class_ind)[[2]] <- c("2004", "2008", "2012")
png("/Users/nkmstar/Dropbox/secplot_recovery/ESJ_2017_poster/sp_class_ind.png", width=650, height=600)
par(mar=c(5.1, 4.5, 4.1, 2.1))
barplot(sp_class_ind, ylim=c(0,2500), ylab="No. individuals", main="No. individuals", legend=c("exotic species", "primary forest species", "primary + secondary forest species", "secondary forest species"), args.legend = list(x=3, y=2600, bty="n", y.intersp = 0.8, cex=2), cex.axis = 2, cex.names = 2, cex.lab=2, cex.main=2)
dev.off()
# stem density END

# basal area per ha
sum(sec1.alive$ba, na.rm=T) / 40 * 25
sum(sec2.alive$ba, na.rm=T) / 40 * 25
sum(sec3.alive$ba, na.rm=T) / 40 * 25
# 99% CI by bootstrapping quadrats
ba1.quad <- tapply(sec1.alive$ba, sec1.alive$quadrat, sum, na.rm=T) * 25
ba2.quad <- tapply(sec2.alive$ba, sec2.alive$quadrat, sum, na.rm=T) * 25
ba3.quad <- tapply(sec3.alive$ba, sec3.alive$quadrat, sum, na.rm=T) * 25
ba1.sim <- numeric()
ba2.sim <- numeric()
ba3.sim <- numeric()
for (i in 1:1000){
  temp <- sample(1:40, 40, replace=T)
  ba1.sim[i] <- mean(ba1.quad[temp])
  ba2.sim[i] <- mean(ba2.quad[temp])
  ba3.sim[i] <- mean(ba3.quad[temp])
}
quantile(ba1.sim, c(0.005, 0.995)); mean(ba1.sim)
quantile(ba2.sim, c(0.005, 0.995)); mean(ba2.sim)
quantile(ba3.sim, c(0.005, 0.995)); mean(ba3.sim)
# basal area by species per ha
sort(tapply(sec1.alive$ba, sec1.alive$sp, sum, na.rm=T)) / 40 * 25
sort(tapply(sec2.alive$ba, sec2.alive$sp, sum, na.rm=T)) / 40 * 25
sort(tapply(sec3.alive$ba, sec3.alive$sp, sum, na.rm=T)) / 40 * 25
# basal area by species class
tapply(sec1.alive$ba, sec1.alive$sp_class, sum, na.rm=T) / 40 * 25
tapply(sec2.alive$ba, sec2.alive$sp_class, sum, na.rm=T) / 40 * 25
tapply(sec3.alive$ba, sec3.alive$sp_class, sum, na.rm=T) / 40 * 25
# poster figure: basal area in each species class in 2004, 2008 and 2012 
sp_class_ba <- cbind(tapply(sec1.alive$ba, sec1.alive$sp_class, sum, na.rm=T) / 40 * 25, tapply(sec2.alive$ba, sec2.alive$sp_class, sum, na.rm=T) / 40 * 25, tapply(sec3.alive$ba, sec3.alive$sp_class, sum, na.rm=T) / 40 * 25 )
dimnames(sp_class_ba)[[2]] <- c("2004", "2008", "2012")
png("/Users/nkmstar/Dropbox/secplot_recovery/ESJ_2017_poster/sp_class_ba.png", width=650, height=600)
par(mar=c(5.1, 5.4, 4.1, 2.1))
barplot(sp_class_ba, ylim=c(0,30), ylab=expression(paste("Basal area", (m^2/ha))), main="Basal area", legend=c("exotic species", "primary forest species", "primary + secondary forest species", "secondary forest species"), args.legend = list(x=3, y=30, bty="n", y.intersp = 0.8, cex=2), cex.axis = 2, cex.names = 2, cex.lab=2, cex.main=2)
dev.off()
# basal area END

# mortality and recruitment
mortality(sec1, sec2)
mortality(sec2, sec3)
recruitment(sec1, sec2)
recruitment(sec2, sec3)
# 99% CI by bootstrapping quadrats
mort12.quad <- mortality(sec1, sec2, split1 = sec1$quadrat)$rate
mort23.quad <- mortality(sec2, sec3, split1 = sec2$quadrat)$rate
recr12.quad <- recruitment(sec1, sec2, split1 = sec2$quadrat)$rate
recr23.quad <- recruitment(sec2, sec3, split1 = sec3$quadrat)$rate
mort12.sim <- numeric()
mort23.sim <- numeric()
recr12.sim <- numeric()
recr23.sim <- numeric()
for (i in 1:1000){
  temp <- sample(1:40, 40, replace=T)
  mort12.sim[i] <- mean(mort12.quad[temp])
  mort23.sim[i] <- mean(mort23.quad[temp])
  recr12.sim[i] <- mean(recr12.quad[temp])
  recr23.sim[i] <- mean(recr23.quad[temp])
}
quantile(mort12.sim, c(0.005, 0.995))
quantile(mort23.sim, c(0.005, 0.995))
quantile(recr12.sim, c(0.005, 0.995))
quantile(recr23.sim, c(0.005, 0.995))
# by species
mortality(sec1, sec2, split1 = sec1$sp)$rate
mortality(sec2, sec3, split1 = sec2$sp)$rate
recruitment(sec1, sec2, split1 = sec1$sp)$rate
recruitment(sec2, sec3, split1 = sec2$sp)$rate
# by species class
mortality(sec1, sec2, split1 = sec1$sp_class)$rate
mortality(sec2, sec3, split1 = sec2$sp_class)$rate
recruitment(sec1, sec2, split1 = sec1$sp_class)$rate
recruitment(sec2, sec3, split1 = sec2$sp_class)$rate
# mortality and recruitment END


# dominant species
splist <- read.delim("/Users/nkmstar/Dropbox/Front Royal 2013/Hokkaido/splist.txt")
splist <- read.delim("/Users/KangMin/Dropbox/Front Royal 2013/Hokkaido/splist.txt")
sort(tapply(sec1$ba, sec1$sp, sum, na.rm=T), decreasing=T)[1:10]
top10_2004 <- names(sort(tapply(sec1$ba, sec1$sp, sum, na.rm=T), decreasing=T)[1:10])
table(subset(sec1.alive, sp %in% top10_2004)$sp)
tapply(subset(sec3, sp %in% top10_2004)$ba, subset(sec3, sp %in% top10_2004)$sp, sum, na.rm=T)
table(subset(sec3.alive, sp %in% top10_2004)$sp)
mortality(sec1, sec3, split1=sec1$sp)$rate[names(mortality(sec1, sec3, split1=sec1$sp)$rate) %in% top10_2004]
recruitment(sec1, sec3, split1=sec1$sp)$rate[names(recruitment(sec1, sec3, split1=sec1$sp)$rate) %in% top10_2004]
# plot of 4 dominant species abundance in census 1 and 3
four_dom_sp <- rbind(subset(sec1.alive, sp %in% c("CAMPAU", "DILLSU", "ADINDU", "RHODCI")), subset(sec3.alive, sp %in% c("CAMPAU", "DILLSU", "ADINDU", "RHODCI")))
four_dom_sp_mat1 <- tapply(four_dom_sp$size, list(four_dom_sp$size, four_dom_sp$CensusID, four_dom_sp$sp), length)
adindu_mat <- matrix(four_dom_sp_mat1[1:6], nrow=3, byrow=F, dimnames=list(c("[10,50)", "[50,100)", "[100,1e+04]"), c("2004", "2012")))
campau_mat <- matrix(four_dom_sp_mat1[7:12], nrow=3, byrow=F, dimnames=list(c("[10,50)", "[50,100)", "[100,1e+04]"), c("2004", "2012")))
dillsu_mat <- matrix(four_dom_sp_mat1[13:18], nrow=3, byrow=F, dimnames=list(c("[10,50)", "[50,100)", "[100,1e+04]"), c("2004", "2012")))
rhodci_mat <- matrix(four_dom_sp_mat1[19:24], nrow=3, byrow=F, dimnames=list(c("[10,50)", "[50,100)", "[100,1e+04]"), c("2004", "2012")))
png("/Users/KangMin/Dropbox/secplot_recovery/figures/dom_sp_abun.png", width=835, height=560)
par(mfrow=c(2,2))
barplot(adindu_mat, beside=T, ylab="No. individuals", ylim=c(0,400), main=expression(italic("Adinandra dumosa")), cex.axis = 1.4, cex.names = 1.4, cex.lab=1.4, cex.main=1.6)
barplot(campau_mat, beside=T, ylab="No. individuals", ylim=c(0,400), main=expression(italic("Campnosperma auriculata")), cex.axis = 1.4, cex.names = 1.4, cex.lab=1.4, cex.main=1.6)
barplot(dillsu_mat, beside=T, ylab="No. individuals", ylim=c(0,400), main=expression(italic("Dillenia suffruticosa")), cex.axis = 1.4, cex.names = 1.4, cex.lab=1.4, cex.main=1.6)
barplot(rhodci_mat, beside=T, ylab="No. individuals", ylim=c(0,400), main=expression(italic("Rhodamnia cinerea")), legend.text = c("1.0-4.9 cm", "5.0-9.9 cm", expression("">="10 cm")), args.legend = list(x=4, y=400, bty="n", y.intersp=1.4, cex=1.4), cex.axis = 1.4, cex.names = 1.4, cex.lab=1.4, cex.main=1.6)
dev.off()
# figure for poster 
png("/Users/nkmstar/Dropbox/secplot_recovery/ESJ_2017_poster/dom_sp_abun.png", width=1000, height=700)
par(mfrow=c(2,2), mar=c(5.1, 4.5, 4.1, 2.1))
barplot(adindu_mat, beside=T, ylab="No. individuals", ylim=c(0,400), main=expression(italic("Adinandra dumosa")), cex.axis = 2, cex.names = 2, cex.lab=2, cex.main=2.4)
barplot(campau_mat, beside=T, ylab="No. individuals", ylim=c(0,400), main=expression(italic("Campnosperma auriculata")), cex.axis = 2, cex.names = 2, cex.lab=2, cex.main=2.4)
barplot(dillsu_mat, beside=T, ylab="No. individuals", ylim=c(0,400), main=expression(italic("Dillenia suffruticosa")), cex.axis = 2, cex.names = 2, cex.lab=2, cex.main=2.4)
barplot(rhodci_mat, beside=T, ylab="No. individuals", ylim=c(0,400), main=expression(italic("Rhodamnia cinerea")), legend.text = c("1.0-4.9 cm", "5.0-9.9 cm", expression("">="10 cm")), args.legend = list(x=4, y=400, bty="n", y.intersp=2, cex=2), cex.axis = 2, cex.names = 2, cex.lab=2, cex.main=2.4)
dev.off()
# dominant species END



# number of species in each size class
table(subset(sec1, size=="[10,20)")$sp); length(table(subset(sec1, size=="[10,20)")$sp))
table(subset(sec1, size=="[20,100)")$sp); length(table(subset(sec1, size=="[20,100)")$sp))-1
table(subset(sec1, size=="[100,1e+04]")$sp); length(table(subset(sec1, size=="[100,1e+04]")$sp))-1
table(subset(sec2, size=="[10,20)")$sp); length(table(subset(sec2, size=="[10,20)")$sp))
table(subset(sec2, size=="[20,100)")$sp); length(table(subset(sec2, size=="[20,100)")$sp))-1
table(subset(sec2, size=="[100,1e+04]")$sp); length(table(subset(sec2, size=="[100,1e+04]")$sp))-1
table(subset(sec3, size=="[10,20)")$sp); length(table(subset(sec3, size=="[10,20)")$sp))
table(subset(sec3, size=="[20,100)")$sp); length(table(subset(sec3, size=="[20,100)")$sp))-1
table(subset(sec3, size=="[100,1e+04]")$sp); length(table(subset(sec3, size=="[100,1e+04]")$sp))-1






# can two spatially random patterns be correlated? 
library(spatstat)
x <- runif(400, 0, 100)
y <- runif(400, 0, 200)
m <- sample(1:2, 400, replace=T)
m <- factor(m, levels=c(1,2))
Ran <- ppp(x, y, c(0,100), c(0,200), marks=m)
plot(Ran)
plot(Kcross(Ran))
plot(Lcross(Ran))



# plot spatial distributions of species in each census
source("/Users/nkmstar/Dropbox/secplot_recovery/R codes/load_species.R")
source("/Users/KangMin/Dropbox/secplot_recovery/R codes/load_species.R")

png("/Users/nkmstar/Desktop/CAMPAU.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(campau.1_alive$gx, campau.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=campau.1_alive$dbh/100, asp=1, main= expression(italic("Campnosperma auriculata")))
points(campau.2_alive$gx, campau.2_alive$gy, cex=campau.2_alive$dbh/150, asp=1, pch=16, col="red")
points(campau.3_alive$gx, campau.3_alive$gy, cex=campau.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=0, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.7)
dev.off()

## for making gif
png("/Users/nkmstar/Desktop/CAMPAU_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(campau.1_alive$gx, campau.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=campau.1_alive$dbh/100, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/CAMPAU_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(campau.2_alive$gx, campau.2_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=campau.2_alive$dbh/100, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/CAMPAU_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(campau.3_alive$gx, campau.3_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=campau.3_alive$dbh/100, asp=1, main="2012")
dev.off()
## gif end

par(xaxs="i", yaxs="i")
plot(dillsu.1_alive$gx, dillsu.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=dillsu.1_alive$dbh/100, asp=1, main= expression(italic("Dillenia suffruticosa")))
points(dillsu.2_alive$gx, dillsu.2_alive$gy, cex=dillsu.2_alive$dbh/150, asp=1, pch=16, col="red")
points(dillsu.3_alive$gx, dillsu.3_alive$gy, cex=dillsu.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

# for making gif
png("/Users/nkmstar/Desktop/DILLSU_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(dillsu.1_alive$gx, dillsu.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=dillsu.1_alive$dbh/100, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/DILLSU_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(dillsu.2_alive$gx, dillsu.2_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=dillsu.2_alive$dbh/100, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/DILLSU_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(dillsu.3_alive$gx, dillsu.3_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=dillsu.3_alive$dbh/100, asp=1, main="2012")
dev.off()
# gif end

par(xaxs="i", yaxs="i")
plot(adindu.1_alive$gx, adindu.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=adindu.1_alive$dbh/100, asp=1, main= expression(italic("Adinandra dumosa")))
points(adindu.2_alive$gx, adindu.2_alive$gy, cex=adindu.2_alive$dbh/150, asp=1, pch=16, col="red")
points(adindu.3_alive$gx, adindu.3_alive$gy, cex=adindu.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

# for making gif
png("/Users/nkmstar/Desktop/ADINDU_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(adindu.1_alive$gx, adindu.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=adindu.1_alive$dbh/100, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/ADINDU_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(adindu.2_alive$gx, adindu.2_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=adindu.2_alive$dbh/100, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/ADINDU_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(adindu.3_alive$gx, adindu.3_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=adindu.3_alive$dbh/100, asp=1, main="2012")
dev.off()
# gif end

par(xaxs="i", yaxs="i")
plot(streel.1_alive$gx, streel.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=streel.1_alive$dbh/100, asp=1, main= expression(italic("Streblus elongatus")))
points(streel.2_alive$gx, streel.2_alive$gy, cex=streel.2_alive$dbh/150, asp=1, pch=16, col="red")
points(streel.3_alive$gx, streel.3_alive$gy, cex=streel.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

# for making gif # NOTE THAT SCALE OF DOTS ARE NOT THE SAME AS CAMPAU, DILLSU & ADINDU
png("/Users/nkmstar/Desktop/STREEL_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(streel.1_alive$gx, streel.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=streel.1_alive$dbh/70, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/STREEL_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(streel.2_alive$gx, streel.2_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=streel.2_alive$dbh/70, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/STREEL_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(streel.3_alive$gx, streel.3_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=streel.3_alive$dbh/70, asp=1, main="2012")
dev.off()
# gif end

par(xaxs="i", yaxs="i")
plot(prunpo.1_alive$gx, prunpo.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=prunpo.1_alive$dbh/100, asp=1, main= expression(italic("Prunus polystachyus")))
points(prunpo.2_alive$gx, prunpo.2_alive$gy, cex=prunpo.2_alive$dbh/150, asp=1, pch=16, col="red")
points(prunpo.3_alive$gx, prunpo.3_alive$gy, cex=prunpo.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(calowa.1_alive$gx, calowa.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=calowa.1_alive$dbh/100, asp=1, main= expression(italic("Calophyllum wallichianum")))
points(calowa.2_alive$gx, calowa.2_alive$gy, cex=calowa.2_alive$dbh/150, asp=1, pch=16, col="red")
points(calowa.3_alive$gx, calowa.3_alive$gy, cex=calowa.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

# for making gif # NOTE THAT SCALE OF DOTS ARE NOT THE SAME AS CAMPAU, DILLSU & ADINDU
png("/Users/nkmstar/Desktop/CALOWA_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(calowa.1_alive$gx, calowa.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=calowa.1_alive$dbh/70, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/CALOWA_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(calowa.2_alive$gx, calowa.2_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=calowa.2_alive$dbh/70, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/CALOWA_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(calowa.3_alive$gx, calowa.3_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=calowa.3_alive$dbh/70, asp=1, main="2012")
dev.off()
# gif end

par(xaxs="i", yaxs="i")
plot(elaems.1_alive$gx, elaems.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=elaems.1_alive$dbh/100, asp=1, main= expression(italic("Elaeocarpus mastersii")))
points(elaems.2_alive$gx, elaems.2_alive$gy, cex=elaems.2_alive$dbh/150, asp=1, pch=16, col="red")
points(elaems.3_alive$gx, elaems.3_alive$gy, cex=elaems.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

# for making gif # NOTE THAT SCALE OF DOTS ARE NOT THE SAME AS CAMPAU, DILLSU & ADINDU
png("/Users/nkmstar/Desktop/ELAEMS_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(elaems.1_alive$gx, elaems.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=elaems.1_alive$dbh/70, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/ELAEMS_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(elaems.2_alive$gx, elaems.2_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=elaems.2_alive$dbh/70, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/ELAEMS_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(elaems.3_alive$gx, elaems.3_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=elaems.3_alive$dbh/70, asp=1, main="2012")
dev.off()
# gif end

par(xaxs="i", yaxs="i")
plot(fagrfr.1_alive$gx, fagrfr.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=fagrfr.1_alive$dbh/100, asp=1, main= expression(italic("Fagraea fragrans")))
points(fagrfr.2_alive$gx, fagrfr.2_alive$gy, cex=fagrfr.2_alive$dbh/150, asp=1, pch=16, col="red")
points(fagrfr.3_alive$gx, fagrfr.3_alive$gy, cex=fagrfr.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(macaba.1_alive$gx, macaba.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=macaba.1_alive$dbh/100, asp=1, main= expression(italic("Macaranga bancana")))
points(macaba.2_alive$gx, macaba.2_alive$gy, cex=macaba.2_alive$dbh/150, asp=1, pch=16, col="red")
points(macaba.3_alive$gx, macaba.3_alive$gy, cex=macaba.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(artoda.1_alive$gx, artoda.1_alive$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=artoda.1_alive$dbh/100, asp=1, main= expression(italic("Artocarpus dadah")))
points(artoda.2_alive$gx, artoda.2_alive$gy, cex=artoda.2_alive$dbh/150, asp=1, pch=16, col="red")
points(artoda.3_alive$gx, artoda.3_alive$gy, cex=artoda.3_alive$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

plot(x=1:3, y=c(nrow(campau.1_alive), nrow(campau.2_alive), nrow(campau.3_alive)), type="l", xlab="", ylab="No. individuals")

# tried to insert inset graphs but failed
library(TeachingDemos)
subplot(plot(x=1:3, y=c(nrow(campau.1_alive), nrow(campau.2_alive), nrow(campau.3_alive)), type="l", xlab="", ylab="No. individuals"), x=0, y=50, size=c(2,0.5))



# spatial distribution of species and their mortality/recruitment (taking 8-year mortality/recruitment 2004-2012)

# NULL MODEL: random mortality hypothesis # Kenkel 1988
# using all trees across different times, and randomly select n_dead_trees from trees alive at census 1 for each species

# plotting observed vs null hypothesis of complete spatial randomness

# 99% CSR envelopes for all trees alive at census 1 # because 99 Monte Carlo simulations were used # Haase 1995, Getzin et al. 2006

# plotting initial and survived trees together
plot(campau.1.env <- envelope(campau.1.ppp, fun=Kest), sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
lines(x=campau.3_surv.env$r, y=sqrt(campau.3_surv.env$obs/pi), lty=3)
legend(x=0, y=40, legend=c("initial", "survivors"), lty=c(1,3), bty="n", y.intersp = 0.4)
plot(dillsu.1.env <- envelope(dillsu.1.ppp, fun=Kest), sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
lines(x=dillsu.3_surv.env$r, y=sqrt(dillsu.3_surv.env$obs/pi), lty=3)
plot(rhodci.1.env <- envelope(rhodci.1.ppp, fun=Kest), sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
lines(x=rhodci.3_surv.env$r, y=sqrt(rhodci.3_surv.env$obs/pi), lty=3)
plot(adindu.1.env <- envelope(adindu.1.ppp, fun=Kest), sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
lines(x=adindu.3_surv.env$r, y=sqrt(adindu.3_surv.env$obs/pi), lty=3)

# 99% envelopes for random mortality hypothesis 

# all surviving trees at census 3 # are surviving trees more regularly spread than expected by random mortality? 
par(mfrow=c(2,2))
plot(campau.3_surv.env, sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
text(x=c(1,1), y=c(38,33), labels=c(expression("N"[initial]), expression("N"[survived])), adj=c(0,NA))
text(x=c(5,5), y=c(38,33), labels=c("=203", "=162"), adj=c(0,NA))
plot(dillsu.3_surv.env, sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
text(x=c(1,1), y=c(38,33), labels=c(expression("N"[initial]), expression("N"[survived])), adj=c(0,NA))
text(x=c(5,5), y=c(38,33), labels=c("=780", "=521"), adj=c(0,NA))
plot(rhodci.3_surv.env, sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
text(x=c(1,1), y=c(38,33), labels=c(expression("N"[initial]), expression("N"[survived])), adj=c(0,NA))
text(x=c(5,5), y=c(38,33), labels=c("=149", "=141"), adj=c(0,NA))
plot(adindu.3_surv.env, sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
text(x=c(1,1), y=c(38,33), labels=c(expression("N"[initial]), expression("N"[survived])), adj=c(0,NA))
text(x=c(5,5), y=c(38,33), labels=c("=471", "=375"), adj=c(0,NA))

# all trees that died by census 3 # are dead trees randomly distributed? 
par(mfrow=c(2,2))
plot(campau.3_dead.env, sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
text(x=c(1,1), y=c(38,33), labels=c(expression("N"[initial]), expression("N"[dead])), adj=c(0,NA))
text(x=c(5,5), y=c(38,33), labels=c("=203", "=41"), adj=c(0,NA))
plot(dillsu.3_dead.env, sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
text(x=c(1,1), y=c(38,33), labels=c(expression("N"[initial]), expression("N"[dead])), adj=c(0,NA))
text(x=c(5,5), y=c(38,33), labels=c("=780", "=280"), adj=c(0,NA))
plot(rhodci.3_dead.env, sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
text(x=c(1,1), y=c(38,33), labels=c(expression("N"[initial]), expression("N"[dead])), adj=c(0,NA))
text(x=c(5,5), y=c(38,33), labels=c("=149", "=9"), adj=c(0,NA))
plot(adindu.3_dead.env, sqrt(./pi) ~ r, ylim=c(0,40), legend=F)
text(x=c(1,1), y=c(38,33), labels=c(expression("N"[initial]), expression("N"[dead])), adj=c(0,NA))
text(x=c(5,5), y=c(38,33), labels=c("=471", "=96"), adj=c(0,NA))









# look at effect of weather 








