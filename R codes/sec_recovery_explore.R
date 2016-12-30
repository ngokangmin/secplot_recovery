
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

# remove quadrats K1, K2, K3, N1, O1, P1, L2, M2, N2 # quadrats with >20% primary forest
exclude <- c("K1", "K2", "K3", "N1", "O1", "P1", "L2", "M2", "N2")
sec1 <- subset(sec1.ori, !(quadrat %in% exclude))
sec2 <- subset(sec2.ori, !(quadrat %in% exclude))
sec3 <- subset(sec3.ori, !(quadrat %in% exclude))

# corrections #
sec1[sec1$tag == "K4-01570",]$dbh <- 39.8 # clearly typo error
# corrections end

sec1.alive <- subset(sec1, status == "A")
sec2.alive <- subset(sec2, status == "A")
sec3.alive <- subset(sec3, status == "A")

# create a combined file of all censuses with unique quadrat names corresponding to each census
# so that it's easy to create an abundance table for each species and each census
sec.alive <- rbind(sec1.alive, sec2.alive, sec3.alive)
sec.alive$quadrat.census <- c(paste(sec1.alive$quadrat, 1, sep = "."), paste(sec2.alive$quadrat, 2, sep = "."), paste(sec3.alive$quadrat, 3, sep = "."))


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


# species abundances
sort(table(sec1.alive$sp)); length(table(sec1.alive$sp))-1 # 126 species
sort(table(sec2.alive$sp)); length(table(sec2.alive$sp))-1 # 140 species (+14)
sort(table(sec3.alive$sp)); length(table(sec3.alive$sp))-1 # 151 species (+11)
# number of species has been increasing at each census
# combine these into one table # refer to sp_abundance_table.xlsx


# stem density
# whole plot
nrow(sec1.alive)
nrow(sec2.alive)
nrow(sec3.alive) # no major changes in stem density, which means that losses in species are balanced by gains in other species


# basal area
sum(pi * (sec1.alive$dbh / 1000 / 2)^2, na.rm=T) / 2
sum(pi * (sec2.alive$dbh / 1000 / 2)^2, na.rm=T) / 2
sum(pi * (sec3.alive$dbh / 1000 / 2)^2, na.rm=T) / 2


# density of trees >= 30 cm dbh in priplot
source("/Users/KangMin/Dropbox/dendro Data analysis/R_scripts/load_main_census_data.R")
source("/Users/nkmstar/Dropbox/dendro Data analysis/R_scripts/load_main_census_data.R")
# use all size classes (because in load_main_census_data.R sizes were calculated for >=5 cm only)
pri1$size <- cut(pri1$dbh, c(10, 100, 200, 400, 10000), include.lowest=T, right=F)
pri2$size <- cut(pri2$dbh, c(10, 100, 200, 400, 10000), include.lowest=T, right=F)
pri4$size <- cut(pri4$dbh, c(10, 100, 200, 400, 10000), include.lowest=T, right=F)
pri5$size <- cut(pri5$dbh, c(10, 100, 200, 400, 10000), include.lowest=T, right=F)
pri6$size <- cut(pri6$dbh, c(10, 100, 200, 400, 10000), include.lowest=T, right=F)
sec1$size <- cut(sec1$dbh, c(10, 100, 200, 400, 10000), include.lowest=T, right=F)
sec2$size <- cut(sec2$dbh, c(10, 100, 200, 400, 10000), include.lowest=T, right=F)
sec3$size <- cut(sec3$dbh, c(10, 100, 200, 400, 10000), include.lowest=T, right=F)
nrow(subset(pri4, dbh>=300))/2
nrow(subset(pri5, dbh>=300))/2
nrow(subset(pri6, dbh>=300))/2
sort(table(sec1.alive$sp))
sort(table(sec2.alive$sp))
sort(table(sec3.alive$sp))
# CAMPAU cluster
campau.quad <- c("L1", "M1", "K4", "K5", "L3", "L4", "L5", "M3", "M4", "M5", "N3", "N4", "N5", "O2", "O3", "O4", "O5")
nrow(subset(sec1, quadrat %in% campau.quad & dbh>=300))/17*25
nrow(subset(sec2, quadrat %in% campau.quad & dbh>=300))/17*25
nrow(subset(sec3, quadrat %in% campau.quad & dbh>=300))/17*25



# plot spatial distributions of species in each census
source("/Users/nkmstar/Dropbox/secplot_recovery/R codes/load_species.R")
source("/Users/KangMin/Dropbox/secplot_recovery/R codes/load_species.R")

# plot dimensions 

png("/Users/nkmstar/Desktop/CAMPAU.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(campau.1$gx, campau.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=campau.1$dbh/100, asp=1, main= expression(italic("Campnosperma auriculata")))
points(campau.2$gx, campau.2$gy, cex=campau.2$dbh/150, asp=1, pch=16, col="red")
points(campau.3$gx, campau.3$gy, cex=campau.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=0, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.7)
dev.off()

## for making gif
png("/Users/nkmstar/Desktop/CAMPAU_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(campau.1$gx, campau.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=campau.1$dbh/100, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/CAMPAU_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(campau.2$gx, campau.2$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=campau.2$dbh/100, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/CAMPAU_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(campau.3$gx, campau.3$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=campau.3$dbh/100, asp=1, main="2012")
dev.off()
## gif end

par(xaxs="i", yaxs="i")
plot(dillsu.1$gx, dillsu.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=dillsu.1$dbh/100, asp=1, main= expression(italic("Dillenia suffruticosa")))
points(dillsu.2$gx, dillsu.2$gy, cex=dillsu.2$dbh/150, asp=1, pch=16, col="red")
points(dillsu.3$gx, dillsu.3$gy, cex=dillsu.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

# for making gif
png("/Users/nkmstar/Desktop/DILLSU_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(dillsu.1$gx, dillsu.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=dillsu.1$dbh/100, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/DILLSU_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(dillsu.2$gx, dillsu.2$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=dillsu.2$dbh/100, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/DILLSU_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(dillsu.3$gx, dillsu.3$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=dillsu.3$dbh/100, asp=1, main="2012")
dev.off()
# gif end

par(xaxs="i", yaxs="i")
plot(adindu.1$gx, adindu.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=adindu.1$dbh/100, asp=1, main= expression(italic("Adinandra dumosa")))
points(adindu.2$gx, adindu.2$gy, cex=adindu.2$dbh/150, asp=1, pch=16, col="red")
points(adindu.3$gx, adindu.3$gy, cex=adindu.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

# for making gif
png("/Users/nkmstar/Desktop/ADINDU_2004.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(adindu.1$gx, adindu.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=adindu.1$dbh/100, asp=1, main="2004")
dev.off()
png("/Users/nkmstar/Desktop/ADINDU_2008.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(adindu.2$gx, adindu.2$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=adindu.2$dbh/100, asp=1, main="2008")
dev.off()
png("/Users/nkmstar/Desktop/ADINDU_2012.png", width=299, height=550)
par(xaxs="i", yaxs="i")
plot(adindu.3$gx, adindu.3$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=adindu.3$dbh/100, asp=1, main="2012")
dev.off()
# gif end

par(xaxs="i", yaxs="i")
plot(streel.1$gx, streel.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=streel.1$dbh/100, asp=1, main= expression(italic("Streblus elongatus")))
points(streel.2$gx, streel.2$gy, cex=streel.2$dbh/150, asp=1, pch=16, col="red")
points(streel.3$gx, streel.3$gy, cex=streel.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(prunpo.1$gx, prunpo.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=prunpo.1$dbh/100, asp=1, main= expression(italic("Prunus polystachyus")))
points(prunpo.2$gx, prunpo.2$gy, cex=prunpo.2$dbh/150, asp=1, pch=16, col="red")
points(prunpo.3$gx, prunpo.3$gy, cex=prunpo.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(calowa.1$gx, calowa.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=calowa.1$dbh/100, asp=1, main= expression(italic("Calophyllum wallichianum")))
points(calowa.2$gx, calowa.2$gy, cex=calowa.2$dbh/150, asp=1, pch=16, col="red")
points(calowa.3$gx, calowa.3$gy, cex=calowa.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(elaems.1$gx, elaems.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=elaems.1$dbh/100, asp=1, main= expression(italic("Elaeocarpus mastersii")))
points(elaems.2$gx, elaems.2$gy, cex=elaems.2$dbh/150, asp=1, pch=16, col="red")
points(elaems.3$gx, elaems.3$gy, cex=elaems.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(fagrfr.1$gx, fagrfr.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=fagrfr.1$dbh/100, asp=1, main= expression(italic("Fagraea fragrans")))
points(fagrfr.2$gx, fagrfr.2$gy, cex=fagrfr.2$dbh/150, asp=1, pch=16, col="red")
points(fagrfr.3$gx, fagrfr.3$gy, cex=fagrfr.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(macaba.1$gx, macaba.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=macaba.1$dbh/100, asp=1, main= expression(italic("Macaranga bancana")))
points(macaba.2$gx, macaba.2$gy, cex=macaba.2$dbh/150, asp=1, pch=16, col="red")
points(macaba.3$gx, macaba.3$gy, cex=macaba.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

par(xaxs="i", yaxs="i")
plot(artoda.1$gx, artoda.1$gy, xlim=c(0,100), ylim=c(0,200), xlab="", ylab="", cex=artoda.1$dbh/100, asp=1, main= expression(italic("Artocarpus dadah")))
points(artoda.2$gx, artoda.2$gy, cex=artoda.2$dbh/150, asp=1, pch=16, col="red")
points(artoda.3$gx, artoda.3$gy, cex=artoda.3$dbh/200, asp=1, pch=16, col="blue")
legend(x=-10, y=70, legend=c("2004", "2008", "2012"), pch=c(1,16,16), bty="n", col=c("black", "red", "blue"), y.intersp = 0.5)

plot(x=1:3, y=c(nrow(campau.1), nrow(campau.2), nrow(campau.3)), type="l", xlab="", ylab="No. individuals")

# tried to insert inset graphs but failed
library(TeachingDemos)
subplot(plot(x=1:3, y=c(nrow(campau.1), nrow(campau.2), nrow(campau.3)), type="l", xlab="", ylab="No. individuals"), x=0, y=50, size=c(2,0.5))


# use Ryan's method to look at neutral vs environmental change



# look at effect of weather 








