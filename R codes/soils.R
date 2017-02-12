
soils <- read.delim("/Users/nkmstar/Dropbox/secplot_recovery/soil0-10.txt")

rownames(soils) <- soils$quadrat

soils <- soils[,-c(1,2)]

# remove 6 quadrats in primary forest
exclude <- c("K1", "K3", "O1", "P1", "L2", "N2")
soils <- soils[!(rownames(soils) %in% exclude),]

# clustering of all soil plots to see if there is clear delineation between primary and secondary forest
library(MASS)
soils.dist <- dist(soils)
soils.mds <- isoMDS(soils.dist)
plot(soils.mds$points, type = "n")
text(soils.mds$points, labels = as.character(rownames(soils))) # no clear clusters, N2 is an outlier
plot(hclust(soils.dist)) # quadrats all over the place, N2 is an outlier

# t-test to test for differences between pri and sec plot
soils$for.type <- c(rep("pri", 26), rep("sec", 20))
t.test(pH ~ for.type, data=soils)
t.test(Al ~ for.type, data=soils)
t.test(Ca ~ for.type, data=soils)
t.test(Fe ~ for.type, data=soils)
t.test(K ~ for.type, data=soils)
t.test(Mg ~ for.type, data=soils)
t.test(Mn ~ for.type, data=soils)
t.test(Na ~ for.type, data=soils)
t.test(CEC ~ for.type, data=soils)
t.test(BS ~ for.type, data=soils) # significant
t.test(pH.water ~ for.type, data=soils) # significant
t.test(pH.CaCl2 ~ for.type, data=soils)
t.test(resin.P ~ for.type, data=soils)
t.test(total.C ~ for.type, data=soils) # significant
t.test(total.N ~ for.type, data=soils)
t.test(C.N ~ for.type, data=soils) # significant
t.test(coarse.fragments ~ for.type, data=soils) # significant
t.test(bulk.density ~ for.type, data=soils)
t.test(fine.roots ~ for.type, data=soils) # significant


# mean of pri and sec plots
soils$forest <- c(rep("pri", 26), rep("sec", 20))
# mean and standard error
tapply(soils$pH, soils$forest, mean, na.rm=T); tapply(soils$pH, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$Al, soils$forest, mean, na.rm=T); tapply(soils$Al, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$Ca, soils$forest, mean, na.rm=T); tapply(soils$Ca, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$Fe, soils$forest, mean, na.rm=T); tapply(soils$Fe, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$K, soils$forest, mean, na.rm=T); tapply(soils$K, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$Mg, soils$forest, mean, na.rm=T); tapply(soils$Mg, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$Mn, soils$forest, mean, na.rm=T); tapply(soils$Mn, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$Na, soils$forest, mean, na.rm=T); tapply(soils$Na, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$CEC, soils$forest, mean, na.rm=T); tapply(soils$CEC, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$BS, soils$forest, mean, na.rm=T); tapply(soils$BS, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$pH.water, soils$forest, mean, na.rm=T); tapply(soils$pH.water, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$pH.CaCl2, soils$forest, mean, na.rm=T); tapply(soils$pH.CaCl2, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$resin.P, soils$forest, mean, na.rm=T); tapply(soils$resin.P, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$total.C, soils$forest, mean, na.rm=T); tapply(soils$total.C, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$total.N, soils$forest, mean, na.rm=T); tapply(soils$total.N, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$C.N, soils$forest, mean, na.rm=T); tapply(soils$C.N, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$coarse.fragments, soils$forest, mean, na.rm=T); tapply(soils$coarse.fragments, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$bulk.density, soils$forest, mean, na.rm=T); tapply(soils$bulk.density, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))
tapply(soils$fine.roots, soils$forest, mean, na.rm=T); tapply(soils$fine.roots, soils$forest, sd, na.rm=T)/sqrt(table(soils$forest))




