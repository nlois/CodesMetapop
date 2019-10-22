rm(list=ls())
graphics.off()

#
#### AMOVA through poppr package!
# Based on analysis for Rockhopper ddRAD data from SWAO colonies
#
# Need to install packages 1 by 1 (at least in R 3.5)

library(ape)
library(ade4)
library(adegenet)
library(hierfstat)
library(poppr)

# Before running code!!
# populations makes a *.tsv Structure file, so we need to change it to *.str
# ... mv batch_1.structure_allSNPs.tsv batch_1.structure_allSNPs.str
# also replace 0 for -9 (missing data in structure format!)
# AND!!!! also need to set the structure of pop in 2nd column 
# (may use _ for nesting 1_1, 1_2, 2_1, 3_1, etc... or 1_1_1, etc..)
# Then with splitStrata() you split using pop structure.



# NORTH

setwd("/Users/user/Documents/Rockhopper/Sep2018/NORTH")

# When importing data drom str format you need further arguments to fill in
# the console. 
# You can also write them down in here and command + enter your way.
# individuals, loci, sample name column, pop column, then enter in the
# console, locus name row, finally "no" for same line alleles.

NORTH <- import2genind(file = "./batch_1.structure_allSNPs.str") 
26
4572
1
2

1
n


#Fis por pop
FIS <- boot.ppfis(NORTH,nboot = 10000)
FIS

# Fis
NORTHdf <- data.frame(NORTH@pop,NORTH)
NORTHstats <- basic.stats(NORTH)
names(NORTHstats)
NORTHstats$Fis

# Check if you don't have 0 as an allele (missing data)
NORTH@all.names
NORTH

# Set pop as stratum for AMOVA
strata(NORTH) <- data.frame(NORTH@pop)
NORTH <- as.genclone(NORTH)
NORTH

# split strata according to populations
splitStrata(NORTH) <- ~Pop/Subpop
NORTH

# Check pop definition...
table(strata(NORTH, ~Pop/Subpop))

# Run AMOVA (could add more parameters, corrections, etc!!)
NORTHAMOVA <- poppr.amova(NORTH, ~Pop/Subpop)
NORTHAMOVA

# Significance test through 1000 permutations
set.seed(1999)
NORTHsig   <- randtest(NORTHAMOVA, nrepet = 999)
NORTHsig
plot(NORTHsig)

# RESULTS
# Variations between NORTH.pop   greater  p = 0.001
#
#
#
#


### SOUTH

setwd("/Users/user/Documents/Rockhopper/Sep2018/SOUTH")

SOUTH <- import2genind(file = "./batch_1.structure_allSNPs.str") 
81
4918
1
2

1
n



SOUTH@all.names
SOUTH

strata(SOUTH) <- data.frame(SOUTH@pop)
SOUTH <- as.genclone(SOUTH)
SOUTH

# split strata according to populations
splitStrata(SOUTH) <- ~Pop/Subpop
SOUTH

# Check pop definition...
table(strata(SOUTH, ~Pop/Subpop))

# Run AMOVA (could add more parameters, corrections, etc!!)
SOUTHAMOVA <- poppr.amova(SOUTH, ~Pop/Subpop)
SOUTHAMOVA

set.seed(1999)
SOUTHsig   <- randtest(SOUTHAMOVA, nrepet = 999)
plot(SOUTHsig)
SOUTHsig
#Variations between SOUTH.pop   p = 0.08



## ALL

setwd("/Users/user/Documents/Rockhopper/Sep2018/popGen_dataSET")

ALL <- import2genind(file = "./batch_1.structure_allSNPs.str") 
107
4975
1
2

1
n



head (ALL@all.names)
ALL@pop
ALL@strata


strata(ALL) <- data.frame(ALL@pop)
ALL <- as.genclone(ALL)
ALL

splitStrata(ALL) <- ~Pop/Subpop
ALL

table(strata(ALL, ~Pop/Subpop/Col))

ALLAMOVA <- poppr.amova(ALL, ~Pop/Subpop/Col)
ALLAMOVA

set.seed(1999)
ALLsig   <- randtest(ALLAMOVA, nrepet = 999)
plot(ALLsig)
ALLsig
#Variations between ALL.pop   p = 0.03
#Variations between X within ALL   p = 0.3


## ALL for basic gen stats
setwd("/Users/user/Documents/Rockhopper/Sep2018/popGen_dataSET")
#### Poppr estimations
ALL <- import2genind(file = "./AllSNPs_struct_Pop_Subpop.str") 
107
4975
1
2

1
n

ALL@pop

#Fst (fstat() and pairwise.fst() from hierfstat and adegenet)
fstats <- fstat(x = ALL, fstonly = FALSE)
pairwise.fst(ALL,pop = NULL)

# plot Fis with range from Bootstrap
FisCI <- read.csv("./FisCI95Boot.csv")
FisCI
ggplot() +
  geom_pointrange(data=FisCI,mapping = aes( x=FisCI$X,
                                            y=FisCI$Fis,ymin=FisCI$CI.low,
                                            ymax=FisCI$CI.high))
## options(title("Fis"),ylab("Fis"),xlab(Colonie))

                  

head (ALL@all.names)
ALL
ALL@pop
ALL@strata

strata(ALL) <- data.frame(ALL@pop)
ALL <- as.genclone(ALL)
ALL

splitStrata(ALL) <- ~Pop/Subpop
ALL

# Diversity stats!
ALLtab <- mlg.table(ALL, strata = ~Pop/Subpop)
diversity_stats(ALLtab)
diversity_ci(ALLtab, n = 10000, raw = FALSE,plot = TRUE)
?diversity_ci

#hierfstat basic stats
div <- summary(NORTH)
div
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus",
     cex=0.25)
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of
     observed heterozygosity per locus", cex=0.1)
bartlett.test(list(div$Hexp, div$Hobs))

ALLdf <- data.frame(ALL@pop,ALL)
ALLstats <- basic.stats(ALL)
ALLstats$Ho
names(ALLstats)
ALLstats$Fis
write.csv(ALLstats$Fis,file = "./ALLFis")

# Fis per loc, per pop, seeking mean
mean(ALLstats$Fis[,1], na.rm = TRUE)
mean(ALLstats$Fis[,2], na.rm = TRUE)
mean(ALLstats$Fis[,3], na.rm = TRUE)
mean(ALLstats$Fis[,4], na.rm = TRUE)
mean(ALLstats$Fis[,5], na.rm = TRUE)
ALL@pop
mean_Fis <- data.frame(c("IM/F","IDLEBF","IP","Chi","IDLESJ"), c(mean(ALLstats$Fis[,1], na.rm = TRUE), mean(ALLstats$Fis[,2], na.rm = TRUE),mean(ALLstats$Fis[,3], na.rm = TRUE),mean(ALLstats$Fis[,4], na.rm = TRUE),mean(ALLstats$Fis[,5], na.rm = TRUE)))
mean_Fis

# Hobs per loc, per pop, seeking mean
mean(ALLstats$Ho[,1], na.rm = TRUE)
mean(ALLstats$Ho[,2], na.rm = TRUE)
mean(ALLstats$Ho[,3], na.rm = TRUE)
mean(ALLstats$Ho[,4], na.rm = TRUE)
mean(ALLstats$Ho[,5], na.rm = TRUE)
mean_Ho <- data.frame(c("IM/F","IDLEBF","IP","Chi","IDLESJ"), c(mean(ALLstats$Ho[,1], na.rm = TRUE), mean(ALLstats$Ho[,2], na.rm = TRUE),mean(ALLstats$Ho[,3], na.rm = TRUE),mean(ALLstats$Ho[,4], na.rm = TRUE),mean(ALLstats$Ho[,5], na.rm = TRUE)))

##
#### Poppr estimations for population classified as N vs South
ALL <- import2genind(file = "./AllSNPs_struct_Pop.str") 
107
4975
1
2

1
n

ALL@pop

#Fst (fstat() and pairwise.fst() from hierfstat and adegenet)
fstats <- fstat(x = ALL, fstonly = FALSE)
pairwise.fst(ALL,pop = NULL)

# Basic stats for 
ALLstats <- basic.stats(ALL)
ALLstats$Ho
names(ALLstats)
ALLstats$Fis
write.csv(ALLstats$Fis,file = "./ALLFis")

# Fis per loc, per pop, seeking mean
mean(ALLstats$Fis[,1], na.rm = TRUE)
mean(ALLstats$Fis[,2], na.rm = TRUE)
ALL@pop
mean_Fis <- data.frame(c("NORTH","SOUTH"), c(mean(ALLstats$Fis[,1], na.rm = TRUE), mean(ALLstats$Fis[,2])))
mean_Fis

# Hobs per loc, per pop, seeking mean
mean(ALLstats$Ho[,1], na.rm = TRUE)
mean(ALLstats$Ho[,2], na.rm = TRUE)
mean_Ho <- data.frame(c("NORTH","SOUTH"), c(mean(ALLstats$Ho[,1], na.rm = TRUE), mean(ALLstats$Ho[,2], na.rm = TRUE),mean(ALLstats$Ho[,3], na.rm = TRUE),mean(ALLstats$Ho[,4], na.rm = TRUE),mean(ALLstats$Ho[,5], na.rm = TRUE)))

summary <- poppr(ALL)


## More measures
ALL$pop
#Malv,BF,PD,Chile,SJ
FIS_boot <- boot.ppfis(ALL,nboot = 1000)
FIS_boot
write.csv(FIS$fis.ci, file="./FisCIBoot.csv")
FST_Boot <- boot.ppfst(ALL,nboot = 1000,diploid = TRUE,quant = c(0.05,0.95))
plot(x = FIS$hl, c("M","",3,4,5))

summary2 <- poppr(ALL, strata = Pop/Subpop,quiet=FALSE,minsamp = 4)


# HWE
library("pegas")
hwe.full <- hw.test(ALL, B = 100) # performs 1000 permuatations
hwe.fulldf <- data.frame(hwe.full)
plot(x = hwe.fulldf$chi.2, y = hwe.fulldf$Pr.exact)
hwe.full
# Tendriamos que sacar estos SNPs que no estan en HWE ??

# Pop struct
micdist    <- provesti.dist(ALL)
theTree <- micdist %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
plot(theTree)
add.scale.bar(length = 0.01) # add a scale bar showing 5% difference.
set.seed(999)
aboot(ALL, dist = provesti.dist, sample = 200, tree = "nj",
      cutoff = 50, quiet = TRUE)




##
##
## ALL_DEAD

setwd("/Users/user/Documents/Rockhopper/Sep2018/popGen_dataSET_andDEAD")

ALL_DEAD <- import2genind(file = "./batch_1.structure_allSNPs.str") 

135
4705
1
2
3

1
n



head (ALL_DEAD@all.names)
ALL_DEAD

strata(ALL_DEAD) <- data.frame(ALL_DEAD@pop,ALL_DEAD@other)
ALL_DEAD <- as.genclone(ALL_DEAD)
ALL_DEAD

table(strata(ALL_DEAD, ~ALL_DEAD.pop/X))

ALL_DEADAMOVA <- poppr.amova(ALL_DEAD, ~ALL_DEAD.pop/X)
ALL_DEADAMOVA

set.seed(1999)
ALL_DEADsig   <- randtest(ALL_DEADAMOVA, nrepet = 999)
ALL_DEADsig
plot(ALL_DEADsig)

#Variations between ALL_DEAD.pop   p = 0.02
#Variations between X within ALL_DEAD   p = 0.17






## NORTH_andDEAD

setwd("/Users/user/Documents/Rockhopper/Sep2018/NORTH_andDEAD")

NORTH_andDEAD <- import2genind(file = "./batch_1.structure_allSNPs.str") 

54
4140
1
2
3

1
n



head (NORTH_andDEAD@all.names)
NORTH_andDEAD

strata(NORTH_andDEAD) <- data.frame(NORTH_andDEAD@pop,NORTH_andDEAD@other)
NORTH_andDEAD <- as.genclone(NORTH_andDEAD)
NORTH_andDEAD

table(strata(NORTH_andDEAD, ~NORTH_andDEAD.pop/X))

NORTH_andDEADAMOVA <- poppr.amova(NORTH_andDEAD, ~NORTH_andDEAD.pop/X)
NORTH_andDEADAMOVA

set.seed(1999)
NORTH_andDEADsig   <- randtest(NORTH_andDEADAMOVA, nrepet = 999)
NORTH_andDEADsig
plot(NORTH_andDEADsig)

#Variations between NORTH_andDEAD.pop   p = 0.2
#Variations between X within NORTH_andDEAD   p = 0.17
#
#
#



## NORTH_andDEAD_PD

setwd("/Users/user/Documents/Rockhopper/Sep2018/NORTH_andDEAD_PD")

NORTH_andDEAD_PD <- import2genind(file = "./batch_1.structure_allSNPs.str") 

54
4140
1
2
3

1
n



head (NORTH_andDEAD_PD@all.names)
NORTH_andDEAD_PD

strata(NORTH_andDEAD_PD) <- data.frame(NORTH_andDEAD_PD@pop,
                                       NORTH_andDEAD_PD@other)
NORTH_andDEAD_PD <- as.genclone(NORTH_andDEAD_PD)
NORTH_andDEAD_PD

table(strata(NORTH_andDEAD_PD, ~NORTH_andDEAD_PD.pop/X))

NORTH_andDEAD_PDAMOVA <- poppr.amova(NORTH_andDEAD_PD, ~NORTH_andDEAD_PD.pop/X)
NORTH_andDEAD_PDAMOVA

set.seed(1999)
NORTH_andDEAD_PDsig   <- randtest(NORTH_andDEAD_PDAMOVA, nrepet = 999)
NORTH_andDEAD_PDsig
plot(NORTH_andDEAD_PDsig)

#Variations between NORTH_andDEAD_PD.pop   p = 0.2
#Variations between X within NORTH_andDEAD_PD   p = 0.17


# OUTPUT !
setwd("/Users/user/Documents/Rockhopper/Sep2018")
pdf("AMoVAsig.pdf",width=10,height=10)
par(mfrow=c(1,1))
plot(SOUTHsig)
plot(NORTHsig)
plot(ALLsig)
plot(ALL_DEADsig)
plot(NORTH_andDEADsig)
plot(NORTH_andDEAD_PDsig)

dev.off()


SOUTHAMOVA
NORTHAMOVA
ALLAMOVA
ALL_DEADAMOVA
NORTH_andDEADAMOVA
NORTH_andDEAD_PDAMOVA






###
# Fst histogram 


# ggplot web tutorials !!!
# http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
# Normalizing y axis
# https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
# set scale limits
# https://ggplot2.tidyverse.org/reference/lims.html

library(ggplot2)

setwd("/Users/user/Documents/Rockhopper/Sep2018/AMOVA")
# Load Fst distrib per pop
FstPerLoc<- read.delim("./FstPerLocus.txt")
# Change histogram plot line colors by groups
ggplot(data = FstPerLoc, aes(FstPerLoc$Fst.per.locus)) +
  geom_histogram(fill="black",binwidth=0.01,alpha=1) + 
ylab("Loci")+
  xlab("Fst")





###### Just trying things out... Don't run this part !!

NORTHdf <- genind2df(NORTH,oneColPerall = TRUE)
df2genind(NORTHdf,NA.char = 0,ploidy = 2,type = "codom",sep = "\t")

class?genind

poppr(dat = NORTH)

genotype_curve(gen = NORTH, sample = 1000, quiet = TRUE)

info_table(NORTH, type = "missing", plot = TRUE)

pop <- c("1","1","1","1","1","1","1","1","1","1","1",
         "1","2","2","2","2","2","2","2","2","2","2","2","2","2","2")

addStrata(NORTH, name = "POP") <- pop
NORTH@strata

poppr.amova(NORTH, within = FALSE,dist = NULL,freq = TRUE,
            filter = TRUE, hier=~NORTH.pop)
