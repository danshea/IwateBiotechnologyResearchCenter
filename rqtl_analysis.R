# Author: Dan Shea
# Date: 2019.01.07
# Description: R\qtl analysis of AKxPU F2 population of Lilium x formolongo
# crossed with Lilium longiflorum Marker set consists of GRAS-Di dominant and
# codominant markers. Dominant markers have been re-coded to codominant
# nomenclature. (See: Jupyter lab notebook GRAS-Di_data-combining.ipynb)

# Load the R\qtl library
library('qtl')
# Set the working directory
setwd('/Users/dj-shea/Desktop/Okazaki_lily')
# Load the F2 population marker data
mapthis <- read.cross(format='csv', dir='', file='rqtl.csv', genotypes=c("A", "B", "H"))
# Display a summary of the cross object
summary(mapthis)
# Examine missing genotype data
par(mfrow=c(1,1), las=1)
plotMissing(mapthis)

# The function ntyped() provides the numbers of genotyped markers for each
# individual (or the number of genotyped individuals for each marker). Let us
# plot these. (And note that there is a related function, nmissing(), which
# provides the number of missing genotypes for each individual or marker.)
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")

# From the plots, we see that some genotyping data is missing. However, no
# particular individual is missing large amounts of genotyping data (horizontal
# black lines). Instead, certain markers are largely absent of information
# within the population (vertical black lines).

# NOTE: No need to use subset to remove individuals from the population
# mapthis <- subset(mapthis, ind=(ntyped(mapthis) > cutoff_value))

# **** CHANGED TO ONLY INCLUDE MARKERS THAT GENOTYPED EVERY INDIVIDUAL ****
# We identify markers that are uninformative and remove them
nt.bymar <- ntyped(mapthis, 'mar')
todrop <- names(nt.bymar[nt.bymar < 125])
mapthis <- drop.markers(mapthis, todrop)
# NOTE: This removes "AMP0068746" and "AMP0071726" as less then 100 individuals
# were genotyped with this marker


# Identifying duplicate individuals
cg <- comparegeno(mapthis)
par(mfrow=c(1,1), las=1)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

# Usually, in an F2 population individuals will share genotypes at roughly 40% of
# markers. Here, we see the population as a whole skews towards closer to 65%.
# Since all individuals share less than 80% of genotypes with one another, we
# opt to not remove any individuals from the sample population.

# Identifying duplicate markers
# The function findDupMarkers identifies markers that have matching genotypes
# We may then remove these duplicate markers as they provide no new information.
print(dup <- findDupMarkers(mapthis, exact.only=FALSE, adjacent.only=TRUE))
# Current number of markers = 7,040
totmar(mapthis)
# Dropping duplicate adjacent markers
mapthis <- drop.markers(mapthis, unlist(dup))
# Results in 6,599 remaining markers
totmar(mapthis)

# Examining markers for segregation distortion
gt <- geno.table(mapthis)
# Markers showing segregation distortion
gt[gt$P.value < 0.05/totmar(mapthis),]
# markers that don't show segregation distortion
gt[gt$P.value > 0.05/totmar(mapthis),]
# Pretty much all of the markers fail to follow a 1:2:1 segregation
# This may be due to the fact that L. x formolongi (AK) was derived from
# continuously backcrossing an F1 plant (L. forosonanum x L. longiflorum) with
# L. longiflorum.
# Normally, we would drop markers that show severe distortion, but doing so in
# this case would remove all markers.

# Examining individuals in the population by genotype frequencies
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3) plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))

# Estimating recombination fractions
mapthis <- est.rf(mapthis)

# **** DO NOT RUN THIS PLOT IF YOU HAVE A LOT OF MARKERS! ****
# A plot of the LOD scores against the estimated recombination fractions
# for all marker pairs
# par(mfrow=c(1,1), las=1)
# rf <- pull.rf(mapthis)
# lod <- pull.rf(mapthis, what="lod")
# plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

# *****************************************************************************
# NOTE: This part should be run during the initial analysis when you are trying
# to determine appropriate settings for the maximum recombination fraction and
# log-odds. After that, you can comment it out.
# *****************************************************************************

# The output of formLinkageGroups() is a matrix with two columns: the initial
# linkage group or chromosome for each marker, and then the assigned linkage
# group, as inferred from the pairwise linkage information. The inferred linkage
# groups are numbered in decreasing order of size (so that linkage group 1 has
# the largest number of markers). Here we see that 11 linkage groups were
# inferred, with the last 2 groups having just one marker. One may play with
# max.rf and min.lod until the result is about as expected.

lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6)
table(lg[,2])
#    1    2    3
# 5832  411  356
lg <- formLinkageGroups(mapthis, max.rf=0.40, min.lod=6)
table(lg[,2])
#    1    2    3
# 5832  411  356
lg <- formLinkageGroups(mapthis, max.rf=0.40, min.lod=9)
table(lg[,2])
#    1    2    3
# 5832  411  356
lg <- formLinkageGroups(mapthis, max.rf=0.40, min.lod=3)
table(lg[,2])
#    1    2    3
# 5832  411  356
lg <- formLinkageGroups(mapthis, max.rf=0.25, min.lod=6)
table(lg[,2])
#    1    2    3
# 5832  411  356
lg <- formLinkageGroups(mapthis, max.rf=0.15, min.lod=6)
table(lg[,2])
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14
# 4193  411  373  356  267  258  221  206  196   44   31    7    6    3
#   15   16   17   18   19   20   21   22   23   24   25   26   27   28
#    2    2    2    1    1    1    1    1    1    1    1    1    1    1
#   29   30   31   32   33   34   35   36   37   38
#    1    1    1    1    1    1    1    1    1    1
lg <- formLinkageGroups(mapthis, max.rf=0.20, min.lod=6)
table(lg[,2])
#    1    2    3    4    5
# 5588  411  356  200   44
lg <- formLinkageGroups(mapthis, max.rf=0.15, min.lod=6)
table(lg[,2])
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14
# 4193  411  373  356  267  258  221  206  196   44   31    7    6    3
#   15   16   17   18   19   20   21   22   23   24   25   26   27   28
#    2    2    2    1    1    1    1    1    1    1    1    1    1    1
#   29   30   31   32   33   34   35   36   37   38
#    1    1    1    1    1    1    1    1    1    1

# We arrive at using max.rf = 0.15 and a min.lod = 6
# We can reorganize the markers into these inferred linkage groups. We do so
# with the same function, formLinkageGroups(), via the argument reorgMarkers.
mapthis <- formLinkageGroups(mapthis, max.rf=0.15, min.lod=6, reorgMarkers=TRUE)
# A plot of the pairwise recombination fractions and LOD scores may indicate how
# well this worked.
par(mfrow=c(1,1), las=1)
plotRF(mapthis, alternate.chrid=TRUE)

# Perform initial marker ordering on LGs (chromosomes) 2-13 (from smallest to largest)
# LG13 to LG2
for(i in 13:2) {
  mapthis <- orderMarkers(mapthis, chr=i, window=4, map.function = 'kosambi', verbose=2)
}

# Examine the orders determined for markers on each LG using pull.map()
pull.map(mapthis, chr=13)

# Now We can plot the LGs using plotMap
plotMap(mapthis, chr=c(2:13), show.marker.names=FALSE, xlab='Linkage Groups')
for(i in 2:13) {
  plotMap(mapthis, chr=i, show.marker.names=TRUE, xlab=sprintf("LG%d", i))
}

# Save the QTL experiment to file
write.cross(mapthis, format="csv")

# Final Notes
# This linkage analysis is not ideal, because the markers show severe 
# segregation distortion. This could be due to the nature of the parent,
# possibly the random amplification in GRAS-Di is biased by repetetive 
# regions of the Lily genome, etc. While biological data is usually not as
# clean as data sets used in tutorials to describe methodologies, the results
# obtained from this analysis are not very useful. However, the method used to
# obtain the results shown here is fine, and follows Professor Bromon's R\QTL
# documentation, which does utilize a nice clean data set for illustrative 
# purposes. - DJS
