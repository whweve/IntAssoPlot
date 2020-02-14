## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 7, fig.align = "center")
library(IntAssoPlot)

## -----------------------------------------------------------------------------
head(association)
#the attribute of each column could be viewed as:
str(association$Marker)
str(association$Locus)
str(association$Site)
str(association$p)

## -----------------------------------------------------------------------------
head(gtf)
#the attribute of each column could be viewed as:
str(gtf$V1)
str(gtf$V2)
str(gtf$V3)
str(gtf$V4)
str(gtf$V5)
str(gtf$V6)
str(gtf$V7)
str(gtf$V8)
str(gtf$V9)

## -----------------------------------------------------------------------------
#only 20 column of the genotype markers are shown.
head(zmvpp1_hapmap[,1:20])

## -----------------------------------------------------------------------------
IntRegionalPlot(chr=9,left=94178074-200000,right=94178074+200000,gtf=gtf,association=association,hapmap=hapmap_am368,hapmap_ld=hapmap_am368,threshold=5,leadsnp_size=2)

## -----------------------------------------------------------------------------
IntRegionalPlot(chr=9,left=94178074-200000,right=94178074+200000,gtf=gtf,association=association,hapmap=hapmap_am368,hapmap_ld=hapmap_am368,threshold=5,leadsnp_size=2,colour02 = "gray1",colour04 = "gray21",colour06 = "gray41",colour08 = "gray61",colour10 = "gray81",)

## -----------------------------------------------------------------------------
#get five colors ranging from white to red
pal <- colorRampPalette(c("white", "red"))
IntRegionalPlot(chr=9,left=94178074-200000,right=94178074+200000,gtf=gtf,association=association,hapmap=hapmap_am368,hapmap_ld=hapmap_am368,threshold=5,leadsnp_size=2,colour02 = pal(5)[1],colour04 = pal(5)[2],colour06 = pal(5)[3],colour08 = pal(5)[4],colour10 = pal(5)[5])

## -----------------------------------------------------------------------------
#get five colors ranging from white to red
pal <- colorRampPalette(c("white", "red"))
IntRegionalPlot(chr=9,left=94178074-200000,right=94178074+200000,gtf=gtf,association=association,hapmap=hapmap_am368,hapmap_ld=hapmap_am368,threshold=5,leadsnp_size=2,colour02 = pal(5)[1],colour04 = pal(5)[2],colour06 = pal(5)[3],colour08 = pal(5)[4],colour10 = pal(5)[5],label_gene_name = TRUE)

## -----------------------------------------------------------------------------
IntRegionalPlot(chr=9,left=94178074-100000,right=94178074+100000,gtf=gtf,association=association,hapmap=hapmap_am368,hapmap_ld=hapmap2,threshold=5,leadsnp_size=2)

## -----------------------------------------------------------------------------
IntRegionalPlot(chr=9,left=94178074-2000,right=94178074+5000,gtf=gtf,association=association,hapmap=hapmap_am368,hapmap_ld=hapmap_am368,threshold=5,leadsnp_size=2)

## -----------------------------------------------------------------------------
IntGenicPlot('GRMZM2G170927_T01',gtf,association=zmvpp1_association,hapmap=zmvpp1_hapmap,hapmap_ld = zmvpp1_hapmap,threshold=8,leadsnpLD = FALSE)

## -----------------------------------------------------------------------------
IntGenicPlot('GRMZM2G170927_T01',gtf,association=zmvpp1_association,hapmap=zmvpp1_hapmap,hapmap_ld = zmvpp1_hapmap,threshold=8,up=500,down=600,leadsnpLD = FALSE)

## -----------------------------------------------------------------------------
IntGenicPlot('GRMZM2G170927_T01',gtf,association=zmvpp1_association,hapmap=zmvpp1_hapmap,hapmap_ld = zmvpp1_hapmap,threshold=8,up=500,down=600,leadsnpLD = FALSE,marker2highlight=marker2highlight)

## -----------------------------------------------------------------------------
IntGenicPlot('GRMZM2G170927_T01',gtf,association=zmvpp1_association,hapmap=zmvpp1_hapmap,hapmap_ld = zmvpp1_hapmap,threshold=8,up=500,down=600,leadsnpLD = FALSE,marker2highlight=marker2highlight,link2gene=marker2link,link2LD=marker2link)

## -----------------------------------------------------------------------------
IntGenicPlot('GRMZM2G170927_T01',gtf,association=zmvpp1_association,hapmap=zmvpp1_hapmap,hapmap_ld = zmvpp1_hapmap,threshold=8,up=500,down=600,leadsnpLD = FALSE,marker2highlight=marker2highlight,link2gene=marker2link,link2LD=marker2link,marker2label=marker2link,marker2label_angle=60,marker2label_size=2)

