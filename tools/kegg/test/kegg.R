#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------------------------------------------------
#- R script for pathway analysis using KEGG
#-------------------------------------------------------------------------------------------------------------------------

#----------------------------------
#- Load / install required packages
#----------------------------------

library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)