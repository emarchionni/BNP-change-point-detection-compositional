library(coda)
library(mcclust)
library(mcclust.ext)

source("C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/posterior inference/functions.R")
setwd("C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/posterior inference/MC_calibration")


PSM <- comp.psm(allocation_matrix)
minVI(PSM, allocation_matrix, method=("all"), include.greedy=TRUE)
minbinder(PSM, allocation_matrix, method=("all"), include.greedy=TRUE)