#==================================================================================================#
# Script created by Mark Christie, contact at markchristie1500@gmail.com
# Script created in version R 4.2.1
# This script:  Generates model output for rescue projects
# Usage notes: Set all parameters and then source this file (see markdown for details)
#==================================================================================================#
# set working directory, load packages

#setwd("C:/Users/fishf/Dropbox/manuscripts/rescue/salmon")
setwd("C:/mark/papers/pink_salmon_ibm/")
#library("hierfstat")
#library("inbreedR")

#==================================================================================================#
base.directory <- getwd() 
outdir         <- paste(base.directory,"/output/",sep="")  # directory to save model output  
source(paste(base.directory, "/source/FunctionSourcer.R", sep = '')) #loads packages, sources functions, and sets source directory
data.frame(unlist(parameters)) # show all parameters

n.replicates  <- 100           # number of replicates for each combination of parameters
replicates    <- Replicates(parameters, n.replicates)

for(i in 1:length(replicates)){
#i=1 
 #parameters[["k.adults.final"]] <- replicates[i, 1]  # final population size
  #parameters[["optima"]]         <- replicates[i, 2]  # phenotypic optimum for local adaptation
  parameters[["allele_freq_range"]] <- replicates[i]
  model <- RunModel(n.generations, i)

}

