Output <- function(pops, n, n.adult, output) {
  if(n == 1) {output <- NULL} # initialize output object as NULL
  juveniles <- pops[pops[, 4] == 1, ]
  adults    <- pops[pops[, 4] == 2, ]
  n.loci    <- parameters[["n.loci"]]
  n.neutral.loci    <- parameters[["n.neutral.loci"]]
   
  n.juveniles <- length(juveniles[, 1])
  nadult.1    <- n.adult
  nadult.2    <- nrow(adults) # adults right now in the population (after mortality)
  
  dims        <- ncol(pops)-(n.neutral.loci*2)-n.loci # seprate out adaptive loci first
#  phenotypes  <- pops[, (dims+1):(dims+n.loci)]
#  phenotypes  <- mean(rowSums(phenotypes))  + pops[1, 8] # takes first d.dev value
  gtypes      <- pops[, 10:11]
  
  alleles     <- c(gtypes[, 1], gtypes[, 2]) # ALLELE FREQ IN ALL INDIVIDS INCLUDING JUVS
  allelea     <- length(which(alleles == 1))/length(alleles)
  alleleb     <- length(which(alleles == 2))/length(alleles)
  

#  ne1 <- 1/(sum(1/output[, 2])/nrow(output))  # harmonic mean of adults
#  ne2 <- 1/(sum(1/output[, 3])/nrow(output))  # harmonic mean of adults.2  
#  ne3 <- 1/(sum(1/(output[, 2] + output[, 3]))/nrow(output)) # harmonic mean of adults + juveniles
  
  ne1 <- 1/(sum(1/tail(output[, 2], 5))/5)  # harmonic mean of juveniles
  ne2 <- 1/(sum(1/tail(output[, 3], 5))/5)  # harmonic mean of adults.2  
  ne3 <- 1/(sum(1/(tail(output[, 2],5) + tail(output[, 3],5)))/5) # harmonic mean of adults + juveniles
  
  ne2 <- tail(output[, 3], 5)
  ne2 <- ne2[-which(ne2 == 0)]
  ne2 <- 1/((sum(1/ne2))/length(ne2))
#  adult.ids <- adults[, 1]
#  families  <- unlist(strsplit(adult.ids, "_"))
#  f.years   <- families[seq(1, length(families), 4)] # query here for mean parent age
#  f.fam     <- families[seq(2, length(families), 4)] # query here for mean parent age
#  family    <- paste(f.years, f.fam, sep="_") 
  #hist(table(family))  family sizes
#  f.sizes <- data.frame(table(family))
#  k  <- mean(f.sizes[, 2])
#  vk <- var(f.sizes[, 2])
#  N  <- nrow(f.sizes)
#  ne4 <- ((N*k)-1)/(k-1+(vk/k))
  
  ages    <- as.character(table(pops[, 3]))        # counts of the below age categories
  ages    <- paste(ages,sep="", collapse="/")
  agevals <- names(table(pops[, 3]))               # actual age categories
  agevals <- paste(agevals,sep="", collapse="/")
  
#  k.adults <- parameters[["k.adults"]]
  
  if(n==1) {ne1=ne2=ne3=0} # nes for first generation for formatting
  #output.new <- data.frame(cbind(n, n.juveniles, nadult.1, nadult.2, phenotypes, ho, fis, ne1, ne2, ne3, ne4, inbreed, het2, F_inbreeding, agevals, ages))
  output.new <- data.frame(n, n.juveniles, nadult.1, nadult.2, allelea, alleleb, ne1, ne2, ne3, agevals, ages)
  output     <- rbind(output, output.new)
  return(output)
}
