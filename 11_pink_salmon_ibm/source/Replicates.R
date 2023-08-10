Replicates <- function(parameters, n.replicates) {
  # takes all *.set variables and creates the precise number of replicates for each variable
 
 # parameters for donor populations
 #BC alleles
# allele_freq_range<- c(0.23214286, 0.03448276, 0.96551724, 0.79310345, 0.80357143, 0.91379310, 0.40740741, 0.81034483, 0.93103448, 0.77586207, 0.94827586,0.82142857, 0.84482759, 0.98275862, 0.82758621, 0.00000000)
 
 # alleles for GL after runing BC sims
 allele_freq_range<- c(0.33333333, 0.08333333, 0.91666667, 0.68333333, 0.70000000, 0.85000000, 0.51666667, 0.71666667, 0.86666667, 0.66666667, 0.90000000, 0.73333333, 0.75000000, 0.95000000, 0.00000000)
 
 #replicates <- expand.grid(k.adults.final, optima)
 replicates <- as.data.frame(allele_freq_range)
 replicates <- replicates[rep(seq_len(nrow(replicates)), n.replicates), ]
 
 return(replicates)
}  
