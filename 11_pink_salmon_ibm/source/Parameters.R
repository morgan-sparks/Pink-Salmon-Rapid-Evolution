# Note: this is the default file; all variables are manipulated in Replicates.R
# If a parameter has a value of NA; it is varied directly in Replicates.R

#Input parameters
parameters <- list()          # all parameters must be added to this list (see code below section break : ======)

## Demographic and life history parameters =======================#
pops                      <- 1         # total n. pops
k.adults                  <- 500       # 500
k.juvenile.multiplier     <- 7       # how many times more juveniles should there be than adults
k.juveniles 	            <- k.adults*k.juvenile.multiplier   # set as proportion of adults (easier to adjust population sizes to ne)
k.adults.final            <- NA
juvenile.survival.var 	  <- 0.1       # default  = 100 , now set as proportion of k.juveniles 
adult.survival.var        <- 0.09       # used to be 80, now set as proportion of k.adults

# life history parameters 
maximum.juvenile.age   <- 1     # maximim age for a juvenile lamprey; used in creating juveniles (Juveniles.R) 
maximum.adult.age      <- 2     # maximum age for adult     
semelparity            <- 0     # 1 = semelparous; 2 (or anything but 1) = iteroparous ; if on, adults cannot be donated via rescue
breed.next.year        <- 0.30  # when semelparity is on, what percentage of new "adults" should remain juveniles for 1 extra year to create 4 year old coho
k.juv.year             <- cbind(0.6, 0.3) # must sum to one and must have as many elements as all juvenile year classes; ordered youngest to oldest



#Timing parameters (modified from donor populations; no selection for genetic rescue model, change for ER)
selection.start.year <- 5000    # year that selection starts
ne.change.start      <- 1   # selection.end.year + 1
ne.change.stop       <- ne.change.start + 199
n.years              <- 254   # total number of years in the model
selection.end.year   <- 5000  # 150 is default when 50 is selection.start.year

  
# Growth parameters
##y   <- L * (1-e^(-K * (t-tzero)))
maximum.juvenile.size <- 150 # equal to L in above equation
maximum.adult.size    <- 800 # equal to L in above equation
transform.threshold   <- 80 #20 # size at which all juveniles above will transform into juveniles
adult.threshold       <- 450 # size at which all juveniles will transform into adults
shape.parameter       <- 0.0011 # equal to K in above equation
shape.var             <- 0.0002 # how much within-population, between-individual variation do we want. analagous to k.sd in growth.R
e                     <- exp(1) # equal to e in above equation
tzero                 <- -15   # size at age 0, all sizes in mm, shifts lines left or right along x-axis (could add variation in egg size/provisioning by varying this parameter)


# Density dependent egg production
egg.multiplier <- 1   # 1 and 10 seem like good starting places
egg.addition   <- 5   # 5 is default; 150 make it flat
egg.max        <- 150 # maximum number of eggs per pair ()

#plot
#n.egg           <- 1000/((1:1000)/2)
#n.eggs          <- (n.egg * egg.multiplier) + egg.addition
#plot(1:1000, n.eggs, ylim = c(0, 150))

# adaptive genetic parameters
n.loci        	       <- 1   # n adaptive loci
n.alleles              <- 2    # n. alleles per adaptive locus
dominance.deviation    <- 0    # starting position for dominance deviation (to be updated throughout model run)  
dominance.effect       <- 0.05 # 0.25 default, absolute value of dominance effects for loci (all equal for now, but can update)   
environmental.efffect  <- 1    # strength of environmental contribution to phenotype (std of normal distribution)

# neutral genetic parameters
n.neutral.loci   <- 1
n.neutral.alleles<- 2
allele_freq_range <- 0.9  # set in Replicates (starting allele frequencies to test)

#Selection parameters
#selection.strength   <- 0.4  # proportion that die   (higher = more mortality, range = 0 to 1); currently not used
optima                <- 0   # phenotpyic optimum in new environment; here modeled as increase (Celsius) in temp tolerance. so 2 would mean 2 degrees higher than current temp conditions. dictactes selection strength in each generation. linear.

# rescue parameters
rescue.years  <- 0#c(80:90)   # set to 0 if want no rescue, years to perform rescue
rescue.number <- 20   # number of individuals to move  
rescue.age    <- 3    # stage of individuals to move 

# Inbreeding parameters
inbreeding.fitness.cost <- 0 # inbreeding fitness cost (implemented in Adult Mortality) is turned on (1) or off (not 1, e.g., 0)


# Add all parameters to list===============================================================#
parameters[["n.years"]] <- n.years
parameters[["selection.start.year"]] <- selection.start.year
parameters[["selection.end.year"]]   <- selection.end.year  
parameters[["ne.change.start"]]      <- ne.change.start 
parameters[["ne.change.stop"]]       <- ne.change.stop  

parameters[["pops"]] <- pops

parameters[["k.juveniles"]]           <- k.juveniles 
parameters[["k.juvenile.multiplier"]] <- k.juvenile.multiplier
parameters[["maximum.juvenile.age"]]  <- maximum.juvenile.age 
parameters[["maximum.adult.age"]]     <- maximum.adult.age
parameters[["semelparity"]]           <- semelparity  
parameters[["breed.next.year"]]       <- breed.next.year


parameters[["juvenile.survival.var"]] <- juvenile.survival.var 
parameters[["adult.survival.var"]]    <- adult.survival.var 
parameters[["egg.max"]]               <- egg.max
parameters[["egg.multiplier"]]        <- egg.multiplier
parameters[["egg.addition"]]          <- egg.addition

 
parameters[["maximum.adult.size"]]  <- maximum.adult.size  
parameters[["transform.threshold"]] <- transform.threshold
parameters[["adult.threshold"]]     <- adult.threshold
parameters[["shape.parameter"]]  <- shape.parameter   
parameters[["e"]] <-e                   
parameters[["tzero"]] <- tzero          
parameters[["shape.var"]]<-shape.var

parameters[["n.loci"]] <- n.loci 
parameters[["n.alleles"]] <- n.alleles 
parameters[["dominance.deviation"]] <- dominance.deviation
parameters[["dominance.effect"]] <- dominance.effect
parameters[["n.neutral.loci"]] <- n.neutral.loci
parameters[["n.neutral.alleles"]] <- n.neutral.alleles

parameters[["k.adults"]] <- k.adults

parameters[["egg.multiplier"]] <-  egg.multiplier 
parameters[["egg.addition"]]   <-  egg.addition
parameters[["egg.max"]]        <-  egg.max

#parameters[["selection.strength"]]    <- selection.strength   
parameters[["optima"]]                <- optima   
parameters[["environmental.efffect"]] <- environmental.efffect
parameters[["k.adults.final"]]        <- k.adults.final  

parameters[["inbreeding.fitness.cost"]] <- inbreeding.fitness.cost

parameters[["rescue.years"]]  <- rescue.years
parameters[["rescue.number"]] <- rescue.number
parameters[["rescue.age"]] <- rescue.age

parameters[["k.juv.year"]] <- k.juv.year

parameters[["allele_freq_range"]] <- allele_freq_range 

