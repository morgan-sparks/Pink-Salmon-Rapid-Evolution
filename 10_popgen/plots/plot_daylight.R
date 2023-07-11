# Header ------------------------------------------------------------------

setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/GL_Pink_Salmon/")



library(geosphere); library(tidyverse); library(lubridate); library(ggdist); library(patchwork)



# Get daylight data -------------------------------------------------------



lakelse.r <- c(54.379998, -128.549272)

steel.r <- c(48.777236, -86.886525)

lakelse.r.dl <- daylength(lakelse.r[1], doy = c(1:365))
lakelse.r.dl <- data.frame(cbind(day = c(1:365), day.length = lakelse.r.dl))

steel.r.dl <- daylength(steel.r[1], doy = c(1:365))
steel.r.dl <- data.frame(cbind(day = c(1:365), day.length = steel.r.dl))




ggplot() + 
  geom_line(data = steel.r.dl, aes(x = day, y = day.length, color = "Great Lakes")) +
  geom_line(data = lakelse.r.dl, aes(x = day, y = day.length, color = "British Columbia")) +
  scale_color_manual(values = c( "Great Lakes" = "darkorange2", 
                                 "British Columbia" = "blue4")) +
  theme_classic() 


#### plot as difference 
dl.comp <- steel.r.dl %>% 
  dplyr::select(day) %>% 
  mutate(daylight.diff = steel.r.dl$day.length- lakelse.r.dl$day.length,
         day.date = as_date(day))



ggplot(dl.comp)+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax =1.2, fill = "grey10", alpha = 0.25 ) +
  annotate("rect", xmin = 253, xmax = 283, ymin = -Inf, ymax =1, fill = "forestgreen", alpha = 0.5 ) +
  annotate("rect", xmin = 122, xmax = 152, ymin = -Inf, ymax =1, fill = "blue4", alpha = 0.5 ) +
  annotate("text", x = 182.5, y =1.1, label = "Lake Growth", size = 2.5 ) +
  annotate("text", x = 80, y = .8, label = "Emergence and\nOutmigration", size = 2.5) +
  annotate("text", x = 220, y = .8, label = "Spawning", size = 2.5) +
  geom_line(aes(x = day, y= daylight.diff), size = .8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # scale_x_date(date_labels = "%b") +
  labs(x = " Day of Year", y = "Hours of light difference") +
  theme_classic()


#ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/Fig 5.pdf", fig.5, height = 2.25, width = 4.5, units = "in", dpi = 600)



# make two year fig -------------------------------------------------------
year.1 <- dl.comp

year.2 <- dl.comp %>% 
  mutate(day.date = day.date+ years(1))

two.years <- rbind(year.1, year.2)

emergence <- two.years[122:152,]
lake <- two.years[152:618,]
spawn <- two.years[618:648,]
  
daylight<- ggplot()+
  # geom_rect(data = two.years[122:152,], aes(xmin = min(day.date), xmax = max(day.date), ymin = -Inf, ymax = 1), fill = "dodgerblue", alpha = 0.25) +
  # geom_rect(data = two.years[152:618,], aes(xmin = min(day.date), xmax = max(day.date), ymin = -Inf, ymax = 1), fill = "grey80", alpha = 0.25 ) +
  # geom_rect(data = two.years[618:648,], aes(xmin = min(day.date), xmax = max(day.date), ymin = -Inf, ymax = 1), fill = "green", alpha = 0.25 ) +
  annotate("rect", xmin = min(emergence$day.date), xmax = max(emergence$day.date), ymin = -Inf, ymax = 1, fill = "dodgerblue", alpha = 0.5) +
  annotate("rect", xmin = min(lake$day.date), xmax = max(lake$day.date), ymin = -Inf, ymax = 1, fill = "grey80", alpha = 0.5) +
  annotate("rect", xmin = min(spawn$day.date), xmax = max(spawn$day.date), ymin = -Inf, ymax = 1, fill = "forestgreen", alpha = 0.5) +
  annotate("text", x = two.years[618-233,]$day.date, y =1.28, label = "Lake Phase", size = 2 ) +
  annotate("text", x = two.years[122-15,]$day.date, y =1.28, label = "Emergence and Outmigration\n(May 3 - June 2)", size = 2 ) +
  annotate("text", x = two.years[648+15,]$day.date, y =1.28, label = "Spawning\n(Sept 11 - Oct 10)", size = 2 ) +
  geom_line(data = two.years, aes(x = day.date, y= daylight.diff), size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_date(date_labels = "%b") +
  lims(y = c(-1.1,1.35)) +
  labs(x = " Date", y = "Hours of light difference\nin Great Lakes", title = "c)") +
  theme_classic()

#ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/Fig 5.pdf", fig.5, height = 2.25, width = 4.5, units = "in", dpi = 600)


# Plot per2 allele freq ---------------------------------------------------
per2.AF <- read_tsv("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/per2/per2.allele.freq.txt")
per2.pval <- read.table("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/mark_model_results/p-vals.txt", row.names = NULL, sep = "\t", header = TRUE)
per2.pval.05 <- per2.pval %>% filter(pval <= 0.05 & A1.freq != 0)

per2.obs <- ggplot(data = per2.AF) +
  geom_line(aes(x = pop, y = A1.freq, group = pos), size = 0.25, color = "darkgrey", linetype = "dashed") +
  geom_line(data = per2.AF[per2.AF$pos %in% per2.pval.05$pos, ], aes(x = pop, y = A1.freq, group = pos), size = 0.25, color = "darkgrey") +
  geom_line(data = per2.AF[per2.AF$pos == 14729846,], aes(x = pop, y = A1.freq, group = pos), size = 0.5, color = "red") +
  geom_point(aes(x = pop, y = A1.freq, color = pop, group = pos), size = 0.75) +
  labs(x = NULL, y= NULL, title = NULL,subtitle = "Observed") +
  scale_color_manual(values = c("BC Odd" = "blue4","GL Odd" = "darkorange2"), ) +
  theme_classic() +
  theme(legend.position = "none")

per2.obs.2 <- ggplot(data = per2.AF) +
  geom_line(aes(x = pop, y = A1.freq, group = pos), size = 0.25, color = "darkorchid2", linetype = "dashed") +
  geom_line(data = per2.AF[per2.AF$pos %in% per2.pval.05$pos, ], aes(x = pop, y = A1.freq, group = pos), size = 0.25, color = "darkorchid2") +
  geom_line(data = per2.AF[per2.AF$pos == 14729846,], aes(x = pop, y = A1.freq, group = pos), size = 0.5, color = "red") +
  geom_point(aes(x = pop, y = A1.freq, color = pop, group = pos), size = 0.75) +
  labs(x = NULL, y= NULL, title = NULL,subtitle = "Observed") +
  scale_color_manual(values = c("BC Odd" = "blue4","GL Odd" = "darkorange2"), ) +
  theme_classic() +
  theme(legend.position = "none")


hi <- per2.AF %>%  
  filter(pop == "GL Odd" & A1.freq >= 0.5) %>% 
  dplyr::select(A2.freq) %>% 
  rename(A2.freq = "freq")

lo <- per2.AF %>%  
  filter(pop == "GL Odd" & A1.freq <= 0.5) %>% 
  dplyr::select(A1.freq) %>% 
  rename(A1.freq = "freq")
  
all <- rbind(hi,lo)

mean(all$freq)

### find SNP used in simulation
per2.AF[per2.AF$pop == "BC Odd" & per2.AF$A1.freq <= 0.034, ] 
# it's SNP 14723696

### look at AF in both
per2.AF[per2.AF$pos ==14723696,]

### look at AF in GL

# Plot sim per2 -----------------------------------------------------------

#per2.sim.AF <- read.table("./data/seqs/aligned_reads_OgorEven_v1.0/10_data_analysis/pop_gen/Fst/outliers/per2/Panel_B.txt", row.names = NULL, sep = "\t", header = TRUE)
GL.sim<- read.table("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/mark_model_results/GL_95%CIs.txt", row.names = NULL, sep = "\t", header = TRUE)
BC.sim<- read.table("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/mark_model_results/BC_95%CIs.txt", row.names = NULL, sep = "\t", header = TRUE)
per2.sim.AF <- GL.sim

# make objs of AF >=0.5 and < 0.5
per2.sim.highAF <- per2.sim.AF %>% filter(A1.freq >= 0.5)
per2.sim.lowAF <- per2.sim.AF %>% filter(A1.freq < 0.5)

#create seperate df for BC values and GL values above and below 0.5 to rbind later
BC.sim.AF <- per2.sim.AF[,c("pos","A1.freq", "pop")]

GL.sim.high.AF <- per2.sim.highAF[, c("pos","low95", "pop")]
GL.sim.high.AF$pop <- "GL Odd"
colnames(GL.sim.high.AF)[2] <- "A1.freq"

GL.sim.low.AF <- per2.sim.lowAF[, c("pos","up95", "pop")]
GL.sim.low.AF$pop <- "GL Odd"
colnames(GL.sim.low.AF)[2] <- "A1.freq"

sim.AF <- rbind(BC.sim.AF, GL.sim.high.AF, GL.sim.low.AF)

# calc avg change in allele freq
per2.sim.AF %>% 
  mutate(delta.AF = abs(A1.freq - sim.AF$A1.freq[39:76])) %>% 
  summarise(mean(delta.AF))

per2.sim <- ggplot(data = sim.AF) +
  geom_line(aes(x = pop, y = A1.freq, group = pos), size = 0.25, color = "darkgrey", linetype = "dashed") +
  geom_line(data = sim.AF[sim.AF$pos == 14729846,], aes(x = pop, y = A1.freq, group = pos), size = 0.5, color = "red", linetype = "dashed") +
  geom_point(aes(x = pop, y = A1.freq, color = pop, group = pos), size = 0.75) +
  labs(x = NULL, y= "Reference allele frequency", title = "b)",subtitle = "Simulated") +
  scale_color_manual(values = c("BC Odd" = "blue4","GL Odd" = "darkorange2"), ) +
  theme_classic() +
  theme(legend.position = "none")

per2.sim.2 <- ggplot(data = sim.AF) +
  geom_line(aes(x = pop, y = A1.freq, group = pos), size = 0.25, color = "darkorchid2", linetype = "dashed") +
  geom_line(data = sim.AF[sim.AF$pos == 14729846,], aes(x = pop, y = A1.freq, group = pos), size = 0.5, color = "red", linetype = "dashed") +
  geom_point(aes(x = pop, y = A1.freq, color = pop, group = pos), size = 0.75) +
  labs(x = NULL, y= "Reference allele frequency", title = "b)",subtitle = "Simulated") +
  scale_color_manual(values = c("BC Odd" = "blue4","GL Odd" = "darkorange2"), ) +
  theme_classic() +
  theme(legend.position = "none")


# Allele trajectories -----------------------------------------------------

drift.sim <- read.table("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/mark_model_results/Panel_A_0.034.txt", header = TRUE)
drift.sim$index <- rep(c(1:1000),each = 254)
#drift.sim <- drift.sim[drift.sim$year %% 2 != 0, ]
drift.sim <- na.omit(drift.sim)

pval.95 <-stats::quantile(drift.sim[drift.sim$year == 254, "allelefreq"], probs =c(0.95), na.rm = TRUE)

nrow(drift.sim[drift.sim$year == 254 & drift.sim$allelefreq == 0,])/nrow(drift.sim[drift.sim$year == 1 ,])
nrow(drift.sim[drift.sim$year == 254 & drift.sim$allelefreq >= 0.966,])/nrow(drift.sim[drift.sim$year == 1 ,])

# traj <- ggplot() +
#   geom_line(data = drift.sim %>% filter(year <= 136 ), aes(x = abs(year-200), y = allelefreq, group = index), color = "blue4", linewidth =.2) +
#   geom_line(data = drift.sim %>% filter(year >= 136 ), aes(x = abs(year-200), y = allelefreq, group = index), color = "darkorange2", linewidth =.2) +
#   geom_hline(yintercept = pval.95, linetype = "dashed" ) +
#   labs(x = "Generations ago", y ="Allele frequency") +
#   theme_classic()

## need to have these objecs from plot_GONE.R script loaded so they're subset at 127 generations.
GONE.dat <- rbind(LAO_Ne[LAO_Ne$Generation >31,],STLO_Ne[STLE_Ne$Generation <=31,] )

dat <- read.table("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/gl_GONE.summary.data.txt", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
dat2 <- dat[which(dat[, 5] == "LAO"), ] # BC
bc   <- dat2[26:100, ]    # FOR METHODS TOOK GENS 26 to 100
dat2 <- dat[which(dat[, 5] == "STLO"), ] # BC
gl   <- dat2[1:52, ] # FOR METHODS TOOK GENS 1:52
both    <- rbind(gl, bc)
both[1] <- 1:nrow(both)
dat2    <- both
dat3 <- dat2[order(dat2[, 1], decreasing = TRUE), ]

GONE <- ggplot() +
  geom_line(data = GONE.dat, aes(x =Generation, y = Ne), stat = "summary", fun ="mean", linewidth = 0.75) +
  labs(y  =  bquote(paste(italic('N'['e']*' '))),  title = "a)") +
  #scale_y_continuous(trans='log10') +
  theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())

GONE.2 <- ggplot() +
  geom_line(data = GONE.dat, aes(x =abs(Generation-127), y = Ne), stat = "summary", fun ="mean", linewidth = 0.75) +
  labs(y  =  bquote(paste(italic('N'['e']*' '))),  title = "a)") +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank(),
                                        axis.line.x = element_blank())

GONE.2.abm <- ggplot() +
  geom_line(data = dat3, aes(x = abs(Generation-127), y = mean.ne), linewidth = 0.75) +
  labs(y  =  bquote(paste(italic('N'['e']*' '))),  title = "a)") +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())

traj <- ggplot(drift.sim) +
  geom_line(aes(x = abs(year-254)/2, y = allelefreq, group = index), color = "", linewidth =.1) +
  geom_point(aes(x = 0, y = 0.933), color = "darkorange2", size = 0.5) +
  geom_point(aes(x = 127, y = 0.033), color = "blue4", size = 0.5) +
  geom_hline(yintercept = pval.95, linetype = "dashed", linewidth = 0.5 ) +
  labs(x = "Generations ago", y ="Allele frequency") +
  theme_classic()

dense <-ggplot(drift.sim[drift.sim$year == 254,], aes(x = year, y = allelefreq )) +
  stat_slab(orientation = "vertical", fill = "darkgrey", slab_type = 'pdf' ) + 
  geom_hline(yintercept = pval.95, linetype = "dashed", linewidth = 0.5) +
  theme_void()


traj.2 <- ggplot(drift.sim[drift.sim$year%%5 == 0 ,]) +
  geom_line(aes(x = year/2, y = allelefreq, group = index), color = "darkorchid2", linewidth =.1) +
  geom_point(aes(x = 127, y = 0.933), color = "darkorange2", size = 0.5) +
  geom_point(aes(x = 0, y = 0.033), color = "blue4", size = 0.5) +
  geom_hline(yintercept = pval.95, linetype = "dashed", linewidth = 0.5 ) +
  lims(y = c(0,1.01)) +
  labs(x = "Generations", y ="Allele frequency") +
  theme_classic()

its <-  sample(drift.sim$index, replace = FALSE, size = 200)
traj.2 <- ggplot(drift.sim) +
  geom_line(data = drift.sim[drift.sim$index %in% its,], aes(x = year/2, y = allelefreq, group = index), color = "darkorchid2", linewidth =.1) +
  geom_point(aes(x = 127, y = 0.933), color = "darkorange2", size = 0.5) +
  geom_point(aes(x = 0, y = 0.033), color = "blue4", size = 0.5) +
  geom_hline(yintercept = pval.95, linetype = "dashed", linewidth = 0.5 ) +
  lims(y = c(0,1.01)) +
  labs(x = "Generations", y ="Allele frequency") +
  theme_classic()

#ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/Fig5.traj.png", traj.2, height = 1.8, width = 3.6, units = "in", dpi = 1200)

dense.2 <-ggplot(drift.sim[drift.sim$year == 254,], aes(x = year, y = allelefreq )) +
  stat_slab(orientation = "vertical", fill = "darkorchid2", slab_type = 'pdf' ) + 
  lims(y = c(0,1)) +
  theme_void()

# ggplot() +
# geom_density(data = drift.sim[drift.sim$year == 200,], aes(y = allelefreq))+
#   scale_y_continuous(breaks = seq(0,1,.01)) +
#   theme_void()

small_GONE / traj


layout <- "
AAAAAA#
AAAAAA#
BBBBBBC
BBBBBBC
BBBBBBC
BBBBBBC
BBBBBBC
"
blank.gg <- ggplot() + theme_void()

fig.5.top <-GONE + traj + dense + plot_layout(design = layout)

fig.5.top.2 <-GONE.2.abm + traj.2 + dense.2 + plot_layout(design = layout)
fig.5.top.2

#sims that fixed
nrow(drift.sim[drift.sim$allelefreq == 0,])/nrow(drift.sim)

# starting allele freq
mean(drift.sim[drift.sim$year== 1,"allelefreq"], na.rm = TRUE)
  

# Combine plots -----------------------------------------------------------



fig.5 <-fig.5.top / (per2.sim | per2.obs) / daylight

fig.5.2 <- fig.5.top.2 / (per2.sim.2 | per2.obs.2) / daylight +  plot_annotation(title = 'Figure 5')


ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/Fig 5.pdf", fig.5.2, height = 8.5, width = 4.5, units = "in", dpi = 100)


# SI Plots ----------------------------------------------------------------

### plot of sim points and 95% CI

BC.sim$pop <- "BC Odd"
GL.sim$pop <- "GL Odd"

per2.sim.SI <- ggplot(per2.AF,aes(x = pop, y = A1.freq, color = pop)) +
  geom_point(data = BC.sim, aes (x = pop, y = allele),  color = "darkgrey") +
  geom_point(data = GL.sim, aes (x = pop, y = allele),  color = "darkgrey") +
  geom_errorbar(data = GL.sim, aes(x = pop, ymin = low95, ymax = up95), color = "darkgrey") +
  geom_errorbar(data = BC.sim, aes(x = pop, ymin = low95, ymax = up95), color = "darkgrey") +
  geom_point() +
  geom_text(data = per2.pval.05,aes(x =2.25, y= .05), label = "*", color = "red", size =8) +
  facet_wrap(~pos) +
  scale_color_manual(values = c("BC Odd" = "blue4","GL Odd" = "darkorange2")) +
  labs(title = "SNP position on Chromosome 22", x = NULL, y = "Reference allele frequency") +
  #coord_cartesian(y=c(0,1), clip = "off") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("~/Dropbox/PhD Work/GL Pink Salmon/Manuscripts/ CH 2/Figures/per2.sim.SI.pdf", per2.sim.SI, height = 8, width = 8.5, units = "in", dpi = 600)





