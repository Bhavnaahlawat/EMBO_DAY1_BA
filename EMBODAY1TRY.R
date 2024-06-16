x<-5
y<-2

mean(x,y)

#trial 123


#trying for Day 2 at EMBO POPGEN

#Genetic Drift 

#new exercise =1

# N= population size
#fA Allel frequency 
#rbinom (n, size, frequecnyc)/size
# rbinom (1, 2*N, 0.5)/2*N
#n= generation 


N <- 100
fA <- 0.5

rbinom(1,2*N,0.5)/(2*N)

#frequency of next generation = 0.525


rbinom(1,2*N,0.525)/(2*N)

#frequency of next generation = 0.48

rbinom(1,2*N,0.48)/(2*N)

#frequency of next generation = 0.485

rbinom(1,2*N,0.485)/(2*N)

#frequency of next generation = 0.46


# I did for multiple generation with different frequencies got at the end of each answer
# now i will take the frequency constant but change the number of generation 

N <- 100
fA <- 0.5

rbinom(1,2*N,0.5)/(2*N)

#frequency of next generation = 0.525

rbinom(2,2*N,0.5)/(2*N)

#frequency of 2 generation = 0.535, 0.510

rbinom(3,2*N,0.5)/(2*N)

#frequency of 3 generation = 0.505, 0.495, 0.460

rbinom(4,2*N,0.5)/(2*N)

#frequency of 4 generation = 0.490 0.520 0.555 0.480

rbinom(100,2*N,0.5)/(2*N)
#frequency of 4 generation = 0.485 0.495 0.570 0.505 0.510 0.475 0.505 0.460 0.485 0.490 0.530 0.500 0.545 0.505 0.435 0.530 0.520 0.475 0.570 0.565 0.445 0.450 0.475 0.540 0.490 0.500 0.490 0.430 0.510 0.515 0.555 0.500 0.480 0.565 0.545 0.480 0.450 0.490 0.535 0.520 0.500 0.520 0.485 0.530 0.525 0.495 0.450 0.520 0.570 0.485 0.520 0.480 0.500 0.495 0.505 0.510 0.485 0.445 0.475 0.435 0.515 0.510 0.525 0.460 0.490 0.545 0.550 0.520 0.505 0.485 0.530 0.500 0.540 0.480 0.520 0.485 0.495 0.460 0.510 0.415 0.515 0.450 0.540 0.530 0.475

frequency<- rbinom(100,2*N,0.5)/(2*N)

print(frequency)
plot(frequency)

# new exercise=2

fA <-rep(NA,100)

# at t=0

fA [1] <- 0.50
for (t in 1:99) fA [t+1] <- rbinom(1,2*N, fA [t])/(2*N)

plot (x=1:100, y=fA, type = "l", ylim = c(0,1), lwd=2)


# did this for multiple generation with frequencies from previous generation

# new exercise =3

#two different populations
#populations small=red
#population large= blue 

NR <- 60

fA <-rep(NA,100)

# at t=0

fA [1] <- 0.50
for (t in 1:99) fA [t+1] <- rbinom(1,2*NR, fA [t])/(2*NR)

plot (x=1:100, y=fA, type = "l", ylim = c(0,1), lwd=2)

NB <- 100000
fA <-rep(NA,100)

# at t=0

fA [1] <- 0.50
for (t in 1:99) fA [t+1] <- rbinom(1,2*NB, fA [t])/(2*NB)


plot (x=1:100, y=fA, type = "l", ylim = c(0,1), lwd=2)


#using velentias code

# large population
N = 500

# initial allele frequency for allele A
fA = 0.5

# number of generations to simulate
n_gen = 100

# table to store the allele frequency over time
large = data.frame(Generation = 1:n_gen, Allele_frequency = fA)

for(i in 1:n_gen){
  
  fA = rbinom(1, 2*N, fA) / (2*N)
  large[i, 2] = fA 
  
}

# small population
N = 10

# table to store the allele frequency over time
small = data.frame(Generation = 1:n_gen, Allele_frequency = fA)

for(i in 1:n_gen){
  
  fA = rbinom(1, 2*N, fA) / (2*N)
  small[i, 2] = fA 
  
}
small$Type = "small"
large$Type = "large"
data = rbind(large, small)

ggplot(data = data, aes(x = Generation, y = Allele_frequency, col = Type)) + geom_line() + geom_point() + theme_bw() + ylab("Allele frequency") + xlab("Generation") + ggtitle("Initial allele frequency of 0.5")

#calculations 
1/(1+1/2+1/3)

#example file by Matteo

import msprime

# Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.
ts = msprime.sim_ancestry(
  samples=10,
  recombination_rate=1e-8, # as in humans
  sequence_length=10_000,
  population_size=10_000,
  random_seed=1234)

# we can add mutations
mutated_ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=1234)
print(mutated_ts.tables.sites)

for variant in mutated_ts.variants():
  print(variant)

#no didnt work had to do in python, so ran in embo serve, worked fine. 

#Day 3

install.packages("reticulate")

msprime <- reticulate::import("msprime")

#trying task 2 day 3 afternon session 

sims <- read.csv("mosquito-task2.csv", head=T)

# check prior distributions
x11()
par(mfrow=c(2,2))
hist(sims$N1)
hist(sims$N2)
hist(sims$T_split)
hist(sims$MigRate)

# remove simulations with NaN for some summary stats!

# find useful summary stats which correlate with T_split
cor(sims$Fst, sims$T_split)
cor(sims$dxy, sims$T_split)
cor(sims$segsites1, sims$T_split)
cor(sims$segsites2, sims$T_split)
cor(sims$pi1, sims$T_split)
cor(sims$pi2, sims$T_split)
cor(sims$tajima1, sims$T_split)
cor(sims$tajima2, sims$T_split)

# load observed summary stats
obs <- read.csv("mosquito-observed.csv", head=T)

# check if simulated retained summary stats contain the observed one
quantile(sims$Fst); cat(obs$Fst)
quantile(sims$segsites1); cat(obs$segsites1)
quantile(sims$segsites2); cat(obs$segsites2)

# merge obs with retained sims to scale them
sumstats <- scale(rbind(obs[c(1,3,4)],sims[,c(5,7,8)]))

library(abc)

est <- abc(target=sumstats[1,], param=sims$T_split, sumstat=sumstats[-1,], tol=0.05, method="rejection")

# check distances in the acceptance region
hist(est$dist)
abline(v=max(est$dist[which(est$region)]), lty=2)

# posterior
x11()
par(mfrow=c(2,1))
hist(est$unadj.values, freq=FALSE, xlim=range(sims$T_split), col=rgb(0,0,1,1/4), main="Posterior probability", xlab="Split time")
# MAP
map <- mean(est$unadj.values)
abline(v=map, lty=2)

# confidence intervals
hpd <- quantile(x=est$unadj.values, probs=c(0.025,0.975))
abline(v=hpd, lty=3)

# prior
hist(sims$T_split, freq=FALSE, xlim=range(sims$T_split), col=rgb(1,0,0,1/4))



#day 4

library(relater)
library(ggplot2)

install.packages(relater)

test1

#Day 5
#installing packages 


install.packages("tidypopgen",
                 repos = c("https://evolecolgroup.r-universe.dev",
                           "https://cloud.r-project.org"))

install.packages("admixtools",
                 repos = c("https://evolecolgroup.r-universe.dev", "https://cloud.r-project.org"))

?tidypopgen::annotate_group_info()
library(tidypopgen)

library(admixtools) 
library(ggplot2)

modern_gt <- tidypopgen::gen_tibble("./modern_samples.bed", valid_alleles = c("A","T","C","G"),
                                    missing_alleles = c("X"))

modern_gt

modern_gt %>% group_by(population) %>% tally()

loci_report <- qc_report_loci(modern_gt)
autoplot(loci_report, type = "missing")


modern_gt <- modern_gt %>% select_loci_if(loci_missingness(genotypes)<0.04)

ancient_gt <- tidypopgen::gen_tibble("./ancient_samples.vcf",
                                     valid_alleles = c("A","T","C","G","X"),
                                     quiet = TRUE)

install.packages("vcfR")
library(vcfR)

ancient_gt$id

ancient_gt$id[ancient_gt$id == "GB20"] <- "Mota"
ancient_gt$population <- ancient_gt$id


merged_dry <- rbind_dry_run(modern_gt, ancient_gt, flip_strand = TRUE)

merged_gt <- rbind(modern_gt, ancient_gt, flip_strand = TRUE,
                   backingfile = "./merged_samples")

# group test_gt by population
merged_gt <- merged_gt %>% group_by(population)


f2_dir <- "./f2_tidypopgen"

# compute f2
f2_tidypopgen <- gt_extract_f2(merged_gt,
                               outdir = "./f2_tidypopgen",
                               overwrite = TRUE)

f2_blocks = f2_from_precomp("./f2_tidypopgen")

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan", "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")

lbk_f3out <- f3(data = f2_blocks, pop1 = "Mbuti",
                pop2 = "LBK",
                pop3 = lbk_modern_panel)

lbk_f3out

lbk_f3out %>% arrange(desc(est))


ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) + labs(y = "Shared drift with LBK", x = "populations") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


lbk_f3out$pop3<-factor(lbk_f3out$pop3, levels = lbk_f3out$pop3[order(lbk_f3out$est)])

ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) + labs(y = "Shared drift with LBK", x = "populations") + theme(axis.text.x = element_text(angle = 90, hjust = 1))



aa_f3admix <- f3(data = f2_blocks, pop1 = "AA",
                 pop2 = "Yoruba",
                 pop3 = "French")

aa_f3admix

eurasian_sources <- c("French","Spanish","Sardinian","LBK")

somali_f3admix <- f3(data = f2_blocks, pop1 = "Somali",
                     pop2 = eurasian_sources,
                     pop3 = "Mota")


somali_f3admix

neand_french_f4 <- f4(data = f2_blocks,
                      pop1 = "pan_troglodytes", 
                      pop2 = "AltaiNea",
                      pop3 = "Mbuti",
                      pop4 = "French")
neand_french_f4



pops <- c("Han", "pan_troglodytes", "AA","Yoruba","French") 

qpf4ratio(data = f2_blocks,
          pops = pops)

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan", "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")
modern_panel_gt <- merged_gt %>% filter (population %in% lbk_modern_panel)

modern_panel_gt <- modern_panel_gt %>% ungroup() %>% select_loci_if(loci_maf(genotypes)>0)


modern_panel_gt <- modern_panel_gt %>% gt_impute_simple(method="mode")


modern_panel_gt <- modern_panel_gt %>% select_loci_if(loci_ld_clump(genotypes))

modern_pca <- modern_panel_gt %>% gt_pca_randomSVD()

autoplot(modern_pca, type = "scores")


library(ggplot2)
autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) + labs(color = "population")



lbk_gt <- merged_gt %>% filter(id == "LBK")
lbk_pca_scores <- predict(modern_pca, new_data = lbk_gt, project_method = "least_square")
lbk_pca_scores



autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")+
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black")


autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population") +
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black") + lims(x=c(30, 70), y = c(-10, 15))



#second pratical 

library(admixtools) 
library(tidypopgen)


f2_blocks = f2_from_precomp("./f2_tidypopgen", verbose = FALSE)

neand_euras_wave <- qpwave(data = f2_blocks,
                           left = c("French","Spanish","Tuscan"),
                           right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)


neand_euras_wave

french_adm <- qpadm(data = f2_blocks,
                    left = c("Loschbour", "LBK", "Yamnaya"),
                    right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"), target= "French")

french_adm$popdrop


base_edges <- matrix( c("R", "Mbuti",
                        "R", "eAfr",
                        "eAfr", "Dinka","eAfr", "outAfrica",
                        "outAfrica",    "Han",
                        "outAfrica",    "Loschbour"),
                      ncol=2,
                      byrow = TRUE,
                      dimnames=list(NULL, c("from","to")))

base_edges


base_igraph <- base_edges %>% edges_to_igraph()


is_valid(base_igraph)

base_igraph %>% plot_graph()

base_igraph %>% plotly_graph()

install.packages("plotly")

library(plotly)

base_qpgraph <- qpgraph(data = f2_blocks, graph = base_igraph)

base_qpgraph$f3

base_qpgraph$f3 %>% filter(abs(z)>2)

base_qpgraph$edges %>% plot_graph()

fits = qpgraph_resample_multi(f2_blocks,
                              graphlist = list(base_qpgraph[[1]], base_swapped_qpgraph[[1]]), 
                              nboot = 100) 
                              

compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
                              

base_igraph %>% plot_graph(highlight_unidentifiable = TRUE)                  


lbk_extra_edges <- matrix( 
  c(
  "R",    "Mbuti",
  "R",    "eAfr",
  "eAfr", "pBasalEurasian",
  "eAfr", "Dinka",
  "pBasalEurasian", "BasalEurasian",
  "pBasalEurasian","outAfrica",
  "outAfrica", "Han",
  "outAfrica","wEurasian",
  "wEurasian", "Yamnaya",
  "wEurasian", "pLoschbour",
  "pLoschbour", "Loschbour",
  "pLoschbour","WHG",
  "BasalEurasian", "pLBK",
  "WHG", "pLBK",
  "pLBK","LBK"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

is_valid(lbk_extra_igraph)
lbk_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

lbk_extra_qpgraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)
lbk_extra_qpgraph$edges %>% plot_graph()



#Day 6 pratical by Tabitia 


# Reading a FST file
Pop1_fst_data <- read.table("./AFR_EAS.weir.fst")
Pop2_fst_data <- read.table("./AFR_EUR.weir.fst")
Pop3_fst_data <- read.table("./EAS_EUR.weir.fst")

#Eliminate duplicate positions

install.packages("dplyr")
library(dplyr)

# Remove duplicate rows based on the Position column using base R
Pop1_fst_data_Nodups <- Pop1_fst_data [!duplicated(Pop1_fst_data), ]
Pop2_fst_data_Nodups <- Pop2_fst_data [!duplicated(Pop2_fst_data), ]
Pop3_fst_data_Nodups <- Pop3_fst_data [!duplicated(Pop3_fst_data), ]

# Display the first few rows of the data to verify duplicates are removed
head(Pop1_fst_data_Nodups)
head(Pop2_fst_data_Nodups)
head(Pop3_fst_data_Nodups)


# Assuming your data frame is named fst_data, exclude values with NA values 
clean_Pop1_fst_data_Nodups <- na.omit(Pop1_fst_data_Nodups)
clean_Pop2_fst_data_Nodups <- na.omit(Pop2_fst_data_Nodups)
clean_Pop3_fst_data_Nodups <- na.omit(Pop3_fst_data_Nodups)

# since there were no duplicates we can also use the initial file as well

clean_Pop1_fst_data <- na.omit(Pop1_fst_data)
clean_Pop2_fst_data <- na.omit(Pop2_fst_data)
clean_Pop3_fst_data <- na.omit(Pop3_fst_data)


# Display the first few rows to verify NAs are removed
head(clean_Pop1_fst_data)
head(clean_Pop2_fst_data)
head(clean_Pop3_fst_data)


#filter the file by retaining only  the common SNPs

commonSNPs_Pop1_Pop3 <- intersect(clean_Pop1_fst_data$POS,clean_Pop3_fst_data$POS)
commonSNPs_Pop1_Pop3_dataframe <- tibble(POS= commonSNPs_Pop1_Pop3)

Common_SNPs_to_retain <- intersect(commonSNPs_Pop1_Pop3_dataframe, clean_Pop2_fst_data$POS)



typeof(commonSNPs_Pop1_Pop2)

# this last step did not work out 
# trying with code from Letizia 

#Tàbita_Selection
#Read the files with the Fst estimates (AFR_EUR.weir.fst, AFR_EAS.weir.fst and EAS_EUR.weir.fst)
setwd("/Users/eugenia/Desktop/EMBO_Practical_Course_2024-main")
Pop1_fst_data <- read.table(file = 'AFR_EUR.weir.fst', header = T)
Pop2_fst_data <- read.table(file = 'AFR_EAS.weir.fst', header = T)
Pop3_fst_data <- read.table(file = 'EAS_EUR.weir.fst', header = T)
#Eliminate duplicate positions
Pop1_fst_data_noDup <- Pop1_fst_data[!duplicated(Pop1_fst_data), ]
Pop2_fst_data_noDup <- Pop2_fst_data[!duplicated(Pop2_fst_data), ]
Pop3_fst_data_noDup <- Pop3_fst_data[!duplicated(Pop3_fst_data), ]
#Exclude positions whose FST result was equal to Na
Pop1_fst_data_noDup_noNA <- Pop1_fst_data_noDup[!is.na(Pop1_fst_data_noDup$WEIR_AND_COCKERHAM_FST), ]
Pop2_fst_data_noDup_noNA <- Pop2_fst_data_noDup[!is.na(Pop2_fst_data_noDup$WEIR_AND_COCKERHAM_FST), ]
Pop3_fst_data_noDup_noNA <- Pop3_fst_data_noDup[!is.na(Pop3_fst_data_noDup$WEIR_AND_COCKERHAM_FST), ]
# Step 1: Find intersection between AFR_EUR, AFR_EAS positions
AFR_EUR_AFR_EAS <- intersect(Pop1_fst_data_noDup_noNA$POS, Pop2_fst_data_noDup_noNA$POS)
AFR_EUR_AFR_EAS_df <- data.frame(POS = AFR_EUR_AFR_EAS)
AFR_EUR_AFR_EAS_EAS_EUR <- intersect(AFR_EUR_AFR_EAS, Pop3_fst_data_noDup_noNA$POS)
# Step 4: Convert the final result to a data frame
AFR_EUR_AFR_EAS_EAS_EUR_df <- data.frame(POS = AFR_EUR_AFR_EAS_EAS_EUR)
AFR_EUR.weir.fst.common <- AFR_EUR.weir.fst.noDup.noNA[AFR_EUR.weir.fst.noDup.noNA$POS %in% AFR_EUR_AFR_EAS_EAS_EUR_df$POS, ]
AFR_EAS.weir.fst.common <- AFR_EAS.weir.fst.noDup.noNA[AFR_EAS.weir.fst.noDup.noNA$POS %in% AFR_EUR_AFR_EAS_EAS_EUR_df$POS, ]
EAS_EUR.weir.fst.common <- EAS_EUR.weir.fst.noDup.noNA[EAS_EUR.weir.fst.noDup.noNA$POS %in% AFR_EUR_AFR_EAS_EAS_EUR_df$POS, ]


#Letizias code which worked 
#Read the files with the Fst estimates (AFR_EUR.weir.fst, AFR_EAS.weir.fst and EAS_EUR.weir.fst)
setwd("/Users/eugenia/Desktop/EMBO_Practical_Course_2024-main")
AFR_EUR.weir.fst <- read.table(file = 'AFR_EUR.weir.fst', header = T)
AFR_EAS.weir.fst <- read.table(file = 'AFR_EAS.weir.fst', header = T)
EAS_EUR.weir.fst <- read.table(file = 'EAS_EUR.weir.fst', header = T)
#Eliminate duplicate positions
AFR_EUR.weir.fst.noDup <- AFR_EUR.weir.fst[!duplicated(AFR_EUR.weir.fst), ]
AFR_EAS.weir.fst.noDup <- AFR_EAS.weir.fst[!duplicated(AFR_EAS.weir.fst), ]
EAS_EUR.weir.fst.noDup <- EAS_EUR.weir.fst[!duplicated(EAS_EUR.weir.fst), ]
#Exclude positions whose FST result was equal to Na
AFR_EUR.weir.fst.noDup.noNA <- AFR_EUR.weir.fst.noDup[!is.na(AFR_EUR.weir.fst.noDup$WEIR_AND_COCKERHAM_FST), ]
AFR_EAS.weir.fst.noDup.noNA <- AFR_EAS.weir.fst.noDup[!is.na(AFR_EAS.weir.fst.noDup$WEIR_AND_COCKERHAM_FST), ]
EAS_EUR.weir.fst.noDup.noNA <- EAS_EUR.weir.fst.noDup[!is.na(EAS_EUR.weir.fst.noDup$WEIR_AND_COCKERHAM_FST), ]
# Step 1: Find intersection between AFR_EUR, AFR_EAS positions
AFR_EUR_AFR_EAS <- intersect(AFR_EUR.weir.fst.noDup.noNA$POS, AFR_EAS.weir.fst.noDup.noNA$POS)
AFR_EUR_AFR_EAS_df <- data.frame(POS = AFR_EUR_AFR_EAS)
AFR_EUR_AFR_EAS_EAS_EUR <- intersect(AFR_EUR_AFR_EAS, EAS_EUR.weir.fst.noDup.noNA$POS)
# Step 4: Convert the final result to a data frame
AFR_EUR_AFR_EAS_EAS_EUR_df <- data.frame(POS = AFR_EUR_AFR_EAS_EAS_EUR)
AFR_EUR.weir.fst.common <- AFR_EUR.weir.fst.noDup.noNA[AFR_EUR.weir.fst.noDup.noNA$POS %in% AFR_EUR_AFR_EAS_EAS_EUR_df$POS, ]
AFR_EAS.weir.fst.common <- AFR_EAS.weir.fst.noDup.noNA[AFR_EAS.weir.fst.noDup.noNA$POS %in% AFR_EUR_AFR_EAS_EAS_EUR_df$POS, ]
EAS_EUR.weir.fst.common <- EAS_EUR.weir.fst.noDup.noNA[EAS_EUR.weir.fst.noDup.noNA$POS %in% AFR_EUR_AFR_EAS_EAS_EUR_df$POS, ]


#step 5: Convert positions with estimates from FST < 0 to = 0

AFR_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST <- ifelse(AFR_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST < 0, 0, AFR_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST)
min(AFR_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST)
AFR_EAS.weir.fst.common$WEIR_AND_COCKERHAM_FST <- ifelse(AFR_EAS.weir.fst.common$WEIR_AND_COCKERHAM_FST < 0, 0, AFR_EAS.weir.fst.common$WEIR_AND_COCKERHAM_FST)
min(AFR_EAS.weir.fst.common$WEIR_AND_COCKERHAM_FST)
EAS_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST <- ifelse(EAS_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST < 0, 0, EAS_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST)
min(EAS_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST)


#Step 6: Given value x (SNP position)
x <- 109513601
calculate_window <- function(x) {
  lower_bound <- x - 5000
  upper_bound <- x + 5000
  return(c(lower_bound, upper_bound))
}
bounds <- calculate_window(x)
AFR_EUR_10kb <- AFR_EUR.weir.fst.common %>%
  filter(POS >= bounds[1], POS <= bounds[2])
AFR_EAS_10kb <- AFR_EAS.weir.fst.common %>%
  filter(POS >= bounds[1], POS <= bounds[2])
EAS_EUR_10kb <- EAS_EUR.weir.fst.common %>%
  filter(POS >= bounds[1], POS <= bounds[2])

#Step 7: Plot FST values 

#trying histogram first
pdf("FST_plot_AFR_EUR_10kb.pdf")
hist(AFR_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST, breaks =100, xlab="fst", ylab = "frequecy", col = "blue")
dev.off()
pdf("FST_plot_AFR_EAS_10kb.pdf")
hist(AFR_EAS.weir.fst.common$WEIR_AND_COCKERHAM_FST, breaks =100, xlab="fst", ylab = "frequecy", col = "Red")
dev.off()
pdf("FST_plot_EAS_EUR_10kb.pdf")
hist(EAS_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST, breaks =100, xlab="fst", ylab = "frequecy", col = "Green")
dev.off()




#trying ggplot second
library(ggplot2)
pdf("FST_plot_AFR_EUR_2_10kb.pdf")

ggplot(data = AFR_EUR.weir.fst.common$WEIR_AND_COCKERHAM_FST, aes(x = fst, y = frequency, col = blue)) + geom_line() + geom_point() + theme_bw() + ylab("frequency") + xlab("fst") + ggtitle("FST_plot_AFR_EUR")
dev.off()





#trying PBS 


# Load necessary libraries


# Generate some synthetic SNP data for three populations
set.seed(123)
n <- 100  # Number of individuals per population
m <- 1000  # Number of SNPs

EAS <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n)
EUR <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n)
AFR <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n)

# Combine into a genind object
pop <- factor(rep(c("EAS", "EUR", "AFR"), each = n))
genind_obj <- df2genind(cbind(EAS, EUR, AFR), pop = pop, ploidy = 2, type = "codom")

# Calculate pairwise Fst values
fst_matrix <- pairwise.fst(genind_obj)
print(fst_matrix)

# Define a function to calculate PBS values
calculate_PBS <- function(fst1, fst2, fst3) {
  T1 <- -log(1 - fst1)
  T2 <- -log(1 - fst2)
  T3 <- -log(1 - fst3)
  PBS <- (T1 + T2 - T3) / 2
  return(PBS)
}

# Extract Fst values from the matrix
fst_EAS_EUR <- fst_matrix["EAS", "EUR"]
fst_EAS_AFR <- fst_matrix["EAS", "AFR"]
fst_EUR_AFR <- fst_matrix["EUR", "AFR"]

# Calculate PBS for EAS
PBS_EAS <- calculate_PBS(fst_EAS_EUR, fst_EAS_AFR, fst_EUR_AFR)
print(PBS_EAS)




