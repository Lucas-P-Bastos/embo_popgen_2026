install.packages("rehh")
library(rehh)

setwd("/Users/lucas/EMBO_PopGen_2026/Tabita_Hunemeier/Practical_Session_Selection_EMBO/input/Part_1_HumanDiversity")

## EHH & Furcation Trees ##

# code necessary to: 
#1. Convert the VCF databases to rehh format using data2haplohh(). 
#2. Estimate the EHH decay for rs3827760 in both populations. 
#3. Plot the EHH decay and furcation trees for both AFR and EAS.

# Load the VCF files for both populations
vsfile1 <- "Chr2_EDAR_CHS_500K.recode.vcf"
vsfile2 <- "Chr2_EDAR_LWK_500K.recode.vcf"

# Convert VCF to rehh format
haplohh1 <- data2haplohh(vsfile1, polarize_vcf = FALSE)
haplohh2 <- data2haplohh(vsfile2, polarize_vcf = FALSE)

# Estimate EHH decay for rs3827760 in both populations
ehh_calc_EAS <- calc_ehh(haplohh1, mrk = "rs3827760")
ehh_calc_AFR <- calc_ehh(haplohh2, mrk = "rs3827760")

# Plot EHH decay for both populations

plot(ehh_calc_EAS, main = "EHH decay - EAS (CHS) at rs3827760")
plot(ehh_calc_AFR, main = "EHH decay - AFR (LWK) at rs3827760")

# Furcation trees
frc_EAS <- calc_furcation(haplohh1, mrk = "rs3827760")
frc_AFR <- calc_furcation(haplohh2, mrk = "rs3827760")

plot(frc_EAS, main = "Furcation tree - EAS at rs3827760")
plot(frc_AFR, main = "Furcation tree - AFR at rs3827760")

## iHS & XP-EHH (Window-based) ##
#Write the R code necessary to: 
#1. Perform a genome-wide scan of homozygosity using scan_hh() for AFR and EAS. 
#2. Calculate iHS scores for both populations using ihh2ihs(). 
#3. Check the iHS score at rs3827760 and generate a single-site iHS plot in EAS. 
#4. Create a function to estimate the average absolute iHS in sliding windows (50 SNPs/40 step) and plot the results. 
#5. Estimate cross-population XP-EHH between EAS and AFR using ies2xpehh(), calculate window-based averages, and plot them.


# Perform genome-wide scan of homozygosity
hh_scan_EAS <- scan_hh(haplohh1)
hh_scan_AFR <- scan_hh(haplohh2)

# Calculate iHS scores for both populations using ihh2ihs()
ihs_EAS <- ihh2ihs(hh_scan_EAS, min_maf = 0.02, freqbin = 0.01)
ihs_AFR <- ihh2ihs(hh_scan_AFR, min_maf = 0.02, freqbin = 0.01)

# extract the iHS scores dataframe from the list
ihs_EAS_df <- ihs_EAS$ihs
ihs_AFR_df <- ihs_AFR$ihs

# check the structure and column names
head(ihs_EAS_df)
colnames(ihs_EAS_df)

# basic iHS plot for EAS highlighting rs3827760
target_pos <- 109513601

plot(ihs_EAS_df$POSITION, ihs_EAS_df$IHS,
     col = ifelse(ihs_EAS_df$POSITION == target_pos, "red", "grey50"),
     pch = 19, cex = 0.5,
     main = "iHS scores - EAS (CHS)",
     xlab = "Position (chr2)",
     ylab = "iHS")

# highlight the candidate SNP
abline(v = target_pos, col = "red", lty = 2)

# add a legend
legend("topright", legend = c("rs3827760 (EDAR)", "Other SNPs"),
       col = c("red", "grey50"), pch = 19)

# create a function to estimate the average absolute iHS in sliding windows (50 SNPs/40 step)
# to be continued.....

###  PBS in Native Americans (NAM)
#1. Read the pairwise FST files involving Native Americans (input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst and input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst) and Europeans-East Asians (input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst). 
#2. Filter duplicates, exclude NA values, and align positions. 
#3. Convert negative FST values to zero. 
#4. Estimate PBSNAM using NAM, EAS, and EUR. 5. Check PBS value at rs3827760, check quantiles, and plot the PBS scan.


#read the pairwise FST files
fst_NAM_EAS <- read.table("Chr2_NAM_EAS.weir.fst", header = TRUE)
fst_NAM_EUR <- read.table("Chr2_NAM_EUR.weir.fst", header = TRUE)
fst_EUR_EAS <- read.table("Chr2_EUR_EAS.weir.fst", header = TRUE)

# filter duplicates, exclude NA values, and align positions
clean <- function(df) {
  df <- df[!duplicated(df$POS), ]
  df <- df[!is.na(df$WEIR_AND_COCKERHAM_FST), ]
  return(df)
}
fst_NAM_EAS<- clean(fst_NAM_EAS)
fst_NAM_EUR<- clean(fst_NAM_EUR)
fst_EUR_EAS<- clean(fst_EUR_EAS)

# align positions
common_pos <- Reduce(intersect, list(
  afr_eas$POS, afr_eur$POS, eas_eur$POS
))

fst_NAM_EAS <- fst_NAM_EAS[fst_NAM_EAS$POS %in% common_pos, ]
fst_NAM_EUR <- fst_NAM_EUR[fst_NAM_EUR$POS %in% common_pos, ]
fst_EUR_EAS <- fst_EUR_EAS[fst_EUR_EAS$POS %in% common_pos, ]

# convert negative FST values to zero

# Set negative $F_{ST}$ values to zero.
zero_neg <- function(df) {
  df$WEIR_AND_COCKERHAM_FST[df$WEIR_AND_COCKERHAM_FST < 0] <- 0
  return(df)
}

fst_NAM_EAS <- zero_neg(fst_NAM_EAS)
fst_NAM_EUR <- zero_neg(fst_NAM_EUR)
fst_EUR_EAS <- zero_neg(fst_EUR_EAS)


# Estimate PBSNAM using NAM, EAS, and EUR
plot(fst_NAM_EAS$POS, fst_NAM_EAS$WEIR_AND_COCKERHAM_FST,
     type = "l", col = "blue", lwd = 2,
     main = "Pairwise FST - NAM vs EAS",
     xlab = "Position (chr2)",
     ylab = "FST")  # <-- closing parenthesis was missing
     