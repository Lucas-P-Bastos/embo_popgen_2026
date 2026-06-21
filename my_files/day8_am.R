### R Code Exercise: Pairwise $F_{ST}$ Calculation
#Write the R code necessary to perform the following:
#1. Read the pairwise $F_{ST}$ files from `input/`.
#2. Filter duplicate SNP positions and exclude NA values.
#3. Align the datasets by overlapping positions.
#4. Set negative $F_{ST}$ values to zero.
#5. Check Fst values at position `109513601`.
#6. Calculate distribution quantiles to determine if rs3827760 is an outlier.
#7. Plot pairwise Fst around `109513601` in a 10kb window, highlighting the candidate SNP.

setwd ("/Users/lucas/EMBO_PopGen_2026/Tabita_Hunemeier/Practical_Session_Selection_EMBO/input/Part_1_HumanDiversity")
library(ggplot2)
# Read the pairwise $F_{ST}$ files from `input/`.
afr_eas <- read.table("AFR_EAS.weir.fst", header = TRUE)
afr_eur <- read.table("AFR_EUR.weir.fst", header = TRUE)
eas_eur <- read.table("EAS_EUR.weir.fst", header = TRUE)
nam_eas <- read.table("Chr2_NAM_EAS.weir.fst", header = TRUE)
nam_eur <- read.table("Chr2_NAM_EUR.weir.fst", header = TRUE)

# Filter duplicate SNP positions and exclude NA values.
clean <- function(df) {
  df <- df[!duplicated(df$POS), ]
  df <- df[!is.na(df$WEIR_AND_COCKERHAM_FST), ]
  return(df)
}

afr_eas <- clean(afr_eas)
afr_eur <- clean(afr_eur)
eas_eur <- clean(eas_eur)
nam_eas <- clean(nam_eas)
nam_eur <- clean(nam_eur)

# Align the datasets by overlapping positions.
common_pos <- Reduce(intersect, list(
  afr_eas$POS, afr_eur$POS, eas_eur$POS
))

afr_eas <- afr_eas[afr_eas$POS %in% common_pos, ]
afr_eur <- afr_eur[afr_eur$POS %in% common_pos, ]
eas_eur <- eas_eur[eas_eur$POS %in% common_pos, ]

# Set negative $F_{ST}$ values to zero.
zero_neg <- function(df) {
  df$WEIR_AND_COCKERHAM_FST[df$WEIR_AND_COCKERHAM_FST < 0] <- 0
  return(df)
}

afr_eas <- zero_neg(afr_eas)
afr_eur <- zero_neg(afr_eur)
eas_eur <- zero_neg(eas_eur)
nam_eas <- zero_neg(nam_eas)
nam_eur <- zero_neg(nam_eur)

# Check Fst values at position `109513601`.
target_pos <- 109513601

cat("AFR-EAS FST at rs3827760:\n"); print(afr_eas[afr_eas$POS == target_pos, ])
cat("AFR-EUR FST at rs3827760:\n"); print(afr_eur[afr_eur$POS == target_pos, ])
cat("EAS-EUR FST at rs3827760:\n"); print(eas_eur[eas_eur$POS == target_pos, ])
cat("NAM-EAS FST at rs3827760:\n"); print(nam_eas[nam_eas$POS == target_pos, ])
cat("NAM-EUR FST at rs3827760:\n"); print(nam_eur[nam_eur$POS == target_pos, ])

# Quantiles to determine if rs3827760 is an outlier ────────────────────
probs <- c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999)

cat("\nEAS-EUR quantiles:\n")
print(quantile(eas_eur$WEIR_AND_COCKERHAM_FST, probs = probs))

cat("\nAFR-EAS quantiles:\n")
print(quantile(afr_eas$WEIR_AND_COCKERHAM_FST, probs = probs))

cat("\nAFR-EUR quantiles:\n")
print(quantile(afr_eur$WEIR_AND_COCKERHAM_FST, probs = probs))

# percentile rank of rs3827760 in each comparison
for (pair in list(
  list(df = eas_eur, name = "EAS-EUR"),
  list(df = afr_eas, name = "AFR-EAS"),
  list(df = afr_eur, name = "AFR-EUR")
)) {
  fst_val <- pair$df$WEIR_AND_COCKERHAM_FST[pair$df$POS == target_pos]
  pct     <- round(100 * mean(pair$df$WEIR_AND_COCKERHAM_FST <= fst_val), 2)
  cat(pair$name, "- rs3827760 FST:", fst_val, "| percentile:", pct, "\n")
}

# ── 7. Plot pairwise FST in 10kb window around rs3827760 ────────────────────
window_start <- target_pos - 5000
window_end   <- target_pos + 5000

plot_fst_window <- function(df, pair_name) {
  win <- df[df$POS >= window_start & df$POS <= window_end, ]
  win$is_target <- win$POS == target_pos
  
  ggplot(win, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
    geom_point(aes(color = is_target, size = is_target), alpha = 0.8) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red"),
                       labels = c("Other SNPs", "rs3827760 (EDAR)")) +
    scale_size_manual(values  = c("FALSE" = 1.5,     "TRUE" = 4)) +
    geom_vline(xintercept = target_pos, linetype = "dashed",
               color = "red", alpha = 0.5) +
    theme_minimal() +
    labs(
      title = paste("FST around rs3827760 (EDAR) —", pair_name),
      x     = "Position (chr2)",
      y     = expression(F[ST]),
      color = ""
    ) +
    theme(legend.position = "top")
}

plot_fst_window(eas_eur, "EAS vs EUR")
plot_fst_window(afr_eas, "AFR vs EAS")
plot_fst_window(afr_eur, "AFR vs EUR")

### R Code Exercise: Population Branch Statistics (PBS)
#Write the R code necessary to:
#1. Estimate the Population Branch Statistic for East Asians ($PBS_{EAS}$) using the AFR, EAS, and EUR populations.
# make sure all three datasets are aligned to the same positions
common_pos_pbs <- Reduce(intersect, list(
  afr_eas$POS, afr_eur$POS, eas_eur$POS
))

afr_eas_pbs <- afr_eas[afr_eas$POS %in% common_pos_pbs, ]
afr_eur_pbs <- afr_eur[afr_eur$POS %in% common_pos_pbs, ]
eas_eur_pbs <- eas_eur[eas_eur$POS %in% common_pos_pbs, ]

# FST values
fst_ae  <- afr_eas_pbs$WEIR_AND_COCKERHAM_FST  # AFR vs EAS
fst_aeu <- afr_eur_pbs$WEIR_AND_COCKERHAM_FST  # AFR vs EUR
fst_eu  <- eas_eur_pbs$WEIR_AND_COCKERHAM_FST  # EAS vs EUR

# branch lengths (T values)
# FST of 1 would give log(0) = -Inf so cap at 0.9999
fst_ae  <- pmin(fst_ae,  0.9999)
fst_aeu <- pmin(fst_aeu, 0.9999)
fst_eu  <- pmin(fst_eu,  0.9999)

T_ae  <- -log(1 - fst_ae)
T_aeu <- -log(1 - fst_aeu)
T_eu  <- -log(1 - fst_eu)

PBS_EAS <- (T_ae + T_eu - T_aeu) / 2

# build results dataframe
pbs_df <- data.frame(
  POS     = afr_eas_pbs$POS,
  CHROM   = afr_eas_pbs$CHROM,
  PBS_EAS = PBS_EAS
)
#2. Convert negative branch lengths to zero.
pbs_df$PBS_EAS[pbs_df$PBS_EAS < 0] <- 0
#3. Check the PBS value for the candidate SNP rs3827760.
target_pos <- 109513601
cat("PBS_EAS at rs3827760:\n")
print(pbs_df[pbs_df$POS == target_pos, ])

#4. Calculate distribution quantiles to determine if it is an outlier.
probs <- c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
cat("\nPBS_EAS quantiles:\n")
print(quantile(pbs_df$PBS_EAS, probs = probs))

pbs_val <- pbs_df$PBS_EAS[pbs_df$POS == target_pos]
pct     <- round(100 * mean(pbs_df$PBS_EAS <= pbs_val), 2)
cat("\nrs3827760 PBS_EAS:", pbs_val, "| percentile:", pct, "\n")
#5. Plot PBS values around the candidate SNP in a 10kb window.
window_start <- target_pos - 5000
window_end   <- target_pos + 5000

win <- pbs_df[pbs_df$POS >= window_start & pbs_df$POS <= window_end, ]
win$is_target <- win$POS == target_pos

ggplot(win, aes(x = POS, y = PBS_EAS)) +
  geom_point(aes(color = is_target, size = is_target), alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red"),
                     labels = c("Other SNPs", "rs3827760 (EDAR)")) +
  scale_size_manual(values  = c("FALSE" = 1.5,     "TRUE" = 4)) +
  geom_vline(xintercept = target_pos, linetype = "dashed",
             color = "red", alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "PBS (EAS) around rs3827760 (EDAR)",
    x     = "Position (chr2)",
    y     = "PBS (EAS)",
    color = ""
  ) +
  theme(legend.position = "top")

#Questions for Students

#1. **The estimate of $F_{ST}$ by the Weir and Cockerham metric can sometimes generate negative values and "NA". What does that mean? How can this interfere with the results?**
#  * *Answer: Missing data or low sample sizes can lead to unreliable estimates of genetic differentiation, resulting in negative $F_{ST}$ values or "NA" values. Negative values indicate that the observed genetic variation is less than expected under random mating, which can occur due to sampling error or other factors. These values can interfere with the results by skewing the distribution of $F_{ST}$ values and potentially leading to incorrect conclusions about population differentiation.*
#  2. **The $F_{ST}$ values observed between pairs of populations for the SNP rs3827760 (position 109,513,601) fall within which distribution quantiles of $F_{ST}$ values for the studied chromosome? Can they be considered outliers?**
#  * *Answer*: all negative values should be transformed in 0
#  3. **From the observed $F_{ST}$ values between population pairs and the significance estimates, what can we say about the rs3827760 SNP differentiation between populations?**
#  * *Answer*:Just in one of the cases. Between EAS and Eur. 
#  4. **Discuss how these results justify performing another type of analysis based on PBS (Population Branch Statistics).**
#  * *Answer*:
#  5. **What does the PBS analysis reveal? What is the difference between PBS and $F_{ST}$ analysis?**
#  * *Answer*: Direction. 


### Part 2: Genomic Selection Sweep Scan in Canines### 

# change wd to the canine file
setwd("/Users/lucas/EMBO_PopGen_2026/Tabita_Hunemeier/Practical_Session_Selection_EMBO/input/Part_2_CanidDiversity")

sample_info <- read.delim(
  "sample_info.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

head(sample_info)
names(sample_info)

# Read PLINK PCA output

# Rename PCA Collumns
colnames(eigenvec) <- c(
  "FID",
  "sampleName",
  paste0("PC", 1:(ncol(eigenvec)-2))
)

# Merge with metadata
pca <- merge(eigenvec, sample_info, by = "sampleName")


# Percentage of variance explained
var_exp <- eigenval/sum(eigenval)*100
pc1 <- round(var_exp[1],2)
pc2 <- round(var_exp[2],2)

# PCA plot:
ggplot(pca, aes(x = PC1, y = PC2, color = breed)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of Canine Genomic Data",
    x     = paste0("PC1 (", pc1, "%)"),
    y     = paste0("PC2 (", pc2, "%)"),
    color = "Breed"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# genomic outlier detection using pcadapt
install.packages("pcadapt")
library(pcadapt)
#Load genomic data in PLINK format
x <- read.pcadapt (~/subset_chr15.bed, type = "bed")

pcadapt(
  x,
  K = 2,
  method = "mahalanobis",
  min.maf = 0.05,
  ploidy = 2,
  LD.clumping = TRUE,
  pca.only = FALSE,
  tol = 1e-04
)
