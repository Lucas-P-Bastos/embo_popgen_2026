# EMBO Practical Course: Genomic Diversity & Natural Selection Scan

This tutorial guides students through the analysis of genomic data in humans and canines to identify genomic regions under natural selection. The tutorial is divided into two main parts:

1. **Part 1: Human Genomic Diversity and Natural Selection**: Investigating selection signatures in the candidate gene ***EDAR*** (associated with ectodermal traits in East Asians and Native Americans) using population differentiation ($F_{ST}$, PBS) and haplotype-based metrics (EHH, iHS, XP-EHH).
2. **Part 2: Genomic Selection Scan in Canines**: Identifying the selective sweep at the ***IGF1*** body-size locus by comparing small vs. large dog breeds using PCA, PCAdapt (outlier scan), and haplotype homozygosity methods (XP-nSL and Rsb).

---

# Part 1: Human Genomic Diversity and Natural Selection

## 1. Background and Dataset

### Goal
Our goal is to explore approaches and methods which seek to identify regions of the genome with signatures of natural selection. We will use real genomic data and two classes of tests: one based on population differentiation ($F_{ST}$ / PBS) and another based on extended haplotype homozygosity (EHH / iHS / XP-EHH).

### Dataset
Whole-genome sequencing data from the 1000 Genomes Project Phase III. The full database can be accessed via:
<ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>

### Data Pre-processing
To optimize class time, we will analyze a pre-processed dataset for chromosome 2 corresponding to individuals sampled from the African (AFR: 504 individuals), European (EUR: 503 individuals), and East Asian (EAS: 504 individuals) populations. In this dataset, INDELs, singletons, and SNPs with MAF < 0.05 have been removed. The pairwise $F_{ST}$ was then estimated using `vcftools`.

All data files are located in the `input/` directory:
- `input/Part_1_HumanDiversity/AFR_EAS.weir.fst` (Fst between Africans and East Asians)
- `input/Part_1_HumanDiversity/AFR_EUR.weir.fst` (Fst between Africans and Europeans)
- `input/Part_1_HumanDiversity/EAS_EUR.weir.fst` (Fst between East Asians and Europeans)
- `input/Part_1_HumanDiversity/Chr2_EDAR_LWK_500K.recode.vcf` (Phased African haplotypes around *EDAR*)
- `input/Part_1_HumanDiversity/Chr2_EDAR_CHS_500K.recode.vcf` (Phased East Asian haplotypes around *EDAR*)
- `input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst` (Fst between Native Americans and East Asians in candidate region)
- `input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst` (Fst between Native Americans and Europeans in candidate region)
- `input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst` (Fst between Europeans and East Asians in candidate region)

---

## 2. Genetic Differentiation ($F_{ST}$ and PBS)

### Investigating the Candidate Gene *EDAR*
The human Ectodysplasin A receptor gene, or ***EDAR***, is part of the EDA signaling pathway which specifies prenatally the location, size, and shape of ectodermal appendages (such as hair follicles, teeth, and glands). *EDAR* is a textbook example of positive selection in East Asians, with genomic and functional experiments corroborating it. A specific non-synonymous variant, **rs3827760** (chr2:109,513,601 A>G), results in a Val370Ala substitution and is strongly associated with thicker hair shafts and shovel-shaped incisors. Another hypothesis states that *EDAR* acted along with *FADS* and *VDR* in the Beringia Standstill, allowing Native American ancestors to survive in extreme arctic environments.

### Questions for Students

1. **The estimate of $F_{ST}$ by the Weir and Cockerham metric can sometimes generate negative values and "NA". What does that mean? How can this interfere with the results?**
   * *Answer*: Weir and Cockerham's $F_{ST}$ is an unbiased estimator that accounts for sample size. Negative values occur when the observed variance within populations is greater than the variance between populations (typical for neutral variants with little to no differentiation). Biologically, negative $F_{ST}$ is meaningless and is treated as 0. "NA" values occur when a SNP is monomorphic (zero variance) or has entirely missing data across the compared populations, meaning the denominator of the $F_{ST}$ equation becomes zero. This can interfere with analyses if not filtered, as R will produce errors or skew statistical distributions.
2. **The $F_{ST}$ values observed between pairs of populations for the SNP rs3827760 (position 109,513,601) fall within which distribution quantiles of $F_{ST}$ values for the studied chromosome? Can they be considered outliers?**
   * *Answer*:
     * **AFR vs. EAS**: $F_{ST} = 0.8729$, which falls in the **99.9%** quantile (the 99% threshold is 0.6521). It is a highly significant outlier.
     * **EAS vs. EUR**: $F_{ST} = 0.8591$, which falls in the **99.9%** quantile (the 99% threshold is 0.4790). It is a highly significant outlier.
     * **AFR vs. EUR**: $F_{ST} = 0.0099$, which is below the median ($F_{ST} \approx 0.081$). It is not an outlier.
3. **From the observed $F_{ST}$ values between population pairs and the significance estimates, what can we say about the rs3827760 SNP differentiation between populations?**
   * *Answer*: The SNP rs3827760 is highly differentiated between East Asians (EAS) and both Africans (AFR) and Europeans (EUR), but has almost zero differentiation between Africans and Europeans. This indicates that the allele frequency changed dramatically specifically in the lineage leading to East Asians after they split from Europeans and Africans.
4. **Discuss how these results justify performing another type of analysis based on PBS (Population Branch Statistics).**
   * *Answer*: Pairwise $F_{ST}$ can show that East Asians are highly different from both Africans and Europeans. However, pairwise metrics alone cannot pinpoint which lineage underwent the change—it could be that East Asians were selected, or both Africans and Europeans changed independently in the same direction. By using PBS, which incorporates a third population (outgroup), we can estimate lineage-specific branch lengths and formally confirm that the evolutionary change occurred specifically along the East Asian branch.
5. **What does the PBS analysis reveal? What is the difference between PBS and $F_{ST}$ analysis?**
   * *Answer*: The PBS analysis reveals an exceptionally long branch for East Asians at the candidate SNP ($PBS_{EAS} = 2.006$), which is in the extreme **99.9%** quantile (99th percentile is 0.6078). This confirms that a massive, population-specific allele frequency shift took place along the East Asian lineage. The difference is that $F_{ST}$ measures pairwise genetic distance between two populations, whereas PBS projects these pairwise distances onto a three-population tree to isolate lineage-specific change.

---

### R Code Exercise: Pairwise $F_{ST}$ Calculation

```R
# 1. Read files with Fst estimates
names_header <- c("CHROM", "POS", "WEIR_AND_COCKERHAM_FST", "NUM", "DEN")

FST_AFR_EAS <- read.table("input/Part_1_HumanDiversity/AFR_EAS.weir.fst", header = FALSE, skip = 1, col.names = names_header, fill = TRUE)
FST_AFR_EUR <- read.table("input/Part_1_HumanDiversity/AFR_EUR.weir.fst", header = FALSE, skip = 1, col.names = names_header, fill = TRUE)
FST_EAS_EUR <- read.table("input/Part_1_HumanDiversity/EAS_EUR.weir.fst", header = FALSE, skip = 1, col.names = names_header, fill = TRUE)

# 2. Remove duplicates
FST_AFR_EAS_filter <- FST_AFR_EAS[!duplicated(FST_AFR_EAS$POS), ]
FST_AFR_EUR_filter <- FST_AFR_EUR[!duplicated(FST_AFR_EUR$POS), ]
FST_EAS_EUR_filter <- FST_EAS_EUR[!duplicated(FST_EAS_EUR$POS), ]

# 3. Exclude NAs
FST_AfrEas_data <- FST_AFR_EAS_filter[!is.na(FST_AFR_EAS_filter[, 3]), ]
FST_AfrEur_data <- FST_AFR_EUR_filter[!is.na(FST_AFR_EUR_filter[, 3]), ]
FST_EasEur_data <- FST_EAS_EUR_filter[!is.na(FST_EAS_EUR_filter[, 3]), ]

# 4. Overlap SNPs across datasets
overlap_AfrEas_AfrEur <- FST_AfrEas_data[FST_AfrEas_data$POS %in% FST_AfrEur_data$POS, ]
overlap_AfrEasEur_EasEur <- overlap_AfrEas_AfrEur[overlap_AfrEas_AfrEur$POS %in% FST_EasEur_data$POS, ]

FST_AfrEas_data_clean <- FST_AfrEas_data[FST_AfrEas_data$POS %in% overlap_AfrEasEur_EasEur$POS, ]
FST_AfrEur_data_clean <- FST_AfrEur_data[FST_AfrEur_data$POS %in% overlap_AfrEasEur_EasEur$POS, ]
FST_EasEur_data_clean <- FST_EasEur_data[FST_EasEur_data$POS %in% overlap_AfrEasEur_EasEur$POS, ]

# Sort by position
FST_AfrEas_data_clean <- FST_AfrEas_data_clean[order(FST_AfrEas_data_clean$POS), ]
FST_AfrEur_data_clean <- FST_AfrEur_data_clean[order(FST_AfrEur_data_clean$POS), ]
FST_EasEur_data_clean <- FST_EasEur_data_clean[order(FST_EasEur_data_clean$POS), ]

# 5. Convert negative values to zero
FST_AfrEas_data_clean[FST_AfrEas_data_clean[, 3] < 0, 3] <- 0
FST_AfrEur_data_clean[FST_AfrEur_data_clean[, 3] < 0, 3] <- 0
FST_EasEur_data_clean[FST_EasEur_data_clean[, 3] < 0, 3] <- 0

# Candidate SNP position
POS <- 109513601

# Calculate quantiles
FST_AfrEas_distrQT <- quantile(FST_AfrEas_data_clean[, 3], c(0.01, 0.05, 0.1, .25, .50, .75, .90, 0.95, .99))
FST_AfrEur_distrQT <- quantile(FST_AfrEur_data_clean[, 3], c(0.01, 0.05, 0.1, .25, .50, .75, .90, 0.95, .99))
FST_EasEur_distrQT <- quantile(FST_EasEur_data_clean[, 3], c(0.01, 0.05, 0.1, .25, .50, .75, .90, 0.95, .99))

# Define adjacent plotting window (10,000 bp around candidate SNP)
SNPfrom_BP <- POS - 10000
SNPto_BP <- POS + 10000

SNPfrom_id_AfrEas <- max(which(FST_AfrEas_data_clean[, 2] <= SNPfrom_BP))
SNPto_id_Afr_Eas <- min(which(FST_AfrEas_data_clean[, 2] >= SNPto_BP))
FSTdata_SNP_AfrEas <- FST_AfrEas_data_clean[SNPfrom_id_AfrEas:SNPto_id_Afr_Eas, ]

SNPfrom_id_AfrEur <- max(which(FST_AfrEur_data_clean[, 2] <= SNPfrom_BP))
SNPto_id_Afr_Eur <- min(which(FST_AfrEur_data_clean[, 2] >= SNPto_BP))
FSTdata_SNP_AfrEur <- FST_AfrEur_data_clean[SNPfrom_id_AfrEur:SNPto_id_Afr_Eur, ]

SNPfrom_id_EasEur <- max(which(FST_EasEur_data_clean[, 2] <= SNPfrom_BP))
SNPto_id_EasEur <- min(which(FST_EasEur_data_clean[, 2] >= SNPto_BP))
FSTdata_SNP_EasEur <- FST_EasEur_data_clean[SNPfrom_id_EasEur:SNPto_id_EasEur, ]

# Plot pairwise Fst around rs3827760 (outliers marked with dashed red lines representing the 95th percentile)
plot(ylim = c(0, 1), x = FSTdata_SNP_AfrEas[, 2], y = FSTdata_SNP_AfrEas[, 3], xlab = 'pos', ylab = 'FST AFR EAS', pch = 20, cex = 1.5)
points(x = FSTdata_SNP_AfrEas[FSTdata_SNP_AfrEas[, 2] == POS, 2], y = FSTdata_SNP_AfrEas[FSTdata_SNP_AfrEas[, 2] == POS, 3], col = 'blue', cex = 2, lwd = 2)
abline(h = FST_AfrEas_distrQT[[8]], lty = 2, col = "red")
```

### Resulting Figures: Pairwise $F_{ST}$ around rs3827760

#### $F_{ST}$ AFR vs. EAS
![Fst AFR EAS](figures/figs_human_diversity/fst_afr_eas.png)

#### $F_{ST}$ AFR vs. EUR
![Fst AFR EUR](figures/figs_human_diversity/fst_afr_eur.png)

#### $F_{ST}$ EAS vs. EUR
![Fst EAS EUR](figures/figs_human_diversity/fst_eas_eur.png)

---

### R Code Exercise: Population Branch Statistics (PBS)

```R
# 1. Perform PBS test, using EAS as candidate population
# PBS = (T_AB + T_AC - T_BC) / 2 where T = -log(1 - Fst)
PBS_EAS <- ((-log(1 - FST_AfrEas_data_clean$WEIR_AND_COCKERHAM_FST)) + 
            (-log(1 - FST_EasEur_data_clean$WEIR_AND_COCKERHAM_FST)) - 
            (-log(1 - FST_AfrEur_data_clean$WEIR_AND_COCKERHAM_FST))) / 2

# Convert negative branch lengths to zero
PBS_EAS[PBS_EAS < 0] <- 0

# Combine into data frame
fst_pbs <- as.data.frame(cbind(FST_EasEur_data_clean$POS, 
                               FST_AfrEas_data_clean$WEIR_AND_COCKERHAM_FST, 
                               FST_AfrEur_data_clean$WEIR_AND_COCKERHAM_FST, 
                               FST_EasEur_data_clean$WEIR_AND_COCKERHAM_FST, 
                               PBS_EAS), stringsAsFactors = FALSE)
names(fst_pbs) <- c("POS", "FST_AFR_EAS", "FST_AFR_EUR", "FST_EAS_EUR", "PBS_EAS")

# Quantile estimation
PBS_distrQT <- quantile(PBS_EAS, c(0.01, 0.05, 0.1, .25, .50, .75, .90, 0.95, .99))

# Define adjacent plotting window (10,000 bp around candidate SNP)
SNP_FROM <- POS - 10000
SNP_TO <- POS + 10000

# Subset region for plotting
SNPfrom_PBS <- max(which(fst_pbs$POS <= SNP_FROM))
SNPto_PBS <- min(which(fst_pbs$POS >= SNP_TO))
subset_fst_PBS <- fst_pbs[SNPfrom_PBS:SNPto_PBS, ]

# Plot PBS
plot(ylim = c(0, 2.5), x = subset_fst_PBS$POS, y = subset_fst_PBS$PBS_EAS, xlab = 'pos', ylab = 'PBS', pch = 20, cex = 1.5)
points(x = subset_fst_PBS[subset_fst_PBS$POS == POS, 1], y = subset_fst_PBS[subset_fst_PBS$POS == POS, 5], col = 'blue', cex = 2, lwd = 2)
abline(h = PBS_distrQT[[9]], lty = 2, col = "red") # 99% threshold
```

### Resulting Figure: PBS scan for EAS
![PBS EAS](figures/figs_human_diversity/pbs_eas.png)

---

## 3. Extended Haplotype Homozygosity (EHH)

### Extended Haplotype Homozygosity (EHH) and Haplotype Sweeps
Different approaches can detect genomic signatures of selection at different timescales. More recent selection signals can be detected from haplotype-based tests. Positive selection causes a rapid rise in the frequency of the selected allele, such that recombination does not have enough time to break down the haplotype on which the mutation arose. This creates a signature of **Extended Haplotype Homozygosity (EHH)** extending over a long physical distance.

### Questions for Students

1. **How is the haplotype profile of genetic variants under recent positive selection?**
   * *Answer*: A variant under recent positive selection shows a very high frequency accompanied by high haplotype homozygosity (long haplotype blocks) that extends over long physical distances. In neutral conditions, high frequency alleles are old, meaning recombination has broken down their haplotype background (short haplotypes). Recent positive selection accelerates an allele to high frequency so rapidly that recombination has not had time to break the haplotype association, resulting in extremely long conserved haplotypes.
2. **What is the profile of ancestral and derived haplotypes of the rs3827760 SNP in AFR and EAS?**
   * *Answer*:
     * **AFR (Africans)**: The derived allele is absent (`FREQ_D = 0`), and only the ancestral allele exists (`FREQ_A = 1`).
     * **EAS (East Asians)**: The derived allele is at high frequency (`FREQ_D = 0.9048`), and the EHH decay plot shows that haplotype homozygosity for the derived allele (`EHH_D`) extends over hundreds of kilobases with very slow decay (`IHH_D = 55,231.78` vs. `IHH_A = 9,979.00`). This reveals a massive selective sweep on the derived allele in East Asians.
3. **The iHS score observed for the SNP rs3827760 falls within which distribution quantiles of iHS values for the studied chromosome? Can it be considered an outlier? How can we make this analysis more robust?**
   * *Answer*: The single-site iHS score for rs3827760 in EAS is `-1.588` (log p-value = 0.95), which is not an extreme outlier (standard threshold is $|iHS| > 2$). This occurs because the derived allele has reached near-fixation (~90% frequency), causing the ancestral allele to be present in very few copies, which reduces the statistical power of the single-site test. We can make this analysis more robust by using a window-based approach, calculating the mean of absolute iHS values across sliding windows of SNPs. The window containing rs3827760 has a mean iHS of `1.146`, which is in the extreme **95th-99th** percentile, revealing the sweep clearly.
4. **What information does the XP-EHH analysis add about natural selection in the candidate SNP?**
   * *Answer*: XP-EHH compares the integrated EHH between two populations. It is designed to detect selective sweeps where the selected allele has approached or achieved fixation in one population but remains polymorphic in the other. While iHS loses power near fixation, XP-EHH retains strong power. The XP-EHH window scan between EAS and AFR shows a highly significant selection signal at the window containing rs3827760 (`mean_xpEHH = 1.956`, which is in the **>99th** percentile), providing strong confirmation of selection in East Asians.

---

### R Code Exercise: EHH & Furcation Trees

```R
library("rehh")

# 1. Load VCF files (data is unpolarized)
data1 <- data2haplohh(hap_file = "input/Part_1_HumanDiversity/Chr2_EDAR_LWK_500K.recode.vcf", polarize_vcf = FALSE, vcf_reader = "data.table")
data2 <- data2haplohh(hap_file = "input/Part_1_HumanDiversity/Chr2_EDAR_CHS_500K.recode.vcf", polarize_vcf = FALSE, vcf_reader = "data.table")

# 2. Calculate EHH for rs3827760
ehh_calc_AFR <- calc_ehh(data1, mrk = "rs3827760")
ehh_calc_EAS <- calc_ehh(data2, mrk = "rs3827760")

# Plot EHH
plot(ehh_calc_AFR)
plot(ehh_calc_EAS)

# 3. Calculate furcation trees
furcation_afr <- calc_furcation(data1, mrk = "rs3827760")
plot(furcation_afr)

furcation_eas <- calc_furcation(data2, mrk = "rs3827760")
plot(furcation_eas)
```

### Resulting Figures: EHH decay & Furcation Trees

#### EHH Decay in AFR
![EHH AFR](figures/figs_human_diversity/ehh_afr.png)

#### EHH Decay in EAS
![EHH EAS](figures/figs_human_diversity/ehh_eas.png)

#### Haplotype Furcation Tree in AFR
![Furcation AFR](figures/figs_human_diversity/furcation_afr.png)

#### Haplotype Furcation Tree in EAS
![Furcation EAS](figures/figs_human_diversity/furcation_eas.png)

---

### R Code Exercise: iHS & XP-EHH (Window-based)

```R
# 1. Scan haplotype homozygosity across chromosome
AFR <- scan_hh(data1)
EAS <- scan_hh(data2)

# 2. Estimate iHS
iHS.AFR <- ihh2ihs(AFR, min_maf = 0.02, freqbin = 0.01)
iHS.EAS <- ihh2ihs(EAS, min_maf = 0.02, freqbin = 0.01)

# Plot single site iHS EAS
plot(iHS.EAS$ihs$POSITION, iHS.EAS$ihs$IHS, col = ifelse(iHS.EAS$ihs$POSITION == 109513601, "red", "black"), pch = 19)
abline(h = c(2, -2), lty = 2)

# 3. Sliding window analysis functions
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(abs(data[spots[i]:(spots[i] + window - 1)]), na.rm = TRUE)
  }
  return(result)
}

slidePos <- function(data, window, step){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- data[spots[i]]
  }
  return(result)
}

# 4. Average iHS over a window of 50 SNPs with steps of 40 SNPs
mean_iHS <- slideFunct(iHS.EAS$ihs$IHS, 50, 40)
pos_wind_Eas <- slidePos(iHS.EAS$ihs$POSITION, 50, 40)
wind_iHS <- as.data.frame(cbind(pos_wind_Eas, mean_iHS), stringsAsFactors = FALSE)

# Identify the window that contains the candidate variant rs3827760 (position 109513601)
Row_WIND_iHS <- wind_iHS[wind_iHS$pos_wind_Eas <= 109513601, ]
POS_WIND_iHS <- max(Row_WIND_iHS$pos_wind_Eas)

windiHS_distrQT <- quantile(wind_iHS$mean_iHS, c(0.01, 0.05, 0.1, .25, .50, .75, .90, 0.95, .99), na.rm = TRUE)

# Plot mean window iHS
plot(ylim = c(0, 1.5), x = wind_iHS[, 1], y = wind_iHS[, 2], xlab = 'pos', ylab = 'iHS windows', pch = 20, cex = 1.5)
points(x = wind_iHS[which(wind_iHS[, 1] == POS_WIND_iHS), 1], y = wind_iHS[which(wind_iHS[, 1] == POS_WIND_iHS), 2], col = 'red', cex = 2)
abline(h = windiHS_distrQT[[7]], lty = 2)

# 5. Calculate cross-population XP-EHH between EAS and AFR
xpEHH.EAS.AFR <- ies2xpehh(EAS, AFR)
mean_xpEHH <- slideFunct(xpEHH.EAS.AFR$XPEHH, 50, 40)
pos_wind_Eas_xp <- slidePos(xpEHH.EAS.AFR$POSITION, 50, 40)
wind_xpEHH <- as.data.frame(cbind(pos_wind_Eas_xp, mean_xpEHH), stringsAsFactors = FALSE)

# Identify the window that contains the candidate variant rs3827760 (position 109513601)
Row_WIND_xpEHH <- wind_xpEHH[wind_xpEHH$pos_wind_Eas_xp <= 109513601, ]
POS_WIND_xpEHH <- max(Row_WIND_xpEHH$pos_wind_Eas_xp)

windxpEHH_distrQT <- quantile(wind_xpEHH$mean_xpEHH, c(0.01, 0.05, 0.1, .25, .50, .75, .90, 0.95, .99), na.rm = TRUE)

# Plot windowed XP-EHH
plot(ylim = c(0, 2.05), x = wind_xpEHH[, 1], y = wind_xpEHH[, 2], xlab = 'pos', ylab = 'xpEHH windows', pch = 20, cex = 1.5)
points(x = wind_xpEHH[which(wind_xpEHH[, 1] == POS_WIND_xpEHH), 1], y = wind_xpEHH[which(wind_xpEHH[, 1] == POS_WIND_xpEHH), 2], col = 'red', cex = 2)
abline(h = windxpEHH_distrQT[[8]], lty = 2)
abline(v = c(109500000, 109605000), col = "red", lty = 2)
```

### Resulting Figures: iHS & XP-EHH Manhattan Plots

#### Single-site iHS in EAS
![iHS EAS](figures/figs_human_diversity/ihs_eas.png)

#### Window-based iHS in EAS
![iHS window EAS](figures/figs_human_diversity/ihs_window_eas.png)

#### Window-based XP-EHH EAS vs. AFR
![XP-EHH window](figures/figs_human_diversity/xpehh_window.png)

---

## 4. Native American Selection Analysis

### Background
Hlusko et al. (2018), using morphological data, found a strong selection signal in the *EDAR* gene in Native Americans. Using the additional database from the 1000 Genomes Project (Peruvian samples with over 95% Native American Ancestry, represented as **NAM**), we evaluate genomic signatures of selection at the functional variant rs3827760.

### Questions for Students

1. **Is the functional allele in East Asian at high frequency in other human populations (e.g. Native Americans)?**
   * *Answer*: Yes, the functional derived allele rs3827760-G is at very high frequency in Native Americans. This is shown by the $F_{ST}$ between Native Americans (NAM) and East Asians (EAS) being practically zero ($F_{ST} = -0.0107$, treated as 0), and the $F_{ST}$ between Native Americans and Europeans (EUR) being extremely high ($F_{ST} = 0.9661$), indicating that both NAM and EAS are almost fixed for the derived allele compared to Europeans.
2. **Can we identify signatures of natural selection on EDAR in Native Americans using PBS?**
   * *Answer*: Yes. By computing the Population Branch Statistic for Native Americans ($PBS_{NAM}$) using the NAM-EAS-EUR triplet, we find a high branch length ($PBS_{NAM} = 0.712$) at the functional variant. This value falls in the top **>99%** quantile of the chromosome 2 candidate region distribution, confirming a strong lineage-specific selection signature.
3. **Is selection targeting the same functional variant?**
   * *Answer*: Yes, the selection signature peaks at the exact same functional variant rs3827760 (position 109,513,601), showing that both East Asian and Native American selections targeted the Val370Ala mutation.
4. **What is your conclusion based on the data generated?**
   * *Answer*: The functional variant rs3827760 in the *EDAR* gene underwent positive selection in the common ancestral population of East Asians and Native Americans (e.g., during the Beringian Standstill) before the split, or underwent parallel selection. This supports the Hlusko et al. (2018) hypothesis, suggesting that this ectodermal variant was adaptive in arctic or sub-arctic environments (possibly related to vitamin D metabolism, glandular function, or UV radiation adaptation).

---

### R Code Exercise: PBS in Native Americans (NAM)

```R
# 1. Read files
FST_NAM_EAS <- read.table("input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst", header = FALSE, skip = 1, col.names = names_header, fill = TRUE)
FST_NAM_EUR <- read.table("input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst", header = FALSE, skip = 1, col.names = names_header, fill = TRUE)
FST_EUR_EAS_reg <- read.table("input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst", header = FALSE, skip = 1, col.names = names_header, fill = TRUE)

# 2. Remove duplicates
FST_NAM_EAS_filter <- FST_NAM_EAS[!duplicated(FST_NAM_EAS$POS), ]
FST_NAM_EUR_filter <- FST_NAM_EUR[!duplicated(FST_NAM_EUR$POS), ]
FST_EUR_EAS_reg_filter <- FST_EUR_EAS_reg[!duplicated(FST_EUR_EAS_reg$POS), ]

# 3. Exclude NAs
FST_NAM_EAS_data <- FST_NAM_EAS_filter[!is.na(FST_NAM_EAS_filter[, 3]), ]
FST_NAM_EUR_data <- FST_NAM_EUR_filter[!is.na(FST_NAM_EUR_filter[, 3]), ]
FST_EUR_EAS_reg_data <- FST_EUR_EAS_reg_filter[!is.na(FST_EUR_EAS_reg_filter[, 3]), ]

# 4. Overlap SNPs
overlap_NAM <- FST_NAM_EAS_data[FST_NAM_EAS_data$POS %in% FST_NAM_EUR_data$POS, ]
overlap_NAM_all <- overlap_NAM[overlap_NAM$POS %in% FST_EUR_EAS_reg_data$POS, ]

FST_NAM_EAS_clean <- FST_NAM_EAS_data[FST_NAM_EAS_data$POS %in% overlap_NAM_all$POS, ]
FST_NAM_EUR_clean <- FST_NAM_EUR_data[FST_NAM_EUR_data$POS %in% overlap_NAM_all$POS, ]
FST_EUR_EAS_reg_clean <- FST_EUR_EAS_reg_data[FST_EUR_EAS_reg_data$POS %in% overlap_NAM_all$POS, ]

# Sort by position
FST_NAM_EAS_clean <- FST_NAM_EAS_clean[order(FST_NAM_EAS_clean$POS), ]
FST_NAM_EUR_clean <- FST_NAM_EUR_clean[order(FST_NAM_EUR_clean$POS), ]
FST_EUR_EAS_reg_clean <- FST_EUR_EAS_reg_clean[order(FST_EUR_EAS_reg_clean$POS), ]

# Convert negative values to zero
FST_NAM_EAS_clean[FST_NAM_EAS_clean[, 3] < 0, 3] <- 0
FST_NAM_EUR_clean[FST_NAM_EUR_clean[, 3] < 0, 3] <- 0
FST_EUR_EAS_reg_clean[FST_EUR_EAS_reg_clean[, 3] < 0, 3] <- 0

# 5. Calculate PBS for NAM
# Topologies: NAM (A), EAS (B), EUR (C)
# PBS_NAM = (T_NAM_EAS + T_NAM_EUR - T_EAS_EUR) / 2
PBS_NAM <- ((-log(1 - FST_NAM_EAS_clean$WEIR_AND_COCKERHAM_FST)) + 
            (-log(1 - FST_NAM_EUR_clean$WEIR_AND_COCKERHAM_FST)) - 
            (-log(1 - FST_EUR_EAS_reg_clean$WEIR_AND_COCKERHAM_FST))) / 2
PBS_NAM[PBS_NAM < 0] <- 0

fst_pbs_nam <- data.frame(
  POS = FST_NAM_EAS_clean$POS,
  FST_NAM_EAS = FST_NAM_EAS_clean$WEIR_AND_COCKERHAM_FST,
  FST_NAM_EUR = FST_NAM_EUR_clean$WEIR_AND_COCKERHAM_FST,
  FST_EAS_EUR = FST_EUR_EAS_reg_clean$WEIR_AND_COCKERHAM_FST,
  PBS_NAM = PBS_NAM
)

PBS_NAM_distrQT <- quantile(PBS_NAM, c(0.01, 0.05, 0.1, .25, .50, .75, .90, 0.95, .99))

# Plot PBS NAM around rs3827760
plot(ylim = c(0, 2.5), x = fst_pbs_nam$POS, y = fst_pbs_nam$PBS_NAM, xlab = 'pos', ylab = 'PBS NAM', pch = 20, cex = 1.5)
points(x = fst_pbs_nam[fst_pbs_nam$POS == POS, 1], y = fst_pbs_nam[fst_pbs_nam$POS == POS, 5], col = 'blue', cex = 2, lwd = 2)
abline(h = PBS_NAM_distrQT[[9]], lty = 2, col = "red")
```

### Resulting Figure: PBS scan for NAM
![PBS NAM](figures/figs_human_diversity/pbs_nam.png)

---
---

# Part 2: Genomic Selection Sweep Scan in Canines

## 1. Background and Dataset

The dataset is sourced from the **Dog10K** consortium ([Download Link](https://dog10k.kiz.ac.cn/Home/Download)). The original genomic dataset is a high-coverage phased BCF file containing 1,929 individuals and over 29 million SNPs:
- Original file: `AutoAndXPAR.Dog10K.phased.bcf`
- Metadata table: `dog10K-alignment-sample-table.2022-02-23-v7.txt`

### Sample Selection
To ensure the analysis runs in seconds during class, we selected a biologically relevant subset of **130 individuals** representing body size extremes:

| Group | Dog Breed (Breed.Type) | Number of Samples |
| :--- | :--- | :---: |
| **Small** (61) | Dachshund | 17 |
| | Toy Fox Terrier | 10 |
| | Pomeranian | 8 |
| | Brussels Griffon | 7 |
| | Yorkshire Terrier | 5 |
| | Shih Tzu | 5 |
| | Maltese | 4 |
| | Pekingese | 3 |
| | Chihuahua | 2 |
| **Large** (69) | Saint Bernard | 13 |
| | Leonberger | 11 |
| | Bernese Mountain Dog | 10 |
| | Greater Swiss Mountain Dog | 10 |
| | Great Pyrenees | 7 |
| | Bullmastiff | 6 |
| | Mastiff | 6 |
| | Newfoundland | 6 |

---

## 2. Preprocessing and Filtering

We filter the massive BCF file to include only our 130 samples and chromosome 15, while removing low-frequency SNPs (MAF < 0.05) that are not informative for breed/size differentiation. The coordinates correspond to the CanFam3 reference genome, where the ***IGF1*** gene is located approximately in the **41.2 Mb to 44.5 Mb** region.

### Bash Code Exercise 1: Extraction & Format Conversion
```bash
# 1. Extract chr15 and filter for samples and MAF >= 0.05
bcftools view \
  -S input/Part_2_CanidDiversity/subset_dogs.txt \
  -r chr15 \
  -q 0.05:minor \
  -O b \
  -o input/Part_2_CanidDiversity/subset_chr15.bcf \
  input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf

# 2. Index the subset BCF
bcftools index input/Part_2_CanidDiversity/subset_chr15.bcf

# 3. Convert BCF to PLINK binary format (.bed/.bim/.fam)
plink1.9 \
  --bcf input/Part_2_CanidDiversity/subset_chr15.bcf \
  --dog \
  --keep-allele-order \
  --make-bed \
  --out output/subset_chr15
```

> **Premise**: The filtered BCF file contains **177,953 SNPs** on chromosome 15 across the 130 samples. The coordinates correspond to the CanFam3 reference genome, where the ***IGF1*** gene is located approximately in the **41.2 Mb to 44.5 Mb** region.

---

## 3. Population Structure Analysis (PCA)

We will first run a Principal Component Analysis (PCA) using **PLINK 1.9** to characterize the genetic structure of our subset and visualize it in R.

### R Code Exercise 2: PCA Visualization
```R
library(ggplot2)

# Load eigenvectors and eigenvalues from PLINK output
eigenvec <- read.table("input/Part_2_CanidDiversity/plink_pca.eigenvec", header = FALSE)
eigenval <- read.table("input/Part_2_CanidDiversity/plink_pca.eigenval", header = FALSE)
sample_info <- read.table("input/Part_2_CanidDiversity/sample_info.txt", header = TRUE, sep = "\t")

names(eigenvec)[1:4] <- c("FID", "sampleName", "PC1", "PC2")
pca_data <- merge(sample_info, eigenvec, by = "sampleName")

pc1_var <- round((eigenval$V1[1] / sum(eigenval$V1)) * 100, 2)
pc2_var <- round((eigenval$V1[2] / sum(eigenval$V1)) * 100, 2)

# Plot PC1 vs PC2
ggplot(pca_data, aes(x = PC1, y = PC2, color = breed, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Canine Genotypes PCA (Chromosome 15)",
    x = paste0("PC1 (", pc1_var, "% variance explained)"),
    y = paste0("PC2 (", pc2_var, "% variance explained)"),
    color = "Breed",
    shape = "Body Size Group"
  )
```

### Resulting Figure: PCA Plot
![PCA Plot](figures/figs_canid_diversity/pca_plot.png)

### Questions for Students
1. **What pattern do you observe along the first principal component (PC1)?**
   * *Answer*: PC1 clearly separates small dog breeds (on the right) from large/giant dog breeds (on the left), suggesting that body size differentiation corresponds to the main axis of genetic variation on this chromosome.
2. **Why does PC1 capture body size differences in this particular dataset?**
   * *Answer*: Because we intentionally selected samples from phenotypic extremes (very small vs. very large breeds), the primary artificial genetic differentiation in this specific test dataset is driven by genes associated with body size.

---

## 4. Genomic Outlier Detection using PCAdapt

**PCAdapt** is a method designed to find SNPs that are exceptionally related to population structure (PCs) rather than neutral drift.

### The Role of Linkage Disequilibrium (LD) and Clumping
> [!IMPORTANT]
> **Key Concept**: PCA is highly sensitive to Linkage Disequilibrium (LD). If a region contains many highly correlated markers (due to a selective sweep or low recombination), that single region will dominate the principal components, biasing the PCA and masking other genomic signals (such as the *IGF1* gene sweep).
> 
> To resolve this, we must enable **LD Clumping** in PCAdapt. Thinning out redundant SNPs in strong LD allows the global genomic structure to be correctly computed and helps locate narrow selection sweeps.

### R Code Exercise 3: PCAdapt with Clumping
```R
library(pcadapt)
library(ggplot2)

# Load genotype data
x <- read.pcadapt("input/Part_2_CanidDiversity/subset_chr15.bed", type = "bed")

# Run pcadapt with K = 2 and LD clumping enabled (size = 500, thr = 0.1)
res <- pcadapt(x, K = 2, LD.clumping = list(size = 500, thr = 0.1))

# Load bim file for positions
bim <- read.table("input/Part_2_CanidDiversity/subset_chr15.bim", col.names = c("chr", "snp_id", "cm", "pos", "a1", "a2"))

# Organize results
pcadapt_results <- data.frame(
  Position = bim$pos,
  Pvalue = res$pvalues,
  log10P = -log10(res$pvalues)
)
pcadapt_results <- na.omit(pcadapt_results)

# Highlight the IGF1 region on chr15 (around 43-45 Mb)
igf1_start <- 43e6
igf1_end <- 45e6
pcadapt_results$is_igf1 <- pcadapt_results$Position >= igf1_start & pcadapt_results$Position <= igf1_end

# Get the peak height for annotation
igf1_peak_y <- max(pcadapt_results$log10P[pcadapt_results$is_igf1], na.rm = TRUE)

# Plot Manhattan Plot
ggplot(pcadapt_results, aes(x = Position / 1e6, y = log10P, color = is_igf1)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E74C3C")) +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  # Annotate the IGF1 peak directly on the plot
  annotate("text", x = 43.4, y = igf1_peak_y + 2, label = "IGF1", color = "#E74C3C", fontface = "bold", size = 5) +
  annotate("segment", x = 43.4, xend = 43.4, y = igf1_peak_y + 1.6, yend = igf1_peak_y + 0.4, color = "#E74C3C", arrow = arrow(length = unit(0.15, "cm"))) +
  theme_minimal() +
  labs(
    title = "PCAdapt Outlier Scan (chr15)",
    subtitle = "Peak at ~43.4 Mb highlights the IGF1 gene after LD clumping",
    x = "Physical Position (Mb)",
    y = "-log10(p-value)"
  ) +
  theme(legend.position = "none")
```

### Resulting Figure: PCAdapt Manhattan Plot
![PCAdapt Manhattan](figures/figs_canid_diversity/pcadapt_manhattan.png)

### Questions for Students
1. **Compare PCAdapt results with and without LD clumping (without clumping, a massive peak at ~61 Mb dominates the plot, hiding other regions). What is the effect of LD pruning on outlier detection?**
   * *Answer*: LD clumping removes redundant genetic markers that are physically linked and share the same evolutionary history. Without clumping, regions with long-range LD (like a local breed bottleneck at ~61 Mb) dominate the PCA eigenvalues and mask the actual selective sweep peak at the *IGF1* locus (~43.4 Mb).
2. **Did you detect the IGF1 locus outlier peak? At what physical position is the top marker located, and what is its significance?**
   * *Answer*: Yes, after LD clumping, the *IGF1* locus is detected. The peak is centered at approximately **43.4 Mb** with a highly significant p-value ($p \approx 8.7 \times 10^{-21}$), showing that this locus is strongly associated with the genetic division between small and large dog breeds on chromosome 15.

---

## 5. Cross-Population Selection Scan using XP-nSL

To confirm that the outlier peak on chromosome 15 is indeed driven by a selective sweep, we will perform a haplotype-based selection scan. Specifically, we will run **XP-nSL (Cross-Population Number of Segregating Sites by Length)** to compare the haplotype homozygosity decay between small dogs and large dogs.

### Key Concepts: nSL and XP-nSL
- **nSL**: A within-population selection scan metric similar to iHS. However, instead of measuring haplotype decay in terms of genetic distance (which requires a genetic map), nSL measures distance by counting the number of segregating sites (segregating site count by length). This makes it highly robust to recombination rate variation and suitable for genomes without well-defined genetic maps.
- **XP-nSL**: A cross-population statistic that compares nSL profiles between a target population and a reference population. A high positive score indicates a selective sweep specific to the target population (longer haplotypes around the derived allele).
- **Phased Mode**: Since our input Dog10K BCF file is already phased (containing haplotype data formatted as `0|0`, `1|0`, etc.), we will perform a phased XP-nSL scan. This utilizes the precise haplotype sequences, which provides a significantly stronger selection signal compared to unphased analyses.

### Outgroup Allele Polarization
Haplotype selection scans require knowing which allele is **ancestral** (original) and which is **derived** (new mutant). `selscan` expects a VCF file where `0` is the ancestral allele and `1` is the derived allele.
To polarize our dataset, we use the **gray wolves** in the Dog10K metadata as an outgroup:
- Gray wolves are the evolutionary ancestor of domestic dogs.
- For each SNP, the most common (major) allele in the gray wolf population is designated as the ancestral allele.
- If the ALT allele in the original VCF is the major allele in wolves, we must physically swap the REF/ALT alleles and swap genotypes (`0` becomes `1`, and `1` becomes `0`) for all individuals.

### Preprocessing and Polarization Pipeline
We will run a script to combine our small and large dog samples with wolf samples, extract chromosome 15, polarize the alleles, and split them back into target (small) and reference (large) VCFs.

### Bash Code Exercise 4: Extracting and Polarizing Alleles
```bash
# 1. Extract wolf samples and combine with dog samples
bcftools query -l input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf | grep -E '^CLUP' > input/Part_2_CanidDiversity/wolves.txt
cat input/Part_2_CanidDiversity/subset_dogs.txt input/Part_2_CanidDiversity/wolves.txt > input/Part_2_CanidDiversity/dogs_and_wolves.txt

# 2. Extract polymorphic sites for both dogs and wolves
bcftools query -f '%CHROM\t%POS\n' input/Part_2_CanidDiversity/subset_chr15.bcf > input/Part_2_CanidDiversity/subset_chr15_positions.txt
bcftools view \
  -S input/Part_2_CanidDiversity/dogs_and_wolves.txt \
  -T input/Part_2_CanidDiversity/subset_chr15_positions.txt \
  -O z \
  -o input/Part_2_CanidDiversity/subset_chr15_with_wolves.vcf.gz \
  input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf

# 3. Polarize alleles using Python (major allele in wolves = 0)
python3 scripts/polarize_by_wolves.py

# 4. Re-compress to block gzip format and index polarized VCF
bcftools view input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz -O z -o output/subset_chr15_polarized_bgzf.vcf.gz
mv output/subset_chr15_polarized_bgzf.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz

# 5. Extract polarized VCFs for small and large dogs separately
bcftools view -S input/Part_2_CanidDiversity/small_dogs.txt -O z -o input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz

bcftools view -S input/Part_2_CanidDiversity/large_dogs.txt -O z -o input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz
```

---

### Selection Scan Execution and Normalization
We will execute the selection scan using the compiled `selscan` binary, comparing small dogs (target) against large dogs (reference) in phased mode, and normalize the raw scores.

### Bash Code Exercise 5: Running XP-nSL and Normalizing
```bash
# 1. Run phased XP-nSL scan (selscan automatically runs in phased mode when VCF contains phased '|' alleles)
/home/patriciopezovalderrama/programs/selscan/src/selscan \
  --xpnsl \
  --vcf input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz \
  --vcf-ref input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz \
  --threads 4 \
  --out output/xpnsl_phased

# 2. Joint normalization of raw scores across alleles
/home/patriciopezovalderrama/programs/selscan/src/selscan norm \
  --xpnsl \
  --files output/xpnsl_phased.xpnsl.out

# 3. Calculate window-based statistics (100 Kb non-overlapping windows)
/home/patriciopezovalderrama/programs/selscan/src/selscan norm \
  --xpnsl \
  --files output/xpnsl_phased.xpnsl.out \
  --bp-win \
  --winsize 100000
```

---

### Plotting results in R
We will read the normalized XP-nSL scores (both at the individual SNP level and the 100 Kb window level) and plot them along chromosome 15 to identify target-specific selection peaks.

### R Code Exercise 6: XP-nSL Haplotype Manhattan Plots
```R
library(ggplot2)

input_dir <- "input/Part_2_CanidDiversity"
figures_dir <- "figures/figs_canid_diversity"

# === 1. Plot SNP-level XP-nSL (-log10 P-value) ===
norm_file <- file.path(input_dir, "xpnsl_phased.xpnsl.out.norm")
xpnsl_data <- read.table(norm_file, header = TRUE, sep = "\t")
xpnsl_data <- na.omit(xpnsl_data)

# Calculate two-tailed p-value from standard normal Z-score
xpnsl_data$p_val <- 2 * pnorm(-abs(xpnsl_data$norm_xpnsl))
xpnsl_data$p_val[xpnsl_data$p_val <= 0] <- 1e-15 # prevent zero p-values
xpnsl_data$log10P <- -log10(xpnsl_data$p_val)

# Highlight the IGF1 region (around 41-45.5 Mb)
igf1_start <- 41e6
igf1_end <- 45.5e6
xpnsl_data$is_igf1 <- xpnsl_data$pos >= igf1_start & xpnsl_data$pos <= igf1_end

# Get the peak values for the sweep region
igf1_peak_log10p <- max(xpnsl_data$log10P[xpnsl_data$is_igf1], na.rm = TRUE)
peak_idx <- which(xpnsl_data$is_igf1 & xpnsl_data$log10P == igf1_peak_log10p)[1]
peak_pos <- xpnsl_data$pos[peak_idx] / 1e6

# Plot SNP-level Manhattan
p_snp <- ggplot(xpnsl_data, aes(x = pos / 1e6, y = log10P, color = is_igf1)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E74C3C")) +
  geom_hline(yintercept = 2, color = "blue", linetype = "dashed") + # p = 0.01 threshold
  annotate("text", x = peak_pos, y = igf1_peak_log10p + 0.5, label = "IGF1", color = "#E74C3C", fontface = "bold", size = 5) +
  annotate("segment", x = peak_pos, xend = peak_pos, y = igf1_peak_log10p + 0.4, yend = igf1_peak_log10p + 0.1, color = "#E74C3C", arrow = arrow(length = unit(0.15, "cm"))) +
  theme_minimal() +
  labs(
    title = "Cross-Population Phased XP-nSL Scan: Small vs Large Dogs (chr15)",
    subtitle = "SNP-level phased selection scan significance",
    x = "Physical Position (Mb)",
    y = "-log10(p-value)"
  ) +
  theme(legend.position = "none")

ggsave(file.path(figures_dir, "xpnsl_manhattan.png"), plot = p_snp, width = 10, height = 5, dpi = 300)

# === 2. Plot Window-level (100 Kb) XP-nSL ===
win_file <- file.path(input_dir, "xpnsl_phased.xpnsl.out.norm.100kb.windows")
win_data <- read.table(win_file, header = TRUE, sep = "\t")
win_data <- na.omit(win_data)
win_data <- win_data[win_data$n_snps >= 10 & win_data$frac_top >= 0, ]

win_data$is_igf1 <- win_data$start >= igf1_start & win_data$end <= igf1_end
win_peak_y <- max(win_data$frac_top[win_data$is_igf1], na.rm = TRUE)
win_peak_idx <- which(win_data$is_igf1 & win_data$frac_top == win_peak_y)[1]
win_peak_pos <- (win_data$start[win_peak_idx] + win_data$end[win_peak_idx]) / 2 / 1e6

p_win <- ggplot(win_data, aes(x = (start + end) / 2 / 1e6, y = frac_top, color = is_igf1)) +
  geom_line(color = "#BDC3C7", alpha = 0.8) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E74C3C")) +
  annotate("text", x = win_peak_pos, y = win_peak_y + 0.08, label = "IGF1", color = "#E74C3C", fontface = "bold", size = 5) +
  annotate("segment", x = win_peak_pos, xend = win_peak_pos, y = win_peak_y + 0.06, yend = win_peak_y + 0.015, color = "#E74C3C", arrow = arrow(length = unit(0.15, "cm"))) +
  theme_minimal() +
  labs(
    title = "Window-Based Phased XP-nSL Scan: Small vs Large Dogs (100 Kb Windows)",
    subtitle = "Proportion of extreme positive SNPs (XP-nSL > 2.0) per window",
    x = "Physical Position (Mb)",
    y = "Fraction of Extreme SNPs"
  ) +
  theme(legend.position = "none")

ggsave(file.path(figures_dir, "xpnsl_window_manhattan.png"), plot = p_win, width = 10, height = 5, dpi = 300)
```

### Resulting Figures: XP-nSL Haplotype Manhattan Plots

#### SNP-level XP-nSL Scores ($-log_{10}$ P-value)
![XP-nSL Manhattan](figures/figs_canid_diversity/xpnsl_manhattan.png)

#### 100 Kb Window-level Proportions
![XP-nSL Window Manhattan](figures/figs_canid_diversity/xpnsl_window_manhattan.png)

---

### Questions for Students
1. **Why do we need to polarize alleles using an outgroup like the gray wolf? What does the `0` vs `1` coding represent in `selscan`?**
   * *Answer*: `selscan` uses EHH-based statistics that track the decay of haplotype homozygosity for ancestral and derived alleles separately. In VCF files, alleles are represented as `0` (REF) and `1` (ALT). Without polarization, the program would treat REF as ancestral and ALT as derived by default, which is incorrect. Using the gray wolf (ancestor/outgroup of dogs) allows us to identify the ancestral allele (the major allele in wolves). We swap the alleles and genotypes in the VCF so that `0` is always the ancestral allele, and `1` is always the derived allele.
2. **Why does the raw SNP-level XP-nSL scan look like a noisy cloud at individual sites? What is the effect of calculating window-based scores (e.g. 100 Kb)?**
   * *Answer*: Haplotype selection statistics calculate the homozygosity decay at each variant individually. Because individual SNPs within a sweep region will have varying frequencies and haplotype contexts, scores naturally fluctuate at the single-marker level, creating a cloud of points. Window-based analysis averages the selection signal across a physical block (e.g. 100 Kb) by calculating the proportion of extreme outliers (`frac_top`). This acts as a spatial low-pass filter, smoothing out local stochastic noise and highlighting the selective sweep with exceptional clarity.

---

## 6. Alternative Haplotype Selection Scan: Rsb using `rehh`

As a complementary approach to XP-nSL, we will run **Rsb**, another widely used cross-population EHH-based statistic.

### Key Concepts & Comparison
- **Rsb**: Compares the integrated Extended Haplotype Homozygosity (iHH) between two populations. It is calculated as $\ln(iES_{pop1} / iES_{pop2})$, where $iES$ is the integrated EHH over physical distance (bp). A high positive score indicates selection in population 1 (small dogs), while a negative score indicates selection in population 2 (large dogs).
- **Difference from XP-nSL**: 
  - **XP-nSL** integrates the nSL metric over the number of segregating sites (SNP count). This makes it highly robust to local recombination rate variation.
  - **Rsb** integrates EHH over physical distance (bp). In species with very strong selective sweeps and long-range Linkage Disequilibrium (like domestic dogs), Rsb can produce exceptionally high, clear peaks at sweep loci like *IGF1*.

---

### R Code Exercise 7: Rsb Selection Scan
```R
library(rehh)
library(ggplot2)

input_dir <- "input/Part_2_CanidDiversity"
figures_dir <- "figures/figs_canid_diversity"

# 1. Load polarized target (small) and reference (large) haplotype files
hh_small <- data2haplohh(
  hap_file = file.path(input_dir, "small_dogs_polarized.vcf.gz"),
  chr = "chr15",
  allele_coding = "map",
  polarize_vcf = FALSE,
  verbose = TRUE
)

hh_large <- data2haplohh(
  hap_file = file.path(input_dir, "large_dogs_polarized.vcf.gz"),
  chr = "chr15",
  allele_coding = "map",
  polarize_vcf = FALSE,
  verbose = TRUE
)

# 2. Scan Haplotype Homozygosity
scan_small <- scan_hh(hh_small, polarized = FALSE, threads = 4)
scan_large <- scan_hh(hh_large, polarized = FALSE, threads = 4)

# 3. Compute cross-population Rsb
rsb_results <- ines2rsb(
  scan_small,
  scan_large,
  popname1 = "small",
  popname2 = "large"
)

# 4. Prepare data for plotting
rsb_df <- rsb_results
rsb_df$log10P <- rsb_df$LOGPVALUE
rsb_df <- na.omit(rsb_df)

# Highlight IGF1 region on chr15 (around 41-45.5 Mb)
igf1_start <- 41e6
igf1_end <- 45.5e6
rsb_df$is_igf1 <- rsb_df$POSITION >= igf1_start & rsb_df$POSITION <= igf1_end

# Find the peak position for annotation
igf1_peak_log10p <- max(rsb_df$log10P[rsb_df$is_igf1], na.rm = TRUE)
peak_idx <- which(rsb_df$is_igf1 & rsb_df$log10P == igf1_peak_log10p)[1]
peak_pos <- rsb_df$POSITION[peak_idx] / 1e6

# 5. Plot Rsb Manhattan Plot
p_rsb <- ggplot(rsb_df, aes(x = POSITION / 1e6, y = log10P, color = is_igf1)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E74C3C")) +
  geom_hline(yintercept = 4, color = "blue", linetype = "dashed") +
  annotate("text", x = peak_pos, y = igf1_peak_log10p + 0.5, label = "IGF1", color = "#E74C3C", fontface = "bold", size = 5) +
  annotate("segment", x = peak_pos, xend = peak_pos, y = igf1_peak_log10p + 0.4, yend = igf1_peak_log10p + 0.1, color = "#E74C3C", arrow = arrow(length = unit(0.15, "cm"))) +
  theme_minimal() +
  labs(
    title = "Cross-Population Rsb Scan: Small vs Large Dogs (chr15)",
    subtitle = "Peak highlights the IGF1 body size locus",
    x = "Physical Position (Mb)",
    y = "-log10(p-value)"
  ) +
  theme(legend.position = "none")

ggsave(file.path(figures_dir, "rsb_manhattan.png"), plot = p_rsb, width = 10, height = 6, dpi = 300)
```

### Resulting Figure: Rsb Manhattan Plot
![Rsb Manhattan](figures/figs_canid_diversity/rsb_manhattan.png)

### Questions for Students
1. **Explain the physical and mathematical difference between Rsb and XP-nSL. Why does Rsb show a much higher, less noisy peak at the *IGF1* locus in dogs compared to XP-nSL?**
   * *Answer*: Rsb integrates EHH decay over physical distance (bp), while XP-nSL integrates nSL over the count of segregating sites (SNPs). Because domestic dogs have very long haplotypes around *IGF1* due to strong selective sweeps and breed bottlenecks, EHH remains high over several megabases of physical distance. Integrating EHH over this massive physical span results in an exceptionally high Rsb score. XP-nSL, on the other hand, measures distance in terms of segregating sites, meaning the high density of SNPs within the sweep region keeps the nSL scale moderate, resulting in a lower and more localized peak.
