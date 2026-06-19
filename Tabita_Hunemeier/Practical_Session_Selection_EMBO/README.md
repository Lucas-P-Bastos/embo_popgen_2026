# EMBO Practical Course: Genomic Diversity & Natural Selection Scan (2026)

This repository contains the practical tutorial and worksheet materials for the **EMBO Practical Course on Genomic Diversity and Natural Selection Scan**. 

The tutorial guides students and instructors through the analysis of real genomic datasets to identify signatures of positive natural selection in both human and domestic dog populations.

---

## Tutorial Overview

The course is divided into two main sections:

### Part 1: Human Genomic Diversity and Natural Selection
- **Objective**: Investigate signatures of positive selection in the candidate gene ***EDAR*** (associated with ectodermal traits like hair thickness and dental morphology in East Asians and Native Americans).
- **Methods**: 
  - Population differentiation metrics: Pairwise $F_{ST}$ (Weir & Cockerham) and **Population Branch Statistics (PBS)**.
  - Haplotype-based tests: **Extended Haplotype Homozygosity (EHH)** decay, furcation trees, **Integrated Haplotype Score (iHS)** (SNP-level and sliding window), and cross-population **XP-EHH** tests.
- **Key Population Comparisons**: African (AFR/LWK), European (EUR), East Asian (EAS/CHS), and Native American (NAM/Peruvian) populations from the 1000 Genomes Project Phase III.

### Part 2: Genomic Selection Scan in Canines
- **Objective**: Locate the selective sweep at the ***IGF1*** body-size locus by comparing small vs. large dog breeds.
- **Methods**:
  - **Principal Component Analysis (PCA)** with PLINK 1.9 to examine genetic structure.
  - **PCAdapt** outlier scan with Linkage Disequilibrium (LD) clumping to isolate localized selection signals.
  - Cross-population haplotype scans: **XP-nSL** (Number of Segregating Sites by Length) with Gray Wolves as an outgroup for ancestral allele polarization, and **Rsb** (Extended Haplotype Homozygosity ratio).
- **Key Dataset**: Dog10K consortium phased BCF genotypes for 130 domestic dogs and gray wolves.

---

## Repository Structure

The project directory is structured as follows:

```
├── tutorial_completo.md      # Instructor Version (Markdown)
├── tutorial_completo.html    # Instructor Version (Self-contained HTML)
├── tutorial_alumnos.md       # Student Version (Markdown)
├── tutorial_alumnos.html     # Student Version (Self-contained HTML)
├── custom.css                # Premium styling configuration
├── input/                    # Organized input directories:
│   ├── Part_1_HumanDiversity/   # Human datasets (VCFs, weir.fst)
│   └── Part_2_CanidDiversity/   # Canine datasets (VCFs, PLINK, XP-nSL)
├── figures/                  # Organized output directories:
│   ├── figs_human_diversity/    # Human selection sweep plots
│   └── figs_canid_diversity/    # Canine selective sweep plots
├── Scripts/                  # R and Python scripts running the selection pipelines
├── output/                   # Directory for active run outputs
└── .gitignore                # Git ignore configuration
```

---

## Document Versions

* **Instructor Version (`tutorial_completo`)**: Contains the fully written Bash/R code snippets, detailed explanations, calculated statistics, and complete answers to all student questions.
* **Student Version (`tutorial_alumnos`)**: Contains the background text, questions, and empty code block skeletons (`# Write your code here`) to be filled out during the practical sessions.

### 💡 Fully Self-Contained HTML Documents
The `.html` versions of both the student and instructor sheets have been compiled using Pandoc's `--embed-resources` flag. **All figures, styles, and web fonts are embedded directly into the HTML code as base64 Data URIs.** 
- This means you can copy, email, or move `tutorial_completo.html` or `tutorial_alumnos.html` to any folder or machine independently, and they will load and render perfectly with all figures intact without needing the `Output/` or `Figures/` folders.

---

## Software & Package Requirements

To run the complete pipelines, you will need the following tools:

### Command Line Tools:
- `bcftools` (v1.18+)
- `plink1.9`
- `selscan` (for XP-nSL calculations)
- `pandoc` (to compile markdown documentation)

### R Packages:
- `rehh` (for EHH, iHS, XP-EHH, and Rsb analyses)
- `data.table` (for fast data manipulation)
- `ggplot2` (for custom premium plots)
- `pcadapt` (for population structure outliers)
