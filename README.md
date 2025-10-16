# Mex-acral-exomes-transcriptomes
This repository contains the scripts and analyses necessary to reproduce results and figures from Basurto-Lozada et al, 2024

## Repository content

This repository contains the following files:

* Admixture_Analysis

    This file contains the code used for the admixture analysis on the mexican patients and plot the proportion of different ancestries per patient. 

* Acral_melanoma_exome_somatic_variants

    This file contains the code and tools used to identify mutational driver genes, and the code used to generate Figure 1a (oncoplot), Figure 1b (lollipop plots) and Supplementary Figures 2a and 2b. 

* Acral_melanoma_exome_copy_number_alterations

    This file contains code and tools used to identify significant regions affected by amplifications and deletions as well as the code used to plot and scrutinize homozygous deletions in NF1 and CDKN2A. It also contains the code used to  generate plots for Figure 2a (significant peaks plot), Figure 2b (copy number heatmap), Supplementary Figure 5 (significant peaks oncoplot), and Supplementary Figure 4 (significant peaks of mutated samples vs QWT samples).

* Acral_melanoma_ancestry_vs_mutations

    This file contains the code used for the association analysis between driver mutational status and ancestry and the code to plot Figure 1c. 

* Acral_melanoma_exome_correlation_analysis

    This file contains the code used to plot and compare the burden of CN alterations by driver mutational status (Figure 2c), the code used to plot and compare the burden of CN alterations between different anatomical sites (Figure 2e) and the code to see the correlation between the burden of CN alterations and TMB (Figure 2d). It also contains code to compare and plot TMB by druver mutational status (Supplementary Figure 6) and to compare the proportion of amerindian ancestry by driver mutational status (Supplementary Figure 3).

* Acral_melanoma_consensus_clustering

    This file contains the code to generate the consesus clustering for RNAseq data and to generate the plot for Figure 4a.

* Cell_types_proportion_by_RNA_cluster

    This file contains the code to compare and plot the presence of different cell types by RNA cluster. 

* Acral and cutaneous comparison Mexican data

    This file contains the code used to generate and compare acral-cutaneous scores between samples with different BRAF mutational status in the mmexican cohort.

* Acral and cutaneous comparison Newell

    This file contains the code used to generate and compare acral-cutaneous scores between samples with different BRAF mutational using data from Newell et al (2020).

* Acral_melanoma_vs_cutaneous_melanomas_comparison

    This file contains the explanation of the data and tests used to select genes for the generation of acral-cutaneous (A:C) scores. 

* Acral_melanoma_exomes_AC_scores_comparison

    This file contains the code used to visualize the comparison of acral-cutaneous scores by BRAF mutational status (Figure 3c), the code to visualize the comparison of normalized expression of classifier genes of melanocytres with induced BRAFV600E and the comparison of A:C scores using TCGA data by BRAF mutational status. 

* Survival_analysis

    This file contains the code used for overall survival and recurrence free survival analysis. 

* Mutational_signature_analysis

    This file contains the code used to identify and plot SBS and CNV signatures present in tumours of the mexican cohort. 
