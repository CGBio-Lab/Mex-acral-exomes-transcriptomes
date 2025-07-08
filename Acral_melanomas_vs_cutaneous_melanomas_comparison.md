# Acral melanoma exomes' analysis
# Acral melanomas vs cutaneous melanomas comparison and generation of A:C classifier


## Comparison of acral and cutaneous melanomas

As reported on the paper, acral and non-acral cutaneous melanomas were collected from which gene expression data was obtained through a custom NanoString nCounter XT CodeSet that included genes differentially expressed between glabrous and non-glabrous melanocytes. Log2 normalised gene expression data was used for a Principal Component Analysis (PCA) using PCA function in PRISM version 10.2.1. 

PCA component values for the samples are available in the data/AC_scores folder. (AM_CM_.csv)

After determining the top differentially expressed genes (AM and CM genes), log2 expression values of the genes were used to generate a multiplicative score and the produce a ratio of acral to cutaneous melanocyte genes. To compare ratios between acral and cutaneous samples, statistical analyses were performed using Prims version 10.2.1, using the Mann-Whitney U test. Visualization of this data was also done using Prism.

Gene expression values for this samples are vailable in the data/AC_scores folder. (PC1_PC2_AM_CM.csv)