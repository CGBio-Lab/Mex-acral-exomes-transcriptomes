# Acral melanoma exomes' analysis
# Variant Calling and driver gene identification 
## Author: Patricia Basurto Lozada

Reference genome used: GRCh38

SNV calling was done using three different tools (Mutect, MuTect2 and Varscan2). Varscan2 and MuTect were used following the respective documentation recomendations. Mutect2 was used using a panel of normals created with all normal samples available from the patients fo this cohort and following GATK recommendations. SNVs that were identified by at least two of these tools were used for further analysis.

Indels were identified using the tool Strelka2, using candidate indels indetified by the structural variant caller manta. 

SNVs and indels used for further analysis were stored as MAF files, one file per sample (vcf files were converted to MAF using vcf2maf). 

# Driver gene identification

Driver gene identification was performed using the tool dndscv using only unique events (snv and indels) per patient as input data (mutations present in different samples from same patients were reported only once). Genes were considered to be driver genes when q-values estimated by dndscv were below 0.1. 

``` R
#Loading library
library(dndscv)
#Reading input data
acral_data <- read.csv("filtered_unique_snp_strelka_indels.csv", header=TRUE, sep='\t')
#Running dndscv
dndsout = dndscv(acral_data, refdb = "RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=NULL)
#Showing significant genes
sel_cv = dndsout$sel_cv
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)
  gene_name   qglobal_cv
1      NRAS 0.000000e+00
2       KIT 4.327078e-08
3      BRAF 1.964711e-06
4       NF1 5.768489e-02

```
# Visualization of mutational profile as an oncoplot

Visualization of the mutational profile was done generating an oncoplot using the R package maftools. One sample per patient was used, giving preference to primary tumours when available. Genes showed were selected based on dndscv results. 

``` R
library(maftools)
library(RColorBrewer)
#Reading mafs containing SNVs (folder contains only one maf per patient)
mafs = list.files(path = "one_per_patient_snp_maf/", pattern = "*.\\.csv$", full.names =TRUE)
#Merging  SNV mafs
all_samples = merge_mafs(mafs = mafs)
#Reading mafs containing indels (folder contains only one maf per patient)
indel_mafs = list.files(path = "one_per_patient_strelka_mafs/", pattern = "*.\\.maf$", full.names =TRUE)
#Merging indek mafs
indel = merge_mafs(mafs = indel_mafs)
#Merging SNV and indel mafs
all_alt = merge_mafs(mafs= c(indel, all_samples))
#Reading file containing clinical data
clindata = read.csv("Supplementary_table_1.csv")

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Accent')
#Matching colors with categories
names(clincolors) = c("FEET", "HAND", "SUBUNGUAL")
clincolors = list(Site = clincolors)

#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Site"), bgCol="white")


#Plotting other annotations

#Gender

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 3,name = 'Dark2')
#Matching colors with categories
names(clincolors) = c("female", "male", "NA")
clincolors = list(Gender = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Gender"), bgCol="white")


#Age

#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Age"), bgCol="white")

#Ulceration

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Set1')
#Matching colors with categories
names(clincolors) = c("yes", "no", "")
clincolors = list(Ulceration = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Ulceration"), bgCol="white")

#Tumor_type
#For this annotation, tumor type was simplified in four categories: primary, metastasis, recurrence and lesion in transit
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Paired')
names(clincolors) = c("primary", "metastasis", "Recurrence","Lesion_in_transit")
clincolors = list(Tumor_type = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Tumor_type"), bgCol="white")


#Stage

#Reading stage data as factors
clindata$Stage <- as.factor(clindata$Stage)
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 5,name = 'BuPu')
names(clincolors) = c("0", "1", "2","3","4")
clincolors = list(Stage = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Stage"), bgCol="white")
```
## Visualizing mutations in specific genes

Mutations affecting genes that were considered relevant were visualized generating lollipop plots using maftools.

``` R
#Lollipop plots
lollipopPlot(maf=all_alt, gene="NRAS",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)
lollipopPlot(maf=all_alt, gene="BRAF",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)
lollipopPlot(maf=all_alt, gene="NF1",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)
lollipopPlot(maf=all_alt, gene="KIT",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)

```
# PLotting oncoplot for sample with no drivers

We generated an oncoplot using one sample per patient (giving priority to primary tumors when available) of all patients that did not have mutations in the genes NRAS, KIT, BRAF and NF1.

``` R
#Loading packages
library(maftools)
library(RColorBrewer)
#Reading MAFs with SNV data
mafs = list.files(path = "one_per_patient_non_driver_maf/", pattern = "*.\\.csv$", full.names =TRUE)
#Merging  SNV mafs
all_samples = merge_mafs(mafs = mafs)
#Reading MAFs with indel data
indel_mafs = list.files(path = "one_per_patient_strelka_non_driver_mafs/", pattern = "*.\\.maf$", full.names =TRUE)
#Merging indel mafs
indel = merge_mafs(mafs = indel_mafs)
#Merging SNV and indel mafs
all_alt = merge_mafs(mafs= c(indel, all_samples))
#Reading file containing clinical data
clindata = read.csv("Supplementary_table_1.csv")

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Accent')
names(clincolors) = c("FEET", "HAND", "SUBUNGUAL")
clincolors = list(Site = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Site"), bgCol="white")


#Gender

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 3,name = 'Dark2')
names(clincolors) = c("female", "male", "NA")
clincolors = list(Gender = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=TRUE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Gender"), bgCol="white", SampleNamefontSize=0.8)


#Age

#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Age"), bgCol="white", SampleNamefontSize=0.8)


#Ulceration

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Set1')
names(clincolors) = c("yes", "no", "")
clincolors = list(Ulceration = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Ulceration"), bgCol="white")


#Tumor_type

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Paired')
names(clincolors) = c("primary", "metastasis", "Recurrence","Lesion_in_transit")
clincolors = list(Tumor_type = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Tumor_type"), bgCol="white")


#Stage

#Choosing colors
clindata$Stage <- as.factor(clindata$Stage)
clincolors = RColorBrewer::brewer.pal(n = 5,name = 'BuPu')
names(clincolors) = c("0", "1", "2","3","4")
clincolors = list(Stage = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Stage"), bgCol="white")
```


## Lollipop plot for low frequency potential drivers (All patients , one sample per patient)

We generated lollipop plots for low frequency potential drivers using all patietns, but using one sample per patient. 

``` R
#Reading SNV mafs
mafs = list.files(path = "one_per_patient_snv_maf/", pattern = "*.\\.csv$", full.names =TRUE)
#Merging SNV mafs
all_samples = merge_mafs(mafs = mafs)
#Reading indel mafs
indel_mafs = list.files(path = "one_per_patient_strelka_mafs/", pattern = "*.\\.maf$", full.names =TRUE)
#Merging indel mafs
indel = merge_mafs(mafs = indel_mafs)
#Merging SNV and indel mafs
all_alt = merge_mafs(mafs= c(indel, all_samples))

#Plotting lollipo plots

#HRAS

lollipopPlot(maf=all_alt, gene="HRAS",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)

#SPHKAP

lollipopPlot(maf=all_alt, gene="SPHKAP",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)

#POU3F3

lollipopPlot(maf=all_alt, gene="POU3F3",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)
```

# Notes

Notes
Sample PD51928d (Lymph node metastasis) was annotated as not having a mutation in NRAS by variant calling tools, although there are reads supporting a mutation in this gene, due to this mutation being annotated as normal artefact. Normal adjacent tissue, which was used as normal, had presence of this same mutation in several reads.
Four patients that had mutations in KIT in their primary tumours are reported here as not having this mutation in their lymph node metastasis. Manual revision and comparison of these samples showed that lymph node metastasis samples had lower cellularity compared to primary samples, as well as lower coverage, which could be the cause that KIT mutations were not detected on these samples.

# References 

Varscan2
Koboldt, D. C., Zhang, Q., Larson, D. E., Shen, D., McLellan, M. D., Lin, L., Miller, C. A., Mardis, E. R., Ding, L., & Wilson, R. K. (2012). VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome research, 22(3), 568–576. https://doi.org/10.1101/gr.129684.111
http://dkoboldt.github.io/varscan/

Mutect2
https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2

Mutect
Cibulskis, K., Lawrence, M., Carter, S. et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnol 31, 213–219 (2013). https://doi.org/10.1038/nbt.2514

Strleka2
Kim, S., Scheffler, K., Halpern, A. L., Bekritsky, M. A., Noh, E., Källberg, M., Chen, X., Kim, Y., Beyter, D., Krusche, P., & Saunders, C. T. (2018). Strelka2: fast and accurate calling of germline and somatic variants. Nature methods, 15(8), 591–594. https://doi.org/10.1038/s41592-018-0051-x
https://github.com/Illumina/strelka

manta 
Chen, X., Schulz-Trieglaff, O., Shaw, R., Barnes, B., Schlesinger, F., Källberg, M., Cox, A. J., Kruglyak, S., & Saunders, C. T. (2016). Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics (Oxford, England), 32(8), 1220–1222. https://doi.org/10.1093/bioinformatics/btv710
https://github.com/Illumina/manta

vcf2maf
https://github.com/mskcc/vcf2maf

dndscv
Martincorena, I., Raine, K. M., Gerstung, M., Dawson, K. J., Haase, K., Van Loo, P., Davies, H., Stratton, M. R., & Campbell, P. J. (2017). Universal Patterns of Selection in Cancer and Somatic Tissues. Cell, 171(5), 1029–1041.e21. https://doi.org/10.1016/j.cell.2017.09.042
github: https://github.com/im3sanger/dndscv?tab=readme-ov-file

maftools
Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Research. PMID: 30341162
github: https://github.com/PoisonAlien/maftools
