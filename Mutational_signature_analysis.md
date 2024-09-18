# Acral melanoma exomes' analysis
# Mutational signature analysis 
## Author: Patricia Basurto Lozada


# Mutational signature analysis

For this analysis we performed a first run and eliminated for firther analyisis all samples that had 50% or more of its mutations associated with artefactual signatures. 

SBS mutational signature analysis (filered samples)

``` python
#Using individual samples, filtering samples with artefacts
source activate sigprofiler1.1.23
python3

from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', bash=True)

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matrices = matGen.SigProfilerMatrixGeneratorFunc("ALM", "GRCh38", "/filtered_samples/",plot=True, exome=False, bed_file="Allexon_v5_Regions.bed", chrom_based=False, tsb_stat=False, seqInfo=True, cushion=100)

from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix", "results", "/output/SBS/ALM.SBS96.region", reference_genome="GRCh38", minimum_signatures=1, maximum_signatures=5, nmf_replicates=100, cpu=-1)
```

Indels

``` python
Extracting indels 
from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix", "results_ID", "/output/ID/ALM.ID83.all", reference_genome="GRCh38", minimum_signatures=1, maximum_signatures=5, nmf_replicates=100, cpu=-1)

```

CNV signatures

```python

from SigProfilerMatrixGenerator.scripts import CNVMatrixGeneratorFunc as scna
file_type="SEQUENZA"
input_file="alm_sequenza_segmentation.csv"
output_path="output"
project="ALM_CNV"
scna.generateCNVMatrix(file_type, input_file, project, output_path)

from SigProfilerExtractor import sigpro as sig
sigProfilerExtractor(“matrix”, “output”, “ALM.CNV48.matrix.tsv”, reference_genome="GRCh38",minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=-1)

```


# Plotting mutational signatures from alm tumours

Generate relative values for mutational signatures

``` python
#For SBS signatures
sbs_activities = pd.read_csv("COSMIC_SBS96_Activities.txt", sep="\t")
sbs_activities["total"] = sbs_activities["SBS1"] + sbs_activities["SBS5"] + sbs_activities["SBS7a"] + sbs_activities["SBS7b"] + sbs_activities["SBS45"]
sbs_activities["SBS1_rel"] = (sbs_activities['SBS1']/(sbs_activities["total"]))*100 
sbs_activities["SBS5_rel"] = (sbs_activities['SBS5']/(sbs_activities["total"]))*100
sbs_activities["SBS7a_rel"] = (sbs_activities['SBS7a']/(sbs_activities["total"]))*100
sbs_activities["SBS7b_rel"] = (sbs_activities['SBS7b']/(sbs_activities["total"]))*100
sbs_activities["SBS45_rel"] = (sbs_activities['SBS45']/(sbs_activities["total"]))*100 
#Generating dataframe with only relative values
rel_sbs = sbs_activities[["Samples","SBS1_rel","SBS5_rel","SBS7a_rel","SBS7b_rel","SBS45_rel"]]
#Saving dataframe with relative activties
rel_sbs.to_csv("Relative_activites_sbs.csv", sep="," , index=False)

#For indels
indel_activities = pd.read_csv("COSMIC_ID83_Activities.txt", sep="\t")
indel_activities["total"] = indel_activities["ID2"] + indel_activities["ID12"]
indel_activities["ID2_rel"] = (indel_activities['ID2']/(indel_activities["total"]))*100
indel_activities["ID12_rel"] = (indel_activities['ID12']/(indel_activities["total"]))*100
#Generating dataframe with only relative values
rel_indel = indel_activities[["Samples","ID2_rel","ID12_rel"]]
#Saving dataframe with relative activties
rel_indel.to_csv("Relative_activites_indel.csv", sep="," , index=False)

#For CNV signatures 
cnv_activities = pd.read_csv("/COSMIC_CNV48_Activities.txt", sep="\t")
cnv_activities["total"] = cnv_activities["CN1"] + cnv_activities["CN2"] + cnv_activities["CN6"] + cnv_activities["CN7"] + cnv_activities["CN9"] + cnv_activities["CN13"] + cnv_activities["CN17"] + cnv_activities["CN20"] + cnv_activities["CN23"] + cnv_activities["CN24"] + cnv_activities["CN25"]
cnv_activities["CN1_rel"] = (cnv_activities['CN1']/(cnv_activities["total"]))*100
cnv_activities["CN2_rel"] = (cnv_activities['CN2']/(cnv_activities["total"]))*100
cnv_activities["CN6_rel"] = (cnv_activities['CN6']/(cnv_activities["total"]))*100
cnv_activities["CN7_rel"] = (cnv_activities['CN7']/(cnv_activities["total"]))*100
cnv_activities["CN9_rel"] = (cnv_activities['CN9']/(cnv_activities["total"]))*100
cnv_activities["CN13_rel"] = (cnv_activities['CN13']/(cnv_activities["total"]))*100
cnv_activities["CN17_rel"] = (cnv_activities['CN17']/(cnv_activities["total"]))*100
cnv_activities["CN20_rel"] = (cnv_activities['CN20']/(cnv_activities["total"]))*100
cnv_activities["CN23_rel"] = (cnv_activities['CN23']/(cnv_activities["total"]))*100
cnv_activities["CN24_rel"] = (cnv_activities['CN24']/(cnv_activities["total"]))*100
cnv_activities["CN25_rel"] = (cnv_activities['CN25']/(cnv_activities["total"]))*100
#Generating dataframe with only relative values
rel_cnv = cnv_activities[["Samples","CN1_rel","CN2_rel","CN6_rel","CN7_rel","CN9_rel", "CN13_rel", "CN17_rel", "CN20_rel", "CN23_rel", "CN24_rel", "CN25_rel"]]
#Saving dataframe with relative activties
rel_cnv.to_csv("Relative_activites_cnv.csv", sep="," , index=False)

```
Generating order of samples for plotting

``` python
#Ordering by percentage of UV signature
UV_ordered_sigs = rel_sbs.sort_values(by=['SBS7b_rel','SBS7a_rel','SBS5_rel','SBS40a_rel', 'SBS1_rel'], ascending=False)
#Getting sample order
sample_order_uvsig = UV_ordered_sigs['Samples']
#Saving sample order 
sample_order_uvsig.to_csv("UV_sample_order.csv", index=False)

```
Plotting the signatures

``` R
#SBS signatures
library(tidyr)
library(ggplot2)
rel_sbs <- read.csv("Relative_activites_sbs.csv", sep=",", header=TRUE)
rel_sbs_pivot = pivot_longer(rel_sbs, cols=2:6, names_to="signature", values_to="count")
p <- p <- ggplot(sbs_mex, aes(x=factor(Samples, level=c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a")), y=count, fill=signature)) + geom_bar(stat="identity", colour="#808080") + scale_fill_brewer(palette="Pastel1")

p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
dev.off()

#Indels 
rel_indel <- read.csv("Relative_activites_indel.csv", sep=",", header=TRUE)
rel_indel_pivot = pivot_longer(rel_indel, cols=2:3, names_to="signature", values_to="count")
indel_mex <- as.data.frame(rel_indel_pivot)

#Ordered bases on sbs (UVfirst)
p <- ggplot(indel_mex, aes(x=factor(Samples, level=c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a")), y=count, fill=signature)) + geom_bar(stat="identity", colour="#808080") + scale_fill_brewer(palette="Pastel2")


p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle=90))
dev.off()

#CNV signatures 

rel_cnv <- read.csv("Relative_activites_cnv.csv", sep=",", header=TRUE)
rel_cnv_pivot = pivot_longer(rel_cnv, cols=2:12, names_to="signature", values_to="count")
cnv_mex <- as.data.frame(rel_cnv_pivot)
# To fet order of samples: sed -z 's/\n/","/g' CN_sample_order.csv 
p <- ggplot(cnv_mex, aes(x=factor(Samples, level=c("PD41026d","PD41917a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD41020a","PD41906a","PD41001a","PD41039a","PD41020d","PD41920c","PD41043a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40988a","PD40972a","PD40990a","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40952a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40970a","PD41921a","PD40968a","PD40961d","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41896a","PD41039d","PD51928a","PD41023f","PD51939a","PD41901a","PD40956d","PD51952a","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD51948a","PD41929d","PD41923f","PD51956a","PD51972a","PD41923g","PD51929a","PD41939d","PD41928d","PD51993a","PD51932a","PD51928d","PD41930d","PD41927d","PD41932d")), y=count, fill=signature)) + geom_bar(stat="identity", colour="gray") + scale_fill_brewer(palette="Paired")

p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle=90))
dev.off()
```

Plotting annotation (tumor mutational burden)

``` R
library(tidyr)
library(ggplot2)
tmb_data <- read.csv("Filtered_correct_snp_indels_strelka_2of3.csv", sep=",", header=TRUE)
q <- ggplot(tmb_data, aes(x=factor(Tumor_Sample_Barcode, level=c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a")), y=snp.indels)) + geom_bar(stat="identity", fill="light blue") 
q + theme(panel.background = element_rect(fill="white"), axis.text.x = element_text(angle=90))
dev.off()

#Plotting other annotations

library(maftools)
library(RColorBrewer)
mafs = list.files(path = "/all_maf/", pattern = "*.\\.csv$", full.names =TRUE)
all_samples = merge_mafs(mafs = mafs)
indel_mafs = list.files(path = "/Strelka_mafs/", pattern = "*.\\.maf$", full.names =TRUE)
indel = merge_mafs(mafs = indel_mafs)
all_alt = merge_mafs(mafs= c(indel, all_samples))
clindata = read.csv("/Supplementary_Table_1.csv")

#Site
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Accent')
names(clincolors) = c("FEET", "HAND", "SUBUNGUAL")
clincolors = list(Site = clincolors)
oncoplot(all_alt, draw_titv = FALSE,removeNonMutated=FALSE, sampleOrder = c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a"), showTumorSampleBarcodes=TRUE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), annotationDat=clindata, annotationColor = clincolors, genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Site"))
dev.off()

#Sex
clincolors = RColorBrewer::brewer.pal(n = 3,name = 'Dark2')
names(clincolors) = c("female", "male", "NA")
clincolors = list(Gender = clincolors)

oncoplot(all_alt, draw_titv = FALSE,removeNonMutated=FALSE, sampleOrder = c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a"), showTumorSampleBarcodes=TRUE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), annotationDat=clindata, annotationColor = clincolors, genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Gender"))
dev.off()

#Age


oncoplot(all_alt, draw_titv = FALSE,removeNonMutated=FALSE, sampleOrder = c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a"), showTumorSampleBarcodes=TRUE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), annotationDat=clindata, , genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Age"))
dev.off()

#Ulceration

clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Set1')
names(clincolors) = c("yes", "no", "")
clincolors = list(Ulceration = clincolors)

oncoplot(all_alt, draw_titv = FALSE,removeNonMutated=FALSE, sampleOrder = c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a"), showTumorSampleBarcodes=TRUE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), annotationDat=clindata, annotationColor = clincolors, genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Ulceration"))
dev.off()

#Tumor_type

clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Paired')
names(clincolors) = c("primary", "metastasis", "Recurrence","Lesion_in_transit")
clincolors = list(Tumor_type = clincolors)

oncoplot(all_alt, draw_titv = FALSE,removeNonMutated=FALSE, sampleOrder = c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a"), showTumorSampleBarcodes=TRUE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), annotationDat=clindata, annotationColor = clincolors, genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Tumor_type"))
dev.off()


#Stage


clindata$Stage <- as.factor(clindata$Stage)
clincolors = RColorBrewer::brewer.pal(n = 5,name = 'BuPu')
names(clincolors) = c("0", "1", "2","3","4")
clincolors = list(Stage = clincolors)

oncoplot(all_alt, draw_titv = TRUE, removeNonMutated=FALSE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE,sampleOrder = c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a"), genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Stage"), bgCol="white")
dev.off()

#mut_status

clincolors = c("#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "gray", "#7ec4cf")
names(clincolors) = c("BRAF", "NRAS", "multihit", "NF1", "Other", "KIT")
clincolors = list(Mutation_status = clincolors)

oncoplot(all_alt, draw_titv = FALSE,removeNonMutated=FALSE, sampleOrder = c("PD40966d","PD41925d","PD41026d","PD41917a","PD41045a","PD41916a","PD41910a","PD40996a","PD41046e","PD41046d","PD41033a","PD41000a","PD41895a","PD40982a","PD51940d","PD40964a","PD41020a","PD41906a","PD40986a","PD41001a","PD41039a","PD41020d","PD41920c","PD40984a","PD51940a","PD41905a","PD41004a","PD41912a","PD41043a","PD41900a","PD41035a","PD41030a","PD41915a","PD41038a","PD41021d","PD40966e","PD41910c","PD40988a","PD40963a","PD40965a","PD40972a","PD40983d","PD40990a","PD40989a","PD41029e","PD40985a","PD41002d","PD41909a","PD41002a","PD41002e","PD40969d","PD40983e","PD40980d","PD40997a","PD40986d","PD41044a","PD41907a","PD41915c","PD40971d","PD40952a","PD40983a","PD41011a","PD40962a","PD40969a","PD40994a","PD40978a","PD40973a","PD40970a","PD41921a","PD40968a","PD40961d","PD40980a","PD40965f","PD40961e","PD40971a","PD41017a","PD41025a","PD40961a","PD41046a","PD41036a","PD40966a","PD40978d","PD41896a","PD40974a","PD41025d","PD41043d","PD40976a","PD41039d","PD51928a","PD41023f","PD51939a","PD41032a","PD41901a","PD40956d","PD41913e","PD51952a","PD40967d","PD41035d","PD41920a","PD40957a","PD41027a","PD41913a","PD40987d","PD51969a","PD40987a","PD40965e","PD40967a"), showTumorSampleBarcodes=TRUE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), annotationDat=clindata, annotationColor = clincolors, genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Mutation_status"))
dev.off()


```
References:

Sigprofiler Matrix Generator
https://github.com/AlexandrovLab/SigProfilerMatrixGenerator

Sigprofiler Extractor
https://github.com/AlexandrovLab/SigProfilerExtractor 

