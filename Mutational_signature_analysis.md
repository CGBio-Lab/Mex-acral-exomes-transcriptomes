# Acral melanoma exomes' analysis
# Mutational signature analysis 
## Author: Patricia Basurto Lozada


# SNV mutational signature analysis

For this analysis we used 116 samples that had a snv count higher than 0.


``` python
source activate sigprofiler1.1.23
python3

from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', bash=True)

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matrices = matGen.SigProfilerMatrixGeneratorFunc("Mex_am", "GRCh38", "/data/Mutational_signatures/",plot=True, exome=True, bed_file="/data/Allexon_v5_Regions.bed", chrom_based=False, tsb_stat=False, seqInfo=True, cushion=100)


from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix", "results", "/data/mut_signatures/output/SBS/Mex_am.SBS96.exome", reference_genome="GRCh38", minimum_signatures=1, maximum_signatures=5, exome=True, nmf_replicates=100, cpu=-1)

```

# CNV signatures

```python

from SigProfilerMatrixGenerator.scripts import CNVMatrixGeneratorFunc as scna
file_type="ASCAT"
input_file="data/Mutational_signatures/acral_ascat_pass_segments.txt"
output_path="cn_matrix_generator_output/"
project="Acral_cn_sigs"
scna.generateCNVMatrix(file_type, input_file, project, output_path)

from SigProfilerExtractor import sigpro as sig
sigProfilerExtractor(“matrix”, “ascat_results”, “Acral_cn_sigs.CNV48.matrix.tsv”, reference_genome="GRCh38",minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=-1)

```


# Plotting mutational signatures from alm tumours (Supplemetary Figure 7)

Generate relative values for mutational signatures

``` python
#For SBS signatures
sbs_activities = pd.read_csv("data/Mutational_signatures/COSMIC_SBS96_Activities.txt", sep="\t")
#Generating relative values per sample
sbs_activities["total"] = sbs_activities["SBS1"] + sbs_activities["SBS5"] + sbs_activities["SBS40a"]
sbs_activities["SBS1_rel"] = (sbs_activities['SBS1']/(sbs_activities["total"]))*100 
sbs_activities["SBS5_rel"] = (sbs_activities['SBS5']/(sbs_activities["total"]))*100
sbs_activities["SBS40a_rel"] = (sbs_activities['SBS40a']/(sbs_activities["total"]))*100
#Generating dataframe with only relative values
rel_sbs = sbs_activities[["Samples","SBS1_rel","SBS5_rel","SBS40a_rel"]]
#Saving dataframe with relative activties
rel_sbs.to_csv("data/Mutational_signatures/Relative_activites_sbs.csv", sep="," , index=False)

#For CNV signatures 
cnv_activities = pd.read_csv("data/Mutational_signatures/COSMIC_CNV48_Activities.txt", sep="\t")
#Generating relative values per sample
cnv_activities["total"] = cnv_activities["CN1"] + cnv_activities["CN2"] + cnv_activities["CN8"] + cnv_activities["CN7"] + cnv_activities["CN9"] + cnv_activities["CN10"] + cnv_activities["CN13"] + cnv_activities["CN17"] + cnv_activities["CN19"] + cnv_activities["CN20"] + cnv_activities["CNV48F"]
cnv_activities["CN1_rel"] = (cnv_activities['CN1']/(cnv_activities["total"]))*100
cnv_activities["CN2_rel"] = (cnv_activities['CN2']/(cnv_activities["total"]))*100
cnv_activities["CN7_rel"] = (cnv_activities['CN7']/(cnv_activities["total"]))*100
cnv_activities["CN8_rel"] = (cnv_activities['CN8']/(cnv_activities["total"]))*100
cnv_activities["CN9_rel"] = (cnv_activities['CN9']/(cnv_activities["total"]))*100
cnv_activities["CN10_rel"] = (cnv_activities['CN10']/(cnv_activities["total"]))*100
cnv_activities["CN13_rel"] = (cnv_activities['CN13']/(cnv_activities["total"]))*100
cnv_activities["CN17_rel"] = (cnv_activities['CN17']/(cnv_activities["total"]))*100
cnv_activities["CN19_rel"] = (cnv_activities['CN19']/(cnv_activities["total"]))*100
cnv_activities["CN20_rel"] = (cnv_activities['CN20']/(cnv_activities["total"]))*100
cnv_activities["CNV48F_rel"] = (cnv_activities['CNV48F']/(cnv_activities["total"]))*100
#Generating dataframe with only relative values
rel_cnv = cnv_activities[["Samples","CN1_rel","CN2_rel","CN7_rel","CN8_rel","CN9_rel", "CN10_rel", "CN13_rel", "CN17_rel", "CN19_rel", "CN20_rel", "CNV48F_rel"]]
#Saving dataframe with relative activties
rel_cnv.to_csv("data/Mutational_signatures/Relative_activites_cnv.csv", sep="," , index=False)


```
Generating order of samples for plotting

``` python
#Ordering by percentage 
ordered_sigs = rel_sbs.sort_values(by=['SBS40a_rel','SBS1_rel','SBS5_rel'], ascending=False)
#Getting sample order
sample_order = ordered_sigs['Samples']
#Saving sample order 
sample_order.to_csv("sample_order.csv", index=False)

```
Plotting the signatures

``` R
#SBS signatures
library(tidyr)
library(ggplot2)
rel_sbs <- read.csv("data/Mutational_signatures/Relative_activites_sbs.csv", sep=",", header=TRUE)
rel_sbs_pivot = pivot_longer(rel_sbs, cols=2:4, names_to="signature", values_to="count")
p <- ggplot(rel_sbs_pivot, aes(x=factor(Samples, level=c('PD40965f','PD40974a','PD40976a','PD40978d','PD40980a','PD40997a','PD41004a','PD41017a','PD41020a','PD41023f','PD41025a','PD41026d','PD41029e','PD41030a','PD41032a','PD41033a','PD41039a','PD41039d','PD41043d','PD41044a','PD41045a','PD41046d','PD41895a','PD41900a','PD41901a','PD41909a','PD41910a','PD41915a','PD41915c','PD41916a','PD41921a','PD41923g','PD51928a','PD51928d','PD51939a','PD51940d','PD51969a','PD51978a','PD51979a','PD41046a','PD41038a','PD41046e','PD40967d','PD40956d','PD41932d','PD41920c','PD40971d','PD41928d','PD41920a','PD41896a','PD40971a','PD40961d','PD40961e','PD41035a','PD51969d','PD41001a','PD40961a','PD41035d','PD40978a','PD41913a','PD41021d','PD40987a','PD41913e','PD51952a','PD41036a','PD40987d','PD40986a','PD40957a','PD40970a','PD40994a','PD40986d','PD40964a','PD40966a','PD41020d','PD40952a','PD51930a','PD51932a','PD41907a','PD40962a','PD40973a','PD41027a','PD41930d','PD41927d','PD40980d','PD40972a','PD41011a','PD41043a','PD40982a','PD51929a','PD41002a','PD41002d','PD40989a','PD40983a','PD40996a','PD51972a','PD41000a','PD51951a','PD41002e','PD40985a','PD40983e','PD40969d','PD40983d','PD40990a','PD41912a','PD51993a','PD40969a','PD41939d','PD40963a','PD40965a','PD40966d','PD40968a','PD40984a','PD40988a','PD41905a','PD41906a','PD51982a')), y=count, fill=signature)) + geom_bar(stat="identity", colour="#808080") + scale_fill_brewer(palette="Pastel1")

p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), )


#CNV signatures 

rel_cnv <- read.csv("data/Mutational_signatures/Relative_activites_cnv.csv", sep=",", header=TRUE)
rel_cnv_pivot = pivot_longer(rel_cnv, cols=2:12, names_to="signature", values_to="count")
cnv_mex <- as.data.frame(rel_cnv_pivot)
p <- ggplot(cnv_mex, aes(x=factor(Samples, level=c('PD40965f','PD40978d','PD40997a','PD41020a','PD41025a','PD41026d','PD41030a','PD41043d','PD41044a','PD41045a','PD41046d','PD41901a','PD41909a','PD41910a','PD41915a','PD41915c','PD41923g','PD51928a','PD51939a','PD51969a','PD51978a','PD41038a','PD41046e','PD40956d','PD41896a','PD40971a','PD40961d','PD40961e','PD41035a','PD40961a','PD41035d','PD41913a','PD41913e','PD40957a','PD40986d','PD41020d','PD40952a','PD51932a','PD40962a','PD41027a','PD40972a','PD41043a','PD40982a','PD51929a','PD41002a','PD41002d','PD40996a','PD51972a','PD41002e','PD40985a','PD40990a','PD51993a','PD40969a','PD41939d','PD41906a','PD40965e','PD40967a','PD41910c','PD41923f','PD51948a')), y=count, fill=signature)) + geom_bar(stat="identity", colour="gray") + scale_fill_brewer(palette="Paired")

p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle=90))


#Plotting SNV count

tmb_data <- read.csv("data/Supplementary_Table_1.csv", sep=",", header=TRUE)
q <- ggplot(tmb_data, aes(x=factor(Tumor_Sample_Barcode, level=c('PD40965f','PD40974a','PD40976a','PD40978d','PD40980a','PD40997a','PD41004a','PD41017a','PD41020a','PD41023f','PD41025a','PD41026d','PD41029e','PD41030a','PD41032a','PD41033a','PD41039a','PD41039d','PD41043d','PD41044a','PD41045a','PD41046d','PD41895a','PD41900a','PD41901a','PD41909a','PD41910a','PD41915a','PD41915c','PD41916a','PD41921a','PD41923g','PD51928a','PD51928d','PD51939a','PD51940d','PD51969a','PD51978a','PD51979a','PD41046a','PD41038a','PD41046e','PD40967d','PD40956d','PD41932d','PD41920c','PD40971d','PD41928d','PD41920a','PD41896a','PD40971a','PD40961d','PD40961e','PD41035a','PD51969d','PD41001a','PD40961a','PD41035d','PD40978a','PD41913a','PD41021d','PD40987a','PD41913e','PD51952a','PD41036a','PD40987d','PD40986a','PD40957a','PD40970a','PD40994a','PD40986d','PD40964a','PD40966a','PD41020d','PD40952a','PD51930a','PD51932a','PD41907a','PD40962a','PD40973a','PD41027a','PD41930d','PD41927d','PD40980d','PD40972a','PD41011a','PD41043a','PD40982a','PD51929a','PD41002a','PD41002d','PD40989a','PD40983a','PD40996a','PD51972a','PD41000a','PD51951a','PD41002e','PD40985a','PD40983e','PD40969d','PD40983d','PD40990a','PD41912a','PD51993a','PD40969a','PD41939d','PD40963a','PD40965a','PD40966d','PD40968a','PD40984a','PD40988a','PD41905a','PD41906a','PD51982a','PD40965e','PD40967a','PD41910c','PD41923f','PD51948a')), y=snv)) + geom_bar(stat="identity", fill="light blue") 

q + theme(panel.background = element_rect(fill="white"), axis.text.x = element_text(angle=90))
dev.off()



```

References:

Sigprofiler Matrix Generator
https://github.com/AlexandrovLab/SigProfilerMatrixGenerator

Sigprofiler Extractor
https://github.com/AlexandrovLab/SigProfilerExtractor 

