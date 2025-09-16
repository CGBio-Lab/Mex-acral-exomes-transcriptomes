# Acral melanoma exomes' analysis
# Correlation analysis 
## Author: Patricia Basurto Lozada

## Correlation between copy number alterations and driver mutational status (Figure 2c)

To compare the burden of copy number alterations between groups of samples with different driver mutational status, we used scores generated with the web application CNApp. This score were generated using ASCAT output data where segment means (seg.mean) were calculated as log2(cn/ploidy) and default parameters in CNApp.

For the analyisis just one sample per patient was used, giving priority to primary samples when available.  

``` python 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Reading file with clinical data, copy number scores and mutation status of all samples
clindata = pd.read_csv("data/Supplementary_Table_1.csv", sep=",")
#Simplifying sample type
clindata["Sample_type"] = clindata["Sample_type"].replace({"LN_recurrence":"Recurrence", "Local_recurrence":"Recurrence", "Pulmonar_metastasis":"metastasis"})
#Filtering out samples with no copy number data 
cn_scores_filtered = clindata.dropna(subset=["GCS","FCS","BCS"])
#Filtering the data frame to get just one sample per patient (giving priotity to primary samples if available)

#Reading sample IDs for one sample per patient
per_patient = pd.read_csv("data/copy_number_alterations/ASCAT_per_patient_pass_samples.txt", header=None)
#Filtering through merge 
per_patient_scores = per_patient.merge(cn_scores_filtered, how="left", left_on=0, right_on="Sample")

#Plotting boxplots comparing GCS (Global Copy Number Scores) by mutation status
#Setting the aesthetics of the figure
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "NRAS":"#D291BC", "multihit":"#FFD758", "NF1":"#A6D472", "QWT":"lightgray", "KIT":"#7ec4cf"}
#Plotting the boxplot of GCS (global copy number scores) by driver mutational status
sns.boxplot(data=per_patient_scores, y="GCS", x="Mutation_status", palette=mypal_mut, showfliers = False, order=["BRAF","NRAS","NF1","KIT","QWT"])
sns.despine(offset=10, trim=False)
my_pal = {"metastasis": "indianred", "primary": "gray", "Recurrence":"blue", "Lesion_in_transit":"green", "LN_metastasis":"purple"}
#Plotting individual data points over the boxplot
sns.stripplot(x="Mutation_status", y="GCS", data=per_patient_scores, size=8, hue="Sample_type", linewidth=0, palette=my_pal,order=["BRAF","NRAS","NF1","KIT","QWT"])
```
Statistical testing

GCS scores by mut

``` python
#Separating GCS scores by mutational status
BRAF_GCS = per_patient_scores.loc[per_patient_scores['Mutation_status'] == "BRAF", 'GCS']
NRAS_GCS = per_patient_scores.loc[per_patient_scores['Mutation_status'] == "NRAS", 'GCS']
KIT_GCS = per_patient_scores.loc[per_patient_scores['Mutation_status'] == "KIT", 'GCS']
NF1_GCS = per_patient_scores.loc[per_patient_scores['Mutation_status'] == "NF1", 'GCS']
QWT_GCS = per_patient_scores.loc[per_patient_scores['Mutation_status'] == "QWT", 'GCS']

#Running Mann Whitney test
from scipy.stats import shapiro
from scipy.stats import mannwhitneyu

mannwhitneyu(BRAF_GCS, KIT_GCS)
mannwhitneyu(BRAF_GCS, NF1_GCS)
mannwhitneyu(BRAF_GCS, NRAS_GCS)
mannwhitneyu(BRAF_GCS, QWT_GCS)
mannwhitneyu(NRAS_GCS, KIT_GCS)
mannwhitneyu(NRAS_GCS, NF1_GCS)
mannwhitneyu(NRAS_GCS, QWT_GCS)
mannwhitneyu(KIT_GCS, NF1_GCS)
mannwhitneyu(KIT_GCS, QWT_GCS)
mannwhitneyu(NF1_GCS, QWT_GCS)
```

## Correlation between copy number alterations and anatomical site (Figure 2e)

We used the same approach to compare copy number alterations between tumors from differente anatomical sites. 


```python

#Setting the figure aesthetics
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette2 = {"#FEB7BB", "#94A6D8", "#94D8A6"}
sns.set_palette(my_palette2)
mypal_mut = {"foot":"#FEB7BB", "subungual":"#94A6D8", "hand":"#94D8A6"}
#Plotting the boxplot of GCS (global copy number scores) by driver mutational status
sns.boxplot(data=per_patient_scores, y="GCS", x="Primary_tumour_site", palette=mypal_mut, order=["foot","hand","subungual"])
my_pal = {"metastasis": "indianred", "primary": "gray", "Recurrence":"blue", "Lesion_in_transit":"green", "LN_metastasis":"purple"}
sns.stripplot(x="Primary_tumour_site", y="GCS", data=per_patient_scores, size=8, color=".3", linewidth=0, hue="Sample_type",palette=my_pal, order=["foot","hand","subungual"] )
sns.despine(offset=10, trim=False)
```
Statistical testing

```python
#Separating GCS scores by site
HAND_GCS = per_patient_scores.loc[per_patient_scores['Primary_tumour_site'] == "hand", 'GCS']
FOOT_GCS = per_patient_scores.loc[per_patient_scores['Primary_tumour_site'] == "foot", 'GCS']
SUBUNGUAL_GCS = per_patient_scores.loc[per_patient_scores['Primary_tumour_site'] == "subungual", 'GCS']

#Running Mann Whitney test
mannwhitneyu(HAND_GCS, FOOT_GCS)
mannwhitneyu(HAND_GCS, SUBUNGUAL_GCS)
mannwhitneyu(FOOT_GCS, SUBUNGUAL_GCS)

```


## Plotting GCS vs burden of TMB(snv+indels) (Figure 2d)

We compared de GCS (Global copy number scores) obtained from CNApp to the snv and indel count of the samples. 

``` python
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
mypal_mut = {"BRAF":"#FF9AA2", "NRAS":"#D291BC", "multihit":"#FFD758", "NF1":"#A6D472", "QWT":"lightgray", "KIT":"#7ec4cf"}
sns_plot2 = sns.scatterplot(x="TMB", y="GCS", data=per_patient_scores, hue="Mutation_status", palette=mypal_mut, s=250, style="Mutation_status", linewidth=1.5, edgecolor="black")

```

## Comparing TMB (snv + indel) by mutation status (Supplementary Figure 6)

``` python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#Reading data 
data = pd.read_csv("data/Supplementary_Table_1.csv", sep=",")
#Simplifying sample type
data["Sample_type"] = data["Sample_type"].replace({"LN_recurrence":"Recurrence", "Local_recurrence":"Recurrence", "Pulmonar_metastasis":"metastasis"})
#Setting up the aesthetics of the plot
sns.set(rc={'figure.figsize':(15,10)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "NRAS":"#D291BC", "multihit":"#FFD758", "NF1":"#A6D472", "QWT":"lightgray", "KIT":"#7ec4cf"}
#Plotting the boxplot
sns.boxplot(x="Mutation_status", y="TMB", data=data, palette=mypal_mut , showfliers = False, order=["BRAF","NRAS","NF1","KIT","multihit","QWT"]  )
sns.despine(offset=10, trim=False)
#Plotting individual data points
my_pal = {"metastasis": "indianred", "primary": "gray", "Recurrence":"blue", "Lesion_in_transit":"green", "LN_metastasis":"purple"}
sns.swarmplot(x="Mutation_status", y="TMB", data=data, size=8, hue="Sample_type", palette=my_pal, order=["BRAF","NRAS","NF1","KIT","multihit","QWT"])

```

## Comparing proportion of amerindian ancestry by driver mutational status (Supplementary Figure 3)

``` python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy

#Reading ancestry data
ancestry_data = pd.read_csv("data/Supplementary_Table_2.csv", sep=",")
#Reading clinical data for all samples
clinical_data = pd.read_csv("data/Supplementary_Table_1.csv", sep=",")
#Reading IDs for one sample per patient
per_patient_ID = pd.read_csv("data/One_sample_per_patient_ID.csv")

#Filtering data to just one sample per patient
clin_data_per_patient = per_patient_ID.merge(clinical_data, how="left", left_on="Tumor_Sample_Barcode", right_on="Sample")

#Removing last character of the ID to match the patients 
ancestry_data['Patient'] = ancestry_data['ID'].str[:-1]

#Mergin ancestry data to clinical data
per_patient_ancestry = clin_data_per_patient.merge(ancestry_data, how="left", on="Patient")

#Filtering data to leave only patients with genotyping data available
per_patient_ancestry_filtered = per_patient_ancestry.dropna(subset=["Q1 (AFR)","Q2 (AMR)","Q3 (SAS)","Q4 (EAS)","Q5 (EUR)"])

##Plotting the boxplot

#Setting the aesthetics of the figure
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "NRAS":"#D291BC", "multihit":"#FFD758", "NF1":"#A6D472", "QWT":"lightgray", "KIT":"#7ec4cf"}
#Plotting the boxplot of Amerindian ancestry (Q2 (AMR))
sns.boxplot(data=per_patient_ancestry_filtered, y="Q2 (AMR)", x="Mutation_status", palette=mypal_mut, showfliers = False, order=["BRAF","NRAS","NF1","KIT","multihit","QWT"])
sns.despine(offset=10, trim=False)
#Plotting individual data points over the boxplot
sns.stripplot(x="Mutation_status", y="Q2 (AMR)", data=per_patient_ancestry_filtered, color="gray", size=8, linewidth=0, order=["BRAF","NRAS","NF1","KIT","multihit","QWT"])

#Grouping data by mutation status

BRAF_patients = per_patient_ancestry_filtered[per_patient_ancestry_filtered["Mutation_status"]=="BRAF"]
NF1_patients = per_patient_ancestry_filtered[per_patient_ancestry_filtered["Mutation_status"]=="NF1"]
KIT_patients = per_patient_ancestry_filtered[per_patient_ancestry_filtered["Mutation_status"]=="KIT"]
NF1_patients = per_patient_ancestry_filtered[per_patient_ancestry_filtered["Mutation_status"]=="NF1"]
NRAS_patients = per_patient_ancestry_filtered[per_patient_ancestry_filtered["Mutation_status"]=="NRAS"]


#Statistical analysis

from scipy.stats import mannwhitneyu
mannwhitneyu(BRAF_patients["Q2 (AMR)"], KIT_patients["Q2 (AMR)"])
mannwhitneyu(BRAF_patients["Q2 (AMR)"], NRAS_patients["Q2 (AMR)"])
mannwhitneyu(NF1_patients["Q2 (AMR)"], KIT_patients["Q2 (AMR)"])
mannwhitneyu(NF1_patients["Q2 (AMR)"], BRAF_patients["Q2 (AMR)"])
mannwhitneyu(NF1_patients["Q2 (AMR)"], NRAS_patients["Q2 (AMR)"])
mannwhitneyu(KIT_patients["Q2 (AMR)"], NRAS_patients["Q2 (AMR)"])
```

