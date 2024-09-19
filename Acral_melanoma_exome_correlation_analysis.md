# Acral melanoma exomes' analysis
# Correlation analysis 
## Author: Patricia Basurto Lozada

## Comparing age at diagnosis between samples with different driver mutational status 

We compared age at diagnosis with driver mutational status using one sample per patient (prioritizing primaries when available).

``` python 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#Reading data using one sample per patient (giving priority to primaries when available)
#One_per_patient_supplementary_table_1.csv was generated filtering Supplementary table 1 to only have one sample per patient prioritizing primaries. 

data = pd.read_csv("One_per_patient_supplementary_table_1.csv", sep=",")

#Setting up the aesthetics of the plot
sns.set(rc={'figure.figsize':(8,20)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "NRAS":"#D291BC", "multihit":"#FFD758", "NF1":"#A6D472", "Other":"lightgray", "KIT":"#7ec4cf"}
#Plotting the boxplot
sns.boxplot(y="Mutation_status", x="Age", data=data, palette=mypal_mut, order=["NRAS","NF1","KIT","BRAF","multihit","Other"])
sns.despine(offset=10, trim=False)
#Plotting individual data points
sns.swarmplot(y="Mutation_status", x="Age", data=data, size=8, color="GRAY", order=["NRAS","NF1","KIT","BRAF","multihit","Other"])
```


## Correlation between copy number alterations and driver mutational status

To compare the burden of copy number alterations between groups of samples with different driver mutational status, we used scores generated with the web application CNApp. This score were generated using sequenza output data as input and default parameters in CNApp.

For the analyisis just one sample per patient was used, giving priority to primary samples when available.  

``` python 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Reading table of samples with copy number data (one sample per patient prioritizing)
#CN_data_one_per_patient_supplementary_table_1.csv was generated filtering Supplementary table 1 to only have samples with copy number data and only one sample per patient prioritizing primaries. 
per_patient_cn_data= pd.read_csv("CN_data_one_per_patient_suppplementary_table_1.csv", sep=",")

#Setting the aesthetics of the figure
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "NRAS":"#D291BC", "multihit":"#FFD758", "NF1":"#A6D472", "Other":"lightgray", "KIT":"#7ec4cf"}
#Plotting the boxplot of GCS (global copy number scores) by driver mutational status
sns.boxplot(data=per_patient_cn_data, y="GCS", x="Mutation_status", palette=mypal_mut, showfliers = False, order=["NRAS","BRAF","NF1","KIT","multihit","Other"])
sns.despine(offset=10, trim=False)
my_pal = {"metastasis": "indianred", "primary": "gray", "Recurrence":"blue", "Lesion_in_transit":"green"}
#Plotting individual data points over the boxplot
sns.stripplot(x="Mutation_status", y="GCS", data=per_patient_cn_data, size=8, hue="Tumor_type", linewidth=0, palette=my_pal,order=["NRAS","BRAF","NF1","KIT","multihit","Other"])
```
Statistical testing

GCS scores by mut

``` python
#Separating GCS scores by mutational status
BRAF_GCS = per_patient_cn_data.loc[per_patient_cn_data['Mutation_status'] == "BRAF", 'GCS']
NRAS_GCS = per_patient_cn_data.loc[per_patient_cn_data['Mutation_status'] == "NRAS", 'GCS']
KIT_GCS = per_patient_cn_data.loc[per_patient_cn_data['Mutation_status'] == "KIT", 'GCS']
NF1_GCS = per_patient_cn_data.loc[per_patient_cn_data['Mutation_status'] == "NF1", 'GCS']
multihit_GCS = per_patient_cn_data.loc[per_patient_cn_data['Mutation_status'] == "multihit", 'GCS']
Other_GCS = per_patient_cn_data.loc[per_patient_cn_data['Mutation_status'] == "Other", 'GCS']

#Running Mann Whitney test
from scipy.stats import shapiro
from scipy.stats import mannwhitneyu

mannwhitneyu(BRAF_GCS, KIT_GCS)
mannwhitneyu(BRAF_GCS, NF1_GCS)
mannwhitneyu(BRAF_GCS, NRAS_GCS)
mannwhitneyu(BRAF_GCS, multihit_GCS)
mannwhitneyu(BRAF_GCS, Other_GCS)
mannwhitneyu(NRAS_GCS, KIT_GCS)
mannwhitneyu(NRAS_GCS, NF1_GCS)
mannwhitneyu(NRAS_GCS, multihit_GCS)
mannwhitneyu(NRAS_GCS, Other_GCS)
mannwhitneyu(KIT_GCS, NF1_GCS)
mannwhitneyu(KIT_GCS, multihit_GCS)
mannwhitneyu(KIT_GCS, Other_GCS)
mannwhitneyu(NF1_GCS, multihit_GCS)
mannwhitneyu(NF1_GCS, Other_GCS)
mannwhitneyu(multihit_GCS, Other_GCS)
```

## Correlation between copy number alterations and anatomical site

We used the same approach to compare copy number alterations between tumors from differente anatomical sites. 


``` python

#Setting the figure aesthetics
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette2 = {"#FEB7BB", "#94A6D8", "#94D8A6"}
sns.set_palette(my_palette2)
mypal_mut = {"FEET":"#FEB7BB", "SUBUNGUAL":"#94A6D8", "HAND":"#94D8A6"}
#Plotting the boxplot of GCS (global copy number scores) by driver mutational status
sns.boxplot(data=per_patient_cn_data, y="GCS", x="Site", palette=mypal_mut, order=["FEET","HAND","SUBUNGUAL"])
my_pal = {"metastasis": "indianred", "primary": "gray", "Recurrence":"blue", "Lesion_in_transit":"green"}
sns.stripplot(x="Site", y="GCS", data=per_patient_cn_data, size=8, color=".3", linewidth=0, hue="Tumor_type",palette=my_pal, order=["FEET","HAND","SUBUNGUAL"] )
sns.despine(offset=10, trim=False)
```
Statistical testing

```python
#Separating GCS scores by site
HAND_GCS = per_patient_cn_data.loc[per_patient_cn_data['Site'] == "HAND", 'GCS']
FEET_GCS = per_patient_cn_data.loc[per_patient_cn_data['Site'] == "FEET", 'GCS']
SUBUNGUAL_GCS = per_patient_cn_data.loc[per_patient_cn_data['Site'] == "SUBUNGUAL", 'GCS']

#Running Mann Whitney test
mannwhitneyu(HAND_GCS, FEET_GCS)
mannwhitneyu(HAND_GCS, SUBUNGUAL_GCS)
mannwhitneyu(FEET_GCS, SUBUNGUAL_GCS)

```

## Plotting GCS vs burden of snv+indels

We compared de GCS (Global copy number scores) obtained from CNApp to tumor mutational burden of the samples. 

``` python
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
mypal_mut = {"BRAF":"#FF9AA2", "NRAS":"#D291BC", "multihit":"#FFD758", "NF1":"#A6D472", "Other":"lightgray", "KIT":"#7ec4cf"}
sns_plot2 = sns.scatterplot(x="snp+indels", y="GCS", data=per_patient_cn_data, hue="Mutation_status", palette=mypal_mut, s=250, style="Mutation_status")

```
