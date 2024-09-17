# Acral melanoma exomes' analysis
# RNA deconvolution cell types by cluster 
## Author: Patricia Basurto Lozada

## Comparing presence of different cell types

Using the outputs of the deconvolution analysis we compared the proprotion of different cell types present in tumors from different transcriptional clusters. 

``` python
import pandas as pd
#Reading file with samples and cluster data
all_clin_data = pd.read_csv("Supplementary_table_1.csv", sep=",")
#Reading file with deconvolution data
deconvolution = pd.read_csv("deconvolution_60samples_primaries_epic.csv", sep=",")

#Merging clinical and cluster data with deconvolution data

cluster_clin_data = pd.merge(all_clin_data, deconvolution, how="right")
cluster_clin_data2 = cluster_clin_data.replace("<1", 0.9)

#Plotting 
import seaborn as sns
import matplotlib.pyplot as plt

#CD4+ T cells by cluster

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=cluster_deconv, y="CD4_Tcells", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
#my_pal = {"metastasis": "indianred", "primary": "gray", "unknown":"blue"}
sns.stripplot(x="RNA_cluster", y="CD4_Tcells", data=cluster_deconv, size=7, linewidth=0, color="gray")


#Endothelial cells by cluster 

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=cluster_deconv, y="Endothelial", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
#my_pal = {"metastasis": "indianred", "primary": "gray", "unknown":"blue"}
sns.stripplot(x="RNA_cluster", y="Endothelial", data=cluster_deconv, size=7, linewidth=0, color="gray")


#CAFs by cluster

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=cluster_deconv, y="CAFs", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
#my_pal = {"metastasis": "indianred", "primary": "gray", "unknown":"blue"}
sns.stripplot(x="RNA_cluster", y="CAFs", data=cluster_deconv, size=7, linewidth=0, color="gray")


#Mitotic index by RNA cluster 

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=cluster_clin_data2, y="Mitotic_index_x", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
#my_pal = {"metastasis": "indianred", "primary": "gray", "unknown":"blue"}
sns.stripplot(x="RNA_cluster", y="Mitotic_index_x", data=cluster_clin_data2, size=7, linewidth=0, color="gray")


```