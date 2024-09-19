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

cluster_deconv = pd.merge(all_clin_data, deconvolution, how="right")
cluster_deconv["RNA_cluster"] = cluster_clin_data["RNA_cluster"].astype("category")
cluster_deconv2 = cluster_clin_data.replace("<1", 0.9)

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
sns.boxplot(data=cluster_deconv, y="Mitotic_index", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
#my_pal = {"metastasis": "indianred", "primary": "gray", "unknown":"blue"}
sns.stripplot(x="RNA_cluster", y="Mitotic_index", data=cluster_deconv2, size=7, linewidth=0, color="gray")


```

Statistical testing

```python
#Comparing CAFs

#Separating data by cluster 
cluster1_data= cluster_deconv.loc[cluster_deconv['cluster'] == 1]
cluster1_CAFs = cluster1_data["CAFs"]

cluster2_data= cluster_deconv.loc[cluster_deconv['cluster'] == 2]
cluster2_CAFs = cluster2_data["CAFs"]

cluster3_data= cluster_deconv.loc[cluster_deconv['cluster'] == 3]
cluster3_CAFs = cluster3_data["CAFs"]

#Running Mann Whitney test
from scipy.stats import mannwhitneyu

mannwhitneyu(cluster1_CAFs, cluster2_CAFs)
mannwhitneyu(cluster2_CAFs, cluster3_CAFs)
mannwhitneyu(cluster1_CAFs, cluster3_CAFs)

#Comparing endothelial cells

cluster1_Endothelial = cluster1_data["Endothelial"]
cluster2_Endothelial = cluster2_data["Endothelial"]
cluster3_Endothelial = cluster3_data["Endothelial"]

#Running Mann Whitney test

mannwhitneyu(cluster1_Endothelial, cluster3_Endothelial)
mannwhitneyu(cluster1_Endothelial, cluster2_Endothelial)
mannwhitneyu(cluster3_Endothelial, cluster2_Endothelial)

#Comparing CD4+

cluster1_CD4 = cluster1_data["CD4_Tcells"]
cluster2_CD4 = cluster2_data["CD4_Tcells"]
cluster3_CD4 = cluster3_data["CD4_Tcells"]


#Running Mann Whitney test

mannwhitneyu(cluster1_CD4, cluster3_CD4)
mannwhitneyu(cluster1_CD4, cluster2_CD4)
mannwhitneyu(cluster3_CD4, cluster2_CD4)

#Comparing mitotic index

#Separating data by cluster 
cluster1_data= cluster_deconv2.loc[cluster_deconv2['RNA_cluster'] == 1]
cluster1_mitotic_index = cluster1_data["Mitotic_index"]

cluster2_data= cluster_deconv2.loc[cluster_deconv2['RNA_cluster'] == 2]
cluster2_mitotic_index = cluster2_data["Mitotic_index"]

cluster3_data= cluster_deconv2.loc[cluster_deconv2['RNA_cluster'] == 3]
cluster3_mitotic_index = cluster3_data["Mitotic_index"]

#Running Mann Whitney test
mannwhitneyu(cluster1_mitotic_index, cluster3_mitotic_index)
mannwhitneyu(cluster1_mitotic_index, cluster2_mitotic_index)
mannwhitneyu(cluster3_mitotic_index, cluster2_mitotic_index)
```
