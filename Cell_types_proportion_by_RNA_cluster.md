# Acral melanoma exomes' analysis
# RNA deconvolution cell types by cluster 
## Author: Patricia Basurto Lozada

## Comparing presence of different cell types (Figure 4c, d, e, f and Supplementary Figure 13)

Using the outputs of the deconvolution analysis we compared the proprotion of different cell types present in tumors from different transcriptional clusters. 

``` python
import pandas as pd
#Reading file with clinical and RNA cluster data 
clindata = pd.read_csv("data/Supplementary_Table_1.csv", sep=",")
#Reading file with deconvolution data
deconv = pd.read_csv("data/Deconvolution_data/Sample_pass_deconvolution_data.csv", sep="\t")
deconv['DNA_Sample'] = deconv['sample'].str.replace('PR', 'PD')
#Merging clinical and cluster data with deconvolution data
deconv_clusters = deconv.merge(clindata, how="left", left_on="DNA_sample", right_on="Sample")


#Plotting boxplots of proportion of cell types by cluster
import seaborn as sns
import matplotlib.pyplot as plt


#Mitotic index by RNA cluster (Figure 4c)

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=deconv_clusters, y="Mitotic_index", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
sns.stripplot(x="RNA_cluster", y="Mitotic_index", data=deconv_clusters, size=7, linewidth=0, color="gray")


#Bcells cells by cluster (Figure 4d)

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=deconv_clusters, y="Bcells", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
sns.stripplot(x="RNA_cluster", y="Bcells", data=deconv_clusters, size=7, linewidth=0, color="gray")


#CD4+ T cells by cluster (Figure 4e)

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=deconv_clusters, y="CD4_Tcells", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
sns.stripplot(x="RNA_cluster", y="CD4_Tcells", data=deconv_clusters, size=7, linewidth=0, color="gray")

#CAFs by cluster (Figure 4f)

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=deconv_clusters, y="CAFs", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
sns.stripplot(x="RNA_cluster", y="CAFs", data=deconv_clusters, size=7, linewidth=0, color="gray")

#### Supplementary Figure 13 ##### 

#Endothelial cells by cluster 

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=deconv_clusters, y="Endothelial", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
sns.stripplot(x="RNA_cluster", y="Endothelial", data=deconv_clusters, size=7, linewidth=0, color="gray")

#CD8+ Tcells

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=deconv_clusters, y="CD8_Tcells", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
sns.stripplot(x="RNA_cluster", y="CD8_Tcells", data=deconv_clusters, size=7, linewidth=0, color="gray")

#Macrophages cells by cluster

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=deconv_clusters, y="Macrophages", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
sns.stripplot(x="RNA_cluster", y="Macrophages", data=deconv_clusters, size=7, linewidth=0, color="gray")

#otherCells cells by cluster

sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.set_palette(my_palette)
mypal_cluster = {"#A99AEC", "#EC99C2", "#F1955C"}
sns.boxplot(data=deconv_clusters, y="otherCells", x="RNA_cluster", palette=mypal_cluster, showfliers = False)
sns.set_palette(mypal_cluster)
sns.despine(offset=10, trim=False)
sns.stripplot(x="RNA_cluster", y="otherCells", data=deconv_clusters, size=7, linewidth=0, color="gray")

```

Statistical testing

```python

#Separating data by cluster 
cluster1_data= deconv_clusters.loc[deconv_clusters['cluster'] == 1]
cluster2_data= deconv_clusters.loc[deconv_clusters['cluster'] == 2]
cluster3_data= deconv_clusters.loc[deconv_clusters['cluster'] == 3]

#Comparing CAFs

cluster1_CAFs = cluster1_data["CAFs"]
cluster2_CAFs = cluster2_data["CAFs"]
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

#Comparing CD8+ Tcells

cluster1_CD8 = cluster1_data["CD8_Tcells"]
cluster2_CD8 = cluster2_data["CD8_Tcells"]
cluster3_CD8 = cluster3_data["CD8_Tcells"]

#Running Mann Whitney test

mannwhitneyu(cluster1_CD8, cluster2_CD8)
mannwhitneyu(cluster1_CD8, cluster3_CD8)
mannwhitneyu(cluster2_CD8, cluster3_CD8)

#Comparing B cells

cluster1_Bcells = cluster1_data["Bcells"]
cluster2_Bcells = cluster2_data["Bcells"]
cluster3_Bcells = cluster3_data["Bcells"]

#Running Mann Whitney test

mannwhitneyu(cluster1_Bcells, cluster3_Bcells)
mannwhitneyu(cluster1_Bcells, cluster2_Bcells)
mannwhitneyu(cluster3_Bcells, cluster2_Bcells)

#Comparing macrophages

cluster1_macrophages = cluster1_data["Macrophages"]
cluster2_macrophages = cluster2_data["Macrophages"]
cluster3_macrophages = cluster3_data["Macrophages"]

#Running Mann Whitney test

mannwhitneyu(cluster1_macrophages, cluster2_macrophages)
mannwhitneyu(cluster1_macrophages, cluster3_macrophages)
mannwhitneyu(cluster2_macrophages, cluster3_macrophages)

#Comparing other cells

cluster1_other = cluster1_data["otherCells"]
cluster2_other = cluster2_data["otherCells"]
cluster3_other = cluster3_data["otherCells"]

#Running Mann Whitney test

mannwhitneyu(cluster1_other, cluster2_other)
mannwhitneyu(cluster1_other, cluster3_other)
mannwhitneyu(cluster2_other, cluster3_other)

#Comparing mitotic index

cluster1_mitotic_index = cluster1_data["mitotic_index"]
cluster2_mitotic_index = cluster2_data["mitotic_index"]
cluster3_mitotic_index = cluster3_data["mitotic_index"]

#Running Mann Whitney test
mannwhitneyu(cluster1_mitotic_index, cluster3_mitotic_index)
mannwhitneyu(cluster1_mitotic_index, cluster2_mitotic_index)
mannwhitneyu(cluster3_mitotic_index, cluster2_mitotic_index)
```