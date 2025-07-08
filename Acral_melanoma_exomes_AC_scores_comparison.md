# Acral melanoma exomes' analysis
# Acral/Cutaneous score by BRAF mutational status
## Author: Patricia Basurto Lozada


## Visualization of the comparison of acral:cutaneous scores by BRAF mutational status (Figure 3c)

We visualized the comparison acral/cutaneous scores generated from transcriptomic data of 77 primary samples by BRAF mutational status. 


```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import pyplot

#Importing score data for 77 samples, one sample per patient
ratios = pd.read_csv("/data/AC_scores/ratios_ac_score.csv", sep=","
#Changing PR for PD to match DNA ids
ratios['Sample'] = ratios['Sample'].str.replace('PR', 'PD')
#Reading clincal data 
clindata = pd.read_csv("/data/Supplementary_Table_1", sep=",")
#Mergin ratio data with clinical data
ratios_clin = ratios.merge(clindata, how="left", on="Sample")
#Generating a new column to annotate if sample is BRAF mutated or BRAF wt 
ratios_clin['AC_BRAF_status'] = ratios_clin['Mutation_status'].apply(lambda x: 'BRAF' if x == 'BRAF' else 'WT')

#Generating the boxplot comparing BRAF mutated and BRAF wildtype scores
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "WT":"lightgray", "NRAS":"#D291BC"}
sns.boxplot(data=ratios_clin, y="Value", x="AC_BRAF_status", palette=mypal_mut, showfliers = False, order=["BRAF","WT"] )
sns.despine(offset=10, trim=False)
sns.stripplot(x="AC_BRAF_status", y="Value", data=ratios_clin, color="gray", size=8, linewidth=0, order=["BRAF","WT"])
plt.yscale('log')
plt.ylim(10**-2, 10**1)

#Separating samples by mutation status

BRAF_data = ratios_clin.loc[ratios_clin['AC_BRAF_status'] == "BRAF"]
BRAF_data = BRAF_data["Value"]

WT_data = ratios_clin.loc[ratios_clin['AC_BRAF_status'] == "WT"]
WT_data = WT_data["Value"]

#Statistical analysis

from scipy.stats import mannwhitneyu

mannwhitneyu(BRAF_data, WT_data, alternative="less")
```

We also compared scores obtained from transcriptomic data of 63 samples from Newell et al (2020) study by BRAF mutational status. Samples were considered as BRAF mutated if they had an activating mutation but were considered as WT if the sample only had structural variants (SV) or amplifications 

```python 
#Importing score data from Newell et al, 2020
newell_data = pd.read_csv("data/AC_scores/ratios_newell_all.csv", sep=",")

#Changing samples that had SV or amplifications in BRAF to WT
newell_data["only_BRAFmut_status"] = newell_data["BRAF_status"]
newell_data.loc[(newell_data['ID'] == 'MELA_0278') & (newell_data['only_BRAFmut_status'] == 'BRAF'), 'only_BRAFmut_status'] = 'WT'
newell_data.loc[(newell_data['ID'] == 'MELA_0268') & (newell_data['only_BRAFmut_status'] == 'BRAF'), 'only_BRAFmut_status'] = 'WT'
newell_data.loc[(newell_data['ID'] == 'MELA_0270') & (newell_data['only_BRAFmut_status'] == 'BRAF'), 'only_BRAFmut_status'] = 'WT'


#Generating the boxplot comparing BRAF mutated vs BRAF wildtype scores
ssns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "WT":"lightgray", "NRAS":"#D291BC"}
sns.boxplot(data=newell_data, y="values", x="only_BRAFmut_status", palette=mypal_mut, showfliers = False, order=["BRAF","WT"] )
sns.despine(offset=10, trim=False)
sns.stripplot(x="only_BRAFmut_status", y="values", data=newell_data, color="gray", size=8, linewidth=0, order=["BRAF","WT"])
plt.yscale('log')
plt.ylim(10**-2, 10**1)

#Separating samples by mutation status

BRAFmut_newell= newell_data.loc[newell_data['only_BRAFmut_status'] == "BRAF"]
BRAFmut_newell = BRAFmut_newell["values"]

WTmut_newell = newell_data.loc[newell_data['only_BRAFmut_status'] == "WT"]
WTmut_newell = WTmut_newell["values"]

#Statistical analysis

from scipy.stats import mannwhitneyu

mannwhitneyu(BRAFmut_newell, WTmut_newell, alternative="less")
```

# McNeal melanocytes cutaneous score (Figure 3c and d)

We compared the cutaneous classifier genes normalized expression and cutaneous score of melanocytes with induced BRAFV600E in PMA conditions and ET1 conditions.  

Plotting of the comparison:


```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import pyplot

# Reading the CM scores 
McNeal_cm_scores = pd.read_csv("data/AC_scores/McNeal_CM_score.csv", sep=",")

#Plotting CM scores as individual replicates
sns.set(rc={'figure.figsize':(10,10)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#F4AE9F", "#C2C1C1", "#F1CDC5", "#989595"}
sns.set_palette(my_palette)
mypal_mut = {"PMA BRAFV600E":"#F4AE9F", "PMA":"#C2C1C1", "ET1 BRAFV600E":"#F1CDC5", "ET1":"#989595" }
sns.stripplot(x="Category", y="Score", data=McNeal_cm_scores, hue="Category", size=20, jitter=1, palette=mypal_mut, linewidth=0, edgecolor="black")
sns.boxplot(showmeans=False,
            meanline=False,
            medianprops={'visible': True},
            whiskerprops={'visible': False},
            zorder=10,
            x="Category",
            y="Score",
            data=McNeal_cm_scores,
            showfliers=False,
            showbox=False,
            showcaps=False)
plt.yscale('log')
plt.ylim(10**10, 10**11)
sns.despine(offset=10, trim=False)


#Gene expression

#Reading gene expression data
genes_data = pd.read_csv("data/AC_scores/McNeal_genes_data.csv", sep=",")
#Formatting dataframe to have genes as a separate column
genes_data_melted = genes_data.melt(id_vars='Category', var_name='Gene', value_name='Expression')

#Plotting data
sns.set(rc={'figure.figsize':(15,10)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#F4AE9F", "#C2C1C1", "#F1CDC5", "#989595"}
sns.set_palette(my_palette)
mypal_mut = {"PMA BRAFV600E":"#F4AE9F", "PMA":"#C2C1C1", "ET1 BRAF V600E":"#F1CDC5", "ET1":"#989595"}
#sns.boxplot(data=McNeal_data, y="Score", x="Category", palette=mypal_mut, showfliers = False,  )
sns.stripplot(x="Gene", y="Expression", data=genes_data_melted, hue="Category", size=12, linewidth=0, jitter=0.2, palette=mypal_mut)
sns.boxplot(showmeans=False,
            meanline=False,
            medianprops={'visible': True},
            whiskerprops={'visible': False},
            zorder=10,
            x="Gene",
            y="Expression",
            data=genes_data_melted,
            showfliers=False,
            showbox=False,
            showcaps=False)
sns.despine(offset=10, trim=False)
plt.ylim(0, 16)

```

# Plotting TCGA data  (Supplementary Figure 9)

We compared the A:C ratios (ratios were calculated with the same strategy as the previous datasets) of BRAF mutated and BRAF WT samples from the TCGA

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import pyplot
from scipy.stats import mannwhitneyu

# Reading ratios for BRAF mutated samples
TCGA_ratios = pd.read_csv("data/AC_scores/TCGA_ratios.csv", sep=",")
#Separating WT and BRAF samples
WT_ratios = TCGA_ratios.loc[TCGA_ratios["Mut_status"] == "WT"]
BRAF_ratios = TCGA_ratios.loc[TCGA_ratios["Mut_status"] == "BRAF"]


#Setting up the aesthetics of the plot
sns.set(rc={'figure.figsize':(15,10)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "NRAS":"#D291BC", "multihit":"#FFD758", "NF1":"#A6D472", "WT":"lightgray", "KIT":"#7ec4cf"}
# Plotting the boxplots
sns.boxplot(data=TCGA_ratios, y="score", x="Mut_status", palette=mypal_mut, showfliers = False, order=["BRAF","WT"] )
# Plotting individual data points 
sns.stripplot(x="Mut_status", y="score", data=TCGA_ratios, color="gray", size=6, linewidth=0, order=["BRAF","WT"], jitter=True)
plt.yscale('log')
plt.ylim(10**-20, 10**10)
sns.despine(offset=10, trim=False)

# Statistical analysis

# Test
mannwhitneyu(BRAF_ratios["score"], WT_ratios["score"])

```

