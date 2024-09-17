# Acral melanoma exomes' analysis
# Acral/Cutaneous score by BRAF/NRAS mutational status
## Author: Patricia Basurto Lozada

## Comparing by BRAF/NRAS mutational status

We compared acral/cutaneous scores generated from transcriptomic data of 80 primary samples by BRAF and NRAS mutational status. 


```python
#Importing score data for 80 samples, one sample per patient
ratio_data_80 = pd.read_csv("80_samples_ratios_BRAFvsNRASvsWT.csv", sep=",")

#Generating the boxplot comparing BRAF mutated, NRAS mutated and BRAF/NRAS wildtype scores
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "WT":"lightgray", "NRAS":"#D291BC"}
sns.boxplot(data=ratio_data_80_NRAS, y="Ratio", x="Mutation_status", palette=mypal_mut, showfliers = False )
sns.despine(offset=10, trim=False)
#my_pal = {"metastasis": "indianred", "primary": "gray", "Recurrence":"blue", "Lesion_in_transit":"green"}
sns.stripplot(x="Mutation_status", y="Ratio", data=ratio_data_80_NRAS, color="gray", size=8, linewidth=0)
plt.yscale('log')
plt.ylim(10**-2, 10**1)

```

We also compared scores obtained from transcriptomic data of 63 samples from Newell et al (2020) study by BRAF mutational status

```python 
#Importing score data from Newell et al, 2020
n_data = pd.read_csv("Newell_data_ratios_bra_vs_brafwt.csv", sep=",")

#Generating the boxplot comparing BRAF mutated vs BRAF wildtype scores
sns.set(rc={'figure.figsize':(9,9)})
sns.set_style("white")
sns.set_context("talk")
my_palette = {"#FF9AA2", "#D291BC", "#FFD758", "#A6D472", "lightgray", "#7ec4cf"}
sns.set_palette(my_palette)
mypal_mut = {"BRAF":"#FF9AA2", "WT":"lightgray"}
sns.boxplot(data=n_data, y="Ratio", x="Mutation", palette=mypal_mut, showfliers = False , order=["WT","BRAF"])
sns.despine(offset=10, trim=False)
#my_pal = {"metastasis": "indianred", "primary": "gray", "Recurrence":"blue", "Lesion_in_transit":"green"}
sns.stripplot(x="Mutation", y="Ratio", data=n_data, color="gray", size=8, linewidth=0, order=["WT","BRAF"])
sns.despine(offset=10, trim=False)
plt.yscale('log')
plt.ylim(10**-2, 10**1)


```