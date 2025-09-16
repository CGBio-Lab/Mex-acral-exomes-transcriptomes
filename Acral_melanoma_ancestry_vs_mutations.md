# Acral melanoma - Ancestry vs Mutations
## Author: Irving Simonin-Wilmer


``` R
library(dplyr)
library(readxl)
library(stringr)
library(readr)
library(tibble)
library(ggplot2)
library(tidyr)
library(modelsummary)
library(forestplot)
library(readxl)
library(boot)
library(pROC)

# Sample IDs, one for each patient

sample_list <- c("PD40952a","PD40956d","PD40957a","PD40961a","PD40962a","PD40963a","PD40964a",
                "PD40965a","PD40966a","PD40967a","PD40968a","PD40969a","PD40970a","PD40971a",
                "PD40972a","PD40973a","PD40974a","PD40976a","PD40978a","PD40980a","PD40982a",
                "PD40983a","PD40984a","PD40985a","PD40986a","PD40987a","PD40988a","PD40989a",
                "PD40990a","PD40994a","PD40996a","PD40997a","PD41000a","PD41001a","PD41002a",
                "PD41004a","PD41011a","PD41017a","PD41020a","PD41021d","PD41023f","PD41025a",
                "PD41026d","PD41027a","PD41029e","PD41030a","PD41032a","PD41033a","PD41035a",
                "PD41036a","PD41038a","PD41039d","PD41043a","PD41044a","PD41045a","PD41046a",
                "PD41895a","PD41896a","PD41900a","PD41901a","PD41905a","PD41906a","PD41907a",
                "PD41909a","PD41910a","PD41912a","PD41913a","PD41915c","PD41916a","PD41920a",
                "PD41921a","PD41923g","PD41927d","PD41928d","PD41930d","PD41932d","PD41939d",
                "PD51928a","PD51929a","PD51930a","PD51932a","PD51939a","PD51940d","PD51948a",
                "PD51951a","PD51952a","PD51969a","PD51972a","PD51978a","PD51979a","PD51982a",
                "PD51993a")

## Covariate Files ##

## Ancestry from ADMIXTURE run with 𝑘=5
ancestry <- read_excel('data/Supplementary_Table_2.xlsx', skip=1)  %>% 
    rename_with(~ str_replace(., "\\s.*", "")) %>%
    mutate(ID=str_sub(Sample,1,-2))


## The SNVs and indel counts by tumor:
snv_indel_data <- read_csv('data/Supplementary_Table_1.csv') %>% 
  mutate(ID=str_sub(Sample,1,-2))

##Clinical data for patients is merged with the counts and only the relevant tumor samples are selected using the sample_list.

data <- snv_indel_data %>%  
  filter(Sample %in% sample_list) %>% 
  select(Sample, ID, Sex, Age, Socioeconomic_status, TMB, Mutation_status, Ulceration_status, Tumour_stage) %>% 
  mutate(values=1) %>%
  pivot_wider(names_from = Mutation_status,values_from = values, values_fill=0) %>% 
  inner_join(ancestry %>% filter(SUPERPOP=="-") %>% select(ID,Q2,Q5)) %>%
  arrange(Q2)

data %>% filter(BRAF==1)

## EUR related ancestry

##Formula to be used with each gene. The Q5 ADMIXTURE cluster is related to the EUR ancestry in the 1000 Genomes Project.
dependent_formula <- "~ Q5 + Age + Sex + TMB"

##Remove rows with missing covariates.
data_complete <- data %>% select(ID, Q5, Age, Sex, TMB ,NRAS ,KIT ,NF1 ,BRAF) %>% drop_na()

data_complete %>% filter(BRAF==1) %>% select(ID)

data_complete <- data_complete %>%
  mutate(Q5 = as.numeric(Q5))

#Logistic regression model for each gene

muts <- c("KIT","BRAF","NRAS","NF1")
results_list <- lapply(1:4, 
    function(x)
    as.data.frame(confint.default(glm(formula(paste0(muts[x],dependent_formula)),
             data=data_complete, family = "binomial"),"Q5")) %>% 
  rename(lower=`2.5 %`,upper=`97.5 %`) %>% 
  mutate(estimate=glm(formula(paste0(muts[x],dependent_formula)),
     data=data_complete, family = "binomial")$coefficients[2],
     p_value=coef(summary(glm(formula(paste0(muts[x],dependent_formula)),
     data=data_complete, family = "binomial")))[2,4],
     mutation=muts[x])
       )  
results_df <- do.call(rbind,results_list)
results_df

        lower    upper  estimate    p_value mutation
Q5  -6.2873964 2.274221 -2.006588 0.35824511      KIT
Q51  0.6639491 9.440615  5.052282 0.02403908     BRAF
Q52 -3.4425955 7.007578  1.782491 0.50373530     NRAS
Q53 -2.5576912 7.562502  2.502406 0.33240804      NF1

## Plotting 

results_df %>% 
  ggplot(aes(y=mutation)) + 
  geom_point(aes(x=estimate), shape=15, size=3) +
  geom_linerange(aes(xmin=lower, xmax=upper)) +
  xlab("Log Odds - Estimate and confidence intervals") +
  geom_text(aes(x=0,y=mutation,label=paste0("p=",round(p_value,2))),nudge_y=0.1)
```