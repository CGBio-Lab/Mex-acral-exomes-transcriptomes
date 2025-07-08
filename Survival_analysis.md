# Survival analysis 
## Author: Timothy Bishop

Overall and recurrence free survival analysis were done using STATA. Supplementary Table 1 was used as input after filtering only patients with primary samples (85 patients). Further filtering was done for each analysis if necessary for patiens with missing data used in the analysis (e.g. ancestry, RNAseq cluster data)

## Recurrence

Multivariate logistic regression to test the association of any driver mutation on recurrnce-free survival (Supplementary Table 22) 

``` 
keep if tumor_type == "primary" 
duplicates drop patient_file, force
generate dxdate = date(date_of_diagnosis,"DMY")
generate condate = date(date_of_recruitment,"MD20Y")
drop if condate == .
generate dxyear = year(dxdate)
generate dthdate = date(date_of_death,"DMY")
generate lndate = date(date_of_last_note,"DMY")
generate recdate = date(recurrence_date,"DMY")
format dxdate dthdate lndate recdate condate %td 
list patient_file dln lndate if lndate != dln
drop dln
list patient_file dxdate recdate condate dthdate lndate
generate death = 0
replace death = 1 if dthdate != .
generate sex = 0 if gender == "female"
replace sex = 1 if gender == "male"
replace mutation = lower(mutation)
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" {
generate i`gn' = .
replace i`gn' = 1 if mutation == "`gn'"
replace i`gn' = 0 if mutation == "other"
tab i`gn' mutation
}
generate imut = 1
replace imut = 0 if mutation_status == "QWT"
generate t = recdate-dxdate
replace t = lndate-dxdate if recdate == .
generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
replace cstage = max(cstage,2)
destring cluster, replace ignore("NA")
Note: imut is any driver mutation
tab imut recurrence, all row

+----------------+
| Key            |
|----------------|
| frequency.     |
| row percentage |
+----------------+
           |     recurrence
imut       | 0             1      |  Total
-----------+----------------------+----------
0          | 27            16     | 43
           | 62.79         37.21  | 100.00
-----------+----------------------+----------
1          | 14            28     | 42
           | 33.33         66.67  | 100.00
-----------+----------------------+----------
Total      | 41            44     | 85
           | 48.24         51.76  | 100.00
Pearson chi2(1) = 7.3839 Pr = 0.007
likelihood-ratio chi2(1) = 7.4967 Pr = 0.006
Cramér’s V = 0.2947
gamma = 0.5429 ASE = 0.160
Kendall’s tau-b = 0.2947 ASE = 0.104

generate dxcat = .
replace dxcat = 0 if dxyear < 2014.5
replace dxcat = 1 if dxyear > 2016.5

logistic recurrence imut sex age cstage dxdate Q2AMR

Logistic regression                                         Number of obs = 73
                                                            LR chi2(6) = 18.11
                                                            Prob > chi2 = 0.0060
Log likelihood = -41.484649                                 Pseudo R2 = 0.1791
------------------------------------------------------------------------------
recurrence   | Odds      Ratio    Std. Err.   z    P>|z|    [95% Conf. Interval]
-------------+----------------------------------------------------------------
imut         | 5.310556  3.325634   2.67.  0.008  1.556293    18.12127
sex          | .9844121  .5578017  -0.03   0.978  .3242328    2.988801
age          | 1.016565  .0244662.  0.68   0.495  .9697256    1.065667
cstage       | 3.545433  1.973259   2.27   0.023  1.191027    10.554
dxdate       | 1.000076  .000339    0.23   0.822  .9994121    1.000741
Q2AMR        | 21.06773  33.44283   1.92   0.055  .9384777    472.9461
_cons        | .0001449  .0010369  -1.24   0.217  1.18e-10    178.0292
------------------------------------------------------------------------------
```

Multivariate logistic regression to test the association of cluster assignmnet on recurrence-free survival (Supplementary Table 23)

``` 
tab cluster recurrence , exact row

+----------------+
| Key            |
|----------------|
| frequency      |
| row percentage |
+----------------+
Enumerating sample-space combinations:
stage 3: enumerations = 1
stage 2: enumerations = 8
stage 1: enumerations = 0
RNA_cluste |      recurrence
r          |    0          1      | Total
-----------+----------------------+----------
1          |    9          5      | 14
           |    64.29      35.71  | 100.00
-----------+----------------------+----------
2          |    3          13     | 16
           |    18.75      81.25  | 100.00
-----------+----------------------+----------
3          |    6          8      | 14
           |    42.86      57.14  | 100.00
-----------+----------------------+----------
Total      |    18         26     | 44
           |    40.91      59.09  | 100.00
Fisher’s exact = 0.039

xi: logistic recurrence i.cluster sex age cstage dxdate

i.cluster       _Icluster_1-3           (naturally coded; _Icluster_1 omitted)
Logistic regression                         Number of obs = 44
                                            LR chi2(6) = 9.30
                                            Prob > chi2 = 0.1573
Log likelihood = -25.116653                 Pseudo R2 = 0.1562
------------------------------------------------------------------------------
recurrence   | Odds Ratio     Std. Err.    z     P>|z|   [.95% Conf. Interval]
-------------+----------------------------------------------------------------
_Icluster_2  | 6.683601       6.598269     1.92  0.054    .9653317   46.2748
_Icluster_3  | 2.373717       2.045131     1.00  0.316    .4385876   12.84699
        sex  | 1.24812        .9355723     0.30  0.767    .2872155   5.42381
        age  | .9963314       .0305672    -0.12  0.905    .9381865   1.05808
     cstage  | 2.591876       1.867251     1.32  0.186    .6315211   10.63752
     dxdate  | 1.000321       .0006131     0.52  0.600    .9991203   1.001524
      _cons  | .0000775       .000994     -0.74  0.461    9.22e-16   6502608
------------------------------------------------------------------------------

``` 

# Overall survival

Log-rank test for homogeneity, any mutation on overall survival and Log-rank test of homogeneity, mutation classification on overall survival
(Supplementary Figure 24 and 25)

``` 
generate t = dthdate
replace t = lndate if dthdate == .
egen tag = tag(Sample)
tab tag
drop tag
generate Mutation = upper(mutation)
replace Mutation = " WT" if mutation == "QWT"

tab Mutation death, row all
+----------------+
| Key            |
|----------------|
| frequency      |
| row percentage |
+----------------+
           |         death
Mutation   |   0             1    | Total
-----------+----------------------+----------
       WT  |   40            3    | 43
           |   93.02         6.98 | 100.00
-----------+----------------------+----------
      BRAF |   7             4    | 11
           |   63.64        36.36 | 100.00
-----------+----------------------+----------
      KIT  |   10            1    | 11
           |   90.91         9.09 | 100.00
-----------+----------------------+----------
  MULTIHIT |   0             1    | 1
           | 0.00          100.00 | 100.00
-----------+----------------------+----------
      NF1  |   5             2    | 7
           |   71.43        28.57 | 100.00
-----------+----------------------+----------
      NRAS |   9             3    | 12
           |   75.00        25.00 | 100.00
-----------+----------------------+----------
     Total |   71           14    | 85
           | 83.53          16.47 | 100.00
Pearson chi2(5) = 12.8676 Pr = 0.025
likelihood-ratio chi2(5) = 11.3010 Pr = 0.046
Cramér’s V = 0.3891
gamma = 0.4118 ASE = 0.159
Kendall’s tau-b = 0.2100 ASE = 0.093

tab imut death, all exact row
+----------------+
| Key            |
|----------------|
| frequency      |
| row percentage |
+----------------+
           |         death
      imut |    0            1    | Total
-----------+----------------------+----------
         0 |   40            3    | 43
           |   93.02         6.98 | 100.00
-----------+----------------------+----------
         1 |   31            11   | 42
           |   73.81         26.19| 100.00
-----------+----------------------+----------
     Total |   71            14   | 85
           |   83.53         16.47| 100.00
Pearson chi2(1) = 5.7013 Pr = 0.017
likelihood-ratio chi2(1) = 5.9920 Pr = 0.014
Cramér’s V = 0.2590
gamma = 0.6510 ASE = 0.200
Kendall’s tau-b = 0.2590 ASE = 0.096
Fisher’s exact = 0.021
1-sided Fisher’s exact = 0.017

Logrank test of homogeneity (group=imut):
Log-rank test for equality of survivor functions
      | Events        Events
imut  | observed      expected
------+-------------------------
    0 |   3            6.91
    1 |  11            7.09
------+-------------------------
Total |  14           14.00
chi2(1) = 4.38
Pr>chi2 = 0.0363


Logrank test of homogeneity (group=Mutation):
Log-rank test for equality of survivor functions
         | Events      Events
Mutation | observed    expected
---------+-------------------------
      WT |  3            6.91
    BRAF |  4            1.88
     KIT |  1            1.74
MULTIHIT |  1            0.04
     NF1 |  2            1.08
    NRAS |  3            2.33
---------+-------------------------
   Total | 14           14.00
chi2(5) = 28.94
Pr>chi2 = 0.0000
``` 
Cox proportional hazards analysis testing the association of any mutation on overall survival. (Supplementary Table 26)

```
generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
replace cstage = max(cstage,2)

stset t, failure(death) origin(time dxdate) enter(time condate)
failure event: death != 0 & death < .
obs. time interval: (origin, t]
enter on or after: time condate
exit on or before: failure
t for analysis: (time-origin)
origin: time dxdate

xi:stcox imut sex age cstage Q2AMR
failure _d: death
analysis time _t: (t-origin)
origin: time dxdate
enter on or after: time condate
Iteration 0: log likelihood = -41.052839
Iteration 1: log likelihood = -35.488137
Iteration 2: log likelihood = -35.332644
Iteration 3: log likelihood = -35.331498
Iteration 4: log likelihood = -35.331498
Refining estimates:
Iteration 0: log likelihood = -35.331498
Cox regression -- no ties
No. of subjects = 73                        Number of obs = 73
No. of failures = 12
Time at risk = 84082
                                            LR chi2(5) = 11.44
Log likelihood = -35.331498                 Prob > chi2 = 0.0433
------------------------------------------------------------------------------
          _t | Haz. Ratio    Std. Err.     z      P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
        imut | 5.067486     4.144197     1.98     0.047     1.020185  25.17134
         sex | 1.855816     1.160979     0.99     0.323      .5445432 6.324665
         age | 1.025344      .033649     0.76     0.446      .9614689 1.093462
      cstage | 3.519382     2.571598     1.72     0.085      .8404214 14.73791
       Q2AMR | .0651891      .1089418   -1.63     0.102      .0024642 1.724543
------------------------------------------------------------------------------

```


Chi squared test of independence, cluster assignment and survival. (Supplemetary Table 27)

``` 

generate ncluster = .
replace ncluster = 1 if cluster == 1
replace ncluster = 2 if cluster == 2
replace ncluster = 3 if cluster == 3

tab ncluster death, row all

+----------------+
| Key            |
|----------------|
| frequency      |
| row percentage |
+----------------+
           |        death
ncluster   |    0            1    | Total
-----------+----------------------+----------
         2 |    10           6    | 16
           |    62.50       37.50 | 100.00
-----------+----------------------+----------
         1 |    14           0    | 14
           |    100.00       0.00 | 100.00
-----------+----------------------+----------
         3 |    11           3    | 14
           |    78.57       21.43 | 100.00
-----------+----------------------+----------
     Total |    35           9    | 44
           |    79.55       20.45 | 100.00
Pearson chi2(2) = 6.4653 Pr = 0.039
likelihood-ratio chi2(2) = 8.8660 Pr = 0.012
Cramér’s V = 0.3833
gamma = -0.3514 ASE = 0.321
Kendall’s tau-b = -0.1732 ASE = 0.162

generate follow = t - dxdate
``` 