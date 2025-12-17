# Survival analysis 
## Author: Timothy Bishop
## Modified by Johana Itzel

Overall and recurrence free survival analysis were done using STATA 19.5. Supplementary Table 1 was used as input after filtering only patients with primary samples (85 patients). Further filtering was done for each analysis if necessary for patiens with missing data used in the analysis (e.g. ancestry, RNAseq cluster data)

## Recurrence

Multivariate logistic regression to test the association of any driver mutation on recurrnce-free survival (Supplementary Table 22) 

```
count
replace last_note = subinstr(last_note,"/",".", 5)
generate dln = date(last_note, "DMY")
format dln %td
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
 
******
*Import ancestry data
******
clear
keep tumor_sample_barcode
generate ID = substr(tumor_sample_barcode,1,7)
save tumor_sample_list.dta, replace
count
clear
import excel "Supplementary_Table_2.xlsx", sheet("Supplementary_Table_2") firstrow
keep if strpos(ID,"PD") > 0
replace ID = substr(ID,1,7)
merge 1:1 ID using tumor_sample_list.dta, generate(mm)
tab mm
keep if mm > 1.5
drop mm
count
merge 1:1 tumor_sample_barcode, generate(mm)
tab mm
drop mm
save amdata.dta, replace
*Note: tumor_sample_list.dta is a list of the saple used per patient.

***DATA READY FOR ANALYSIS
****************************************************************


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
stset t, failure(recurrence) origin(time dxdate) enter(time condate)
xi: logistic recurrence imut sex age cstage dxdate Q2AMR

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
generate follow = t - dxdate

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
       WT  |   39            4    | 43
           |   90.70         9.30 | 100.00
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
     Total |   70           15    | 85
           | 82.35          17.65 | 100.00

          Pearson chi2(5) =  10.9539   Pr = 0.052
 Likelihood-ratio chi2(5) =   9.6104   Pr = 0.087
               Cramér's V =   0.3590
                    gamma =   0.3504  ASE = 0.168
          Kendall's tau-b =   0.1782  ASE = 0.096
           Fisher's exact =                 0.042

tab imut death, all exact row
+----------------+
| Key            |
|----------------|
| frequency      |
| row percentage |
+----------------+

        |         death
      imut |         0          1 |     Total
-----------+----------------------+----------
         0 |        39          4 |        43 
           |     90.70       9.30 |    100.00 
-----------+----------------------+----------
         1 |        31         11 |        42 
           |     73.81      26.19 |    100.00 
-----------+----------------------+----------
     Total |        70         15 |        85 
           |     82.35      17.65 |    100.00 

          Pearson chi2(1) =   4.1698   Pr = 0.041
 Likelihood-ratio chi2(1) =   4.3015   Pr = 0.038
               Cramér's V =   0.2215
                    gamma =   0.5515  ASE = 0.220
          Kendall's tau-b =   0.2215  ASE = 0.100
           Fisher's exact =                 0.050
   1-sided Fisher's exact =                 0.038

Logrank test of homogeneity (group=imut):

Equality of survivor functions
Log-rank test

      |  Observed       Expected
imut  |    events         events
------+-------------------------
    0 |         4           7.43
    1 |        11           7.57
------+-------------------------
Total |        15          15.00

                chi2(1) =   3.13
                Pr>chi2 = 0.0766


Logrank test of homogeneity (group=Mutation):

Equality of survivor functions
Log-rank test

         |  Observed       Expected
Mutation |    events         events
---------+-------------------------
      WT |         4           7.43
    BRAF |         4           2.02
     KIT |         1           1.85
MULTIHIT |         1           0.05
     NF1 |         2           1.17
    NRAS |         3           2.48
---------+-------------------------
   Total |        15          15.00

                   chi2(5) =  21.88
                   Pr>chi2 = 0.0006
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
Failure _d: death
   Analysis time _t: (t-origin)
             Origin: time dxdate
  Enter on or after: time condate

Iteration 0:  Log likelihood = -43.760889
Iteration 1:  Log likelihood = -38.280201
Iteration 2:  Log likelihood = -38.182057
Iteration 3:  Log likelihood = -38.181793
Refining estimates:
Iteration 0:  Log likelihood = -38.181793

Cox regression with no ties

No. of subjects =     73                                Number of obs =     73
No. of failures =     13
Time at risk    = 83,920
                                                        LR chi2(5)    =  11.16
Log likelihood = -38.181793                             Prob > chi2   = 0.0483

------------------------------------------------------------------------------
          _t | Haz. ratio   Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
        imut |   3.185216   2.253241     1.64   0.101      .796148    12.74336
         sex |   1.573245   .9712176     0.73   0.463     .4691595    5.275602
         age |   1.027321    .032233     0.86   0.390     .9660492     1.09248
  _Icstage_3 |   4.051354   2.897533     1.96   0.050     .9972859    16.45813
       Q2AMR |   .0445928   .0713461    -1.94   0.052     .0019382    1.025987
------------------------------------------------------------------------------
i.cstage          _Icstage_2-3        (naturally coded; _Icstage_2 omitted)

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
|   frequency    |
| row percentage |
+----------------+

           |         death
  ncluster |         0          1 |     Total
-----------+----------------------+----------
         1 |         9          7 |        16 
           |     56.25      43.75 |    100.00 
-----------+----------------------+----------
         2 |        14          0 |        14 
           |    100.00       0.00 |    100.00 
-----------+----------------------+----------
         3 |        11          3 |        14 
           |     78.57      21.43 |    100.00 
-----------+----------------------+----------
     Total |        34         10 |        44 
           |     77.27      22.73 |    100.00 

          Pearson chi2(2) =   8.1576   Pr = 0.017
 Likelihood-ratio chi2(2) =  10.6862   Pr = 0.005
               Cramér's V =   0.4306
                    gamma =  -0.4344  ASE = 0.290
          Kendall's tau-b =  -0.2265  ASE = 0.158
```

```
Likelihood-ratio test statistic of homogeneity (group=ncluster):
chi2( 2 ) = 13.783539,   P = .00101611
 
Logrank test of homogeneity (group=ncluster):

Equality of survivor functions
Log-rank test

         |  Observed       Expected
ncluster |    events         events
---------+-------------------------
       1 |         7           3.10
       2 |         0           3.99
       3 |         3           2.91
---------+-------------------------
   Total |        10          10.00

                   chi2(2) =   8.99
                   Pr>chi2 = 0.0111
```
Sex and mutational status association

```
. logistic imut dxdate fsex age Q2AMR cstage
Logistic regression Number of obs = 73
LR chi2(5) = 13.16
Prob > chi2 = 0.0219
Log likelihood = -43.849075 Pseudo R2 = 0.1305
------------------------------------------------------------------------------
imut | Odds ratio Std. err. z P>|z| [95% conf. interval]
-------------+----------------------------------------------------------------
dxdate | 1.000281 .0003449 0.82 0.415 .9996054 1.000957
fsex | 3.817422 2.066026 2.48 0.013 1.321574 11.02678
age | .9646099 .0228448 -1.52 0.128 .9208582 1.01044
Q2AMR | .0625818 .0849547 -2.04 0.041 .0043746 .8952685
cstage | 1.597506 .8650562 0.87 0.387 .5527302 4.617126
_cons | .0277799 .194587 -0.51 0.609 3.03e-08 25471.06
------------------------------------------------------------------------------
```
