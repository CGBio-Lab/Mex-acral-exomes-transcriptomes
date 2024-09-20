# Acral melanoma exomes' analysis
# Survival analysis
## Author: Timothy Bishop

Overall and recurrence free survival analysis were done using STATA

```
import delimited SumplementaryTable.csv", delimiter("", collapse) varnames(1) encoding(ISO-8859-1)clear
generate dln = date(date_of_last_note, "MDY")
format dln %td
keep tumor_sample_barcode patient_file date_of_recruitment dln
save consent.dta, replace

import delimited SumplementaryTable.csv", delimiter("", collapse) varnames(1) encoding(ISO-8859-1)clear
merge 1:1 tumor_sample_barcode patient_file using consent.dta, generate(mm)
tab mm
drop mm
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
replace imut = 0 if mutation_status == "other"
generate delay = condate-dxdate
summ delay
tab1 mutation ibraf ihras ikit ikras inf1 inras
save amdata.dta, replace

**************************************************************
*** PROSPECTIVE ANALYSIS OF RECURRENCE
**************************************************************

log using analysis_prospective_relapse.log, replace
use amdata.dta
replace lndate = dthdate if dthdate != .
generate rectime = recdate-condate
drop if rectime <= 0 & rectime != .
replace rectime = lndate-condate if recurrence == 0
export excel using "dataset", sheet("Recurrence") sheetreplace firstrow(variables)

ltable rectime recurrence, test notable graph xlab(0(400)2000) ylab(0(0.2)1)  xtitle("Since recruitment (days)") ytitle("Prop. Relapse-free") ci
graph export ltable_prec.pdf, as(pdf) replace

ltable rectime recurrence, test notable graph by(imut) overlay xlab(0(400)2000) ylab(0(0.2)1)  xtitle("Since recruitment (days)") ytitle("Prop. Relapse-free") ci
graph export ltable_imut_prec.pdf, as(pdf) replace
ltable rectime recurrence, test notable graph by(mutation) overlay xlab(0(400)2000) ylab(0(0.2)1)  xtitle("Since recruitment (days)") ytitle("Prop. Relapse-free") ci
graph export ltable_mut_prec.pdf, as(pdf) replace
destring rna_cluster, replace ignore("NA")
ltable rectime recurrence, test notable graph by(rna_cluster) overlay xlab(0(400)2000) ylab(0(0.2)1)  xtitle("Since recruitment (days)") ytitle("Prop. Relapse-free") ci
graph export ltable_rna_prec.pdf, as(pdf) replace

generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
replace cstage = max(cstage,2)
bysort recurrence: summ rectime
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" "mut" {
tab1 i`gn'
ltable rectime recurrence, test notable graph overlay xlab(0(400)2000) ylab(0(0.2)1) by(i`gn') xtitle("Since recruitment (days)") ytitle("Prop. Relapse-free") ci
graph export ltable_`gn'_prec.pdf, as(pdf)replace
}
* Now move onto cox regression
stset rectime, failure(recurrence)
stdescribe
stcox imut
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" "mut" {
display "`gn'"
tab1 i`gn'
xi: stcox i`gn' sex age i.cstage
}
destring rna_cluster, replace ignore("NA")
xi: stcox i.rna_cluster
xi: stcox i.rna_cluster sex age i.cstage
log close
graphlog using analysis_prospective_relapse.log, replace keeptex

**************************************************************
*** RETROSPECTIVE ANALYSIS OF RECURRENCE
**************************************************************
use amdata.dta
log using analysis_retrospective_relapse.log

generate diff1 = recdate-condate
sort diff1
list patient_file dxdate recdate condate diff1 dthdate lndate
list lndate dthdate if dthdate != .
replace lndate = dthdate if dthdate != .
display "Consent date preceeds recurrence date"
count if diff1 >= 0 & diff1 != .
display "Consent date AFTER recurrence date "
count if diff1 <0 & diff1 != .
keep if diff1 <= 0
sort dxyear

generate ttr = recdate-dxdate
bysort mutation: summ ttr
bysort imut: summ ttr
ttest ttr, by(imut)
generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
replace cstage = max(cstage,2)
export excel using "dataset", sheet("Retrospective") sheetreplace firstrow(variables)

tab cstage stage
xi:reg ttr imut sex age i.cstage
log close
graphlog using analysis_retrospective_relapse.log, replace keeptex


**************************************************************
*** PROSPECTIVE ANALYSIS OF DEATH
**************************************************************
use amdata.dta
log using analysis_prospective_death.log, replace

drop imut
count
generate imut = 1
replace imut = 0 if mutation_status == "other"
list lndate dthdate if dthdate != .
replace lndate = dthdate if dthdate != .
generate studytime = dthdate-condate
replace studytime = lndate-condate if death == 0
tab mutation death
export excel using "dataset", sheet("Death") sheetreplace firstrow(variables)

ltable studytime death, test notable graph overlay xlab(0(400)2000) ylab(0(0.2)1) by(imut) xtitle("Since recruitment (days)") ytitle("Proportion Alive") ci
graph export ltable_death_imut.pdf, as(pdf) replace

foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" "mut" {
tab1 i`gn'
ltable studytime death, test notable graph overlay xlab(0(400)2000) ylab(0(0.2)1) by(i`gn') xtitle("Since recruitment (days)") ytitle("Proportion Alive") ci
graph export ltable_`gn'_prec.pdf, as(pdf)replace
}

generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
replace cstage = max(cstage,2)
destring rna_cluster, replace ignore("NA")
tab rna_cluster death, all
generate irna_cluster = 3-rna_cluster
xi: logistic death i.irna_cluster sex age i.cstage
ltable studytime death, test notable graph overlay xlab(0(400)2000) ylab(0(0.2)1) by(rna_cluster) xtitle("Since recruitment (days)") ytitle("Proportion Alive") ci
graph export ltable_death_rna.pdf, as(pdf) replace

stset studytime, failure(death)
stcox imut
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" "mut" {
display "`gn'"
xi: stcox i`gn' sex age i.cstage
}
destring rna_cluster, replace ignore("NA")
tab rna_cluster death, all
xi: stcox i.irna_cluster
xi: stcox i.irna_cluster sex age i.cstage

log close
graphlog using analysis_prospective_death.log, replace keeptex

logistic death age sex
use amdata.dta
tab year imut
generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
generate rectime = recdate-dxdate
replace rectime = lndate-dxdate if recurrence == 0
keep if tumor_type == "primary"
/*
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" {
display "`gn'"
logistic i`gn' year
ranksum rectime, by(i`gn')
}
bysort patient_file: egen ct = count(1)
bysort patient_file: egen ord = rank(dxdate)
keep if ord == 1
*drop if dxyear < 2016
bysort recurrence: summ rectime 
ltable rectime recurrence, by(inf1) graph overlay test notable

stset rectime, failure(recurrence)
stcox inf1, strata(dxyear)
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" "mut" {
display "`gn'"
xi: stcox i`gn' sex age i.cstage, strata(dxyear)
}

destring rna_cluster, replace ignore("NA")
xi: stcox i.rna_cluster
xi: stcox i.rna_cluster sex age i.cstage, strata(dxyear)

 program define mylogit
        args lnf Xb
		quietly replace `lnf' = -ln(1+exp(-`Xb')) if $ML_y1==1
		quietly replace `lnf' = -`Xb' - ln(1+exp(-`Xb')) if $ML_y1==0
 end

clear
use amdata.dta  
  ml model lf mylogit (death=age sex)
  ml maximize
  
****Rewrite

***Look at year of diagnosis and recurrence
clear
use amdata.dta
keep if tumor_type == "primary" 
tab year imut
generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
generate rectime = recdate-dxdate
replace rectime = lndate-dxdate if recurrence == 0
keep if tumor_type == "primary"
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" {
display "`gn'"
logistic i`gn' year
ranksum rectime, by(i`gn')
}
bysort patient_file: egen ct = count(1)
bysort patient_file: egen ord = rank(dxdate)
keep if ord == 1
*drop if dxyear < 2016
bysort recurrence: summ rectime 
ltable rectime recurrence, by(inf1) graph overlay test notable
generate condate = date("1/1/2017","DMY")
replace condate = max(condate, dxdate)
generate studydate = lndate
replace studydate = recdate if recurrence == 1
stset studydate, failure(recurrence) origin(dxdate) enter(condate) exit(lndate)
stcox inf1, strata(dxyear)
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" "mut" {
display "`gn'"
xi: stcox i`gn' sex age i.cstage
}

destring rna_cluster, replace ignore("NA")
xi: stcox i.rna_cluster
xi: stcox i.rna_cluster sex age i.cstage

***Look at year of diagnosis and death
clear
use amdata.dta
keep if tumor_type == "primary" 
duplicates drop patient_file, force
tab dxyear mutation, row
generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
generate dthtime = dthdate-dxdate
replace dthtime = lndate-dxdate if death == 0

foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" {
display "`gn'"
logistic i`gn' dxdate
ranksum dxdate, by(i`gn')
}
*drop if dxyear < 2016
bysort death: summ dthtime 
ltable dthtime death, by(inf1) graph overlay test notable
generate condate = date("1/1/2016","DMY")
replace condate = max(condate, dxdate)
format %td condate
generate studydate = lndate
replace studydate = dthdate if death == 1
stset studydate, failure(death) origin(dxdate) enter(condate) exit(lndate)
stcox inf1, strata(dxyear)
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" "mut" {
display "`gn'"
xi: stcox i`gn' sex age i.cstage
}

destring rna_cluster, replace ignore("NA")
xi: stcox i.rna_cluster
xi: stcox i.rna_cluster sex age i.cstage


***Look at time until recurrence
clear
use amdata.dta
keep if tumor_type == "primary" 
duplicates drop patient_file, force
tab dxyear mutation, row
generate rectime = recdate-dxdate
replace rectime = lndate-dxdate if recurrence == 0
keep if recurrence == 1
drop if dxyear < 2016
foreach gn in "braf" "hras" "kit" "kras" "nf1" "nras" "mut" {
display "`gn'"
ranksum rectime, by(i`gn')
}

***Look at survival from 2016

clear
use amdata.dta
keep if tumor_type == "primary" 
duplicates drop patient_file, force
drop if dxyear < 2016
count
generate studytime = (dthdate-dxdate)/365.25
replace studytime = (lndate-dxdate)/365.25 if death == 0
ltable studytime death, graph test notable
destring rna_cluster, replace ignore("NA")
ltable studytime death, by(rna_cluster) overlay graph test notable
generate cstage = stage
destring cstage, replace ignore("NA")
replace cstage = . if cstage == 0
replace cstage = min(cstage,3)
tab cstage rna_cluster, all



clear
use amdata.dta
keep if tumor_type == "primary" 
duplicates drop patient_file, force
drop if dxyear < 2017
count
generate studytime = (dthdate-dxdate)/365.25
replace studytime = (lndate-dxdate)/365.25 if death == 0
ltable studytime death, graph test notable
destring rna_cluster, replace ignore("NA")
ltable studytime death, by(rna_cluster) overlay graph test notables

```