////////////////////
/// SCF Analysis ///
////////////////////

* Date: 16-08-2020
* Author: Lucas Rosso

clear all
set more off, permanently
gl main "C:\Users\Lucas Rosso\Desktop\ME\Tesis\TESIS\Codigos\Stata"

cd "${main}\Data"
use all_tmp.dta

*********************************************
*** EDITING DATA AND GENERATING VARIABLES ***
*********************************************

keep if inrange(age, 21,70) & financial_wealth>0 & income>-1000000 
*keep if inrange(age, 21,70) & income>-1000000 
replace part = part*100
replace wage_risk = wage_risk*100
bys age: egen mean_part = mean(part) // Participation Rate
g risky_share = (Risky_assets/( Risky_assets + Safe_assets ))*100

bys age: egen mean_risky_share = mean(risky_share) // Unconditional Share
by age: egen mean_risky_share_cond = mean(risky_share) if risky_share>0 // Conditional Share
g debt = credit_cards_debt +  other_consumer_loans + education_loan // Debt

g Total_wealth = financial_wealth_houseNH - debt
g W_to_inc = Total_wealth/income
g mean_W_to_inc=.
g mean_wage_risk=.

* Percentiles *
xtile quant_inc = income if income>-1000000, nquantiles(5)
xtile quant_w = financial_wealth, nquantiles(5)
xtile FW_percentile = financial_wealth, nquantiles(100)
* ------------

 quietly forvalues i = 21/70 { 
 sum risky_share [aw=wgt] if age == `i', detail 
 replace mean_risky_share = r(mean) if age == `i' 
 sum risky_share [aw=wgt] if age == `i' & risky_share>0, detail 
 replace mean_risky_share_cond = r(mean) if age == `i' 
 sum part [aw=wgt] if age == `i', detail
 replace mean_part = r(mean) if age == `i'
 sum W_to_inc [aw=wgt] if age==`i' & income>0, detail
 replace mean_W_to_inc = r(p50) if age==`i'
 sum wage_risk [aw=wgt] if age==`i', detail
 replace mean_wage_risk = r(mean) if age==`i'
 }
 
bys age: g ID_age = _n

g age_group=.

qui forval i=1/10 {
replace age_group = `i' if age >= 16+ 5*`i' & age < 21 +5*`i'
}

g risky_share_5 =.
g risky_share_cond_5=.
g part_5=.
g wage_risk_5=.

qui forval i=1/10 {
	sum risky_share [aw=wgt] if age_group == `i', detail
	replace risky_share_5 = r(mean) if age_group == `i'
	sum risky_share [aw=wgt] if age_group==`i', detail
	replace risky_share_cond_5 = r(mean) if age_group == `i'
	sum risky_share [aw=wgt] if age_group == `i' & risky_share>0, detail 
	replace risky_share_cond_5 = r(mean) if age_group == `i'
	sum part [aw=wgt] if age_group == `i', detail 
	replace part_5 = r(mean) if age_group == `i'
	sum wage_risk [aw=wgt] if age_group == `i', detail
	replace wage_risk_5 = r(mean) if age_group == `i'
}
	
qui forval j=1/5{
	g part_q`j' = .
	g part_q`j'_5=.
}

qui forval i=21/70 {
	qui forval j=1/5 {
 		sum part [aw=wgt] if age==`i' & quant_w==`j'
		replace part_q`j' = r(mean) if age == `i'  
	}
}	

qui forval i=1/10 {
	qui forval j=1/5 {
		sum part_q`j' [aw=wgt] if age_group == `i', detail
		replace part_q`j'_5 = r(mean) if age_group == `i'
	}
}

g p75_=.
g p25_=.
g p90_=.
g p10_=.

qui forval i=21/70 {
	sum financial_wealth [aw=wgt] if age==`i', detail 
	replace p75_=r(p75) if age==`i'
	replace p25_=r(p25) if age==`i'
	replace p90_=r(p90) if age==`i'
	replace p10_=r(p10) if age==`i'
}

g p7525=p75_/p25_
g p9010=p90_/p10_ 

* moving average
g p7525_MA =.
g p9010_MA=.
qui forval i=23/68 {
	sum p7525 if inlist(age, `i'-2,`i'-1,`i',`i'+1,`i'+2)
	replace p7525_MA = r(mean) if inlist(age, `i'-1,`i',`i'+1)
	sum p9010 if inlist(age, `i'-2,`i'-1,`i',`i'+1,`i'+2)
	replace p9010_MA = r(mean) if inlist(age, `i'-1,`i',`i'+1)
}

* ----------------------------------------------------------- *

****************************************
***             FIGURES              ***
****************************************

cd "${main}\Figures_LR"

* STOCK SHARE
twoway (line mean_risky_share age, lcolor(green%75) lwidth(medthin)) (line mean_risky_share_cond age, lcolor(gs8%75) lwidth(medthin)) ///
(connected risky_share_cond_5 age if inlist(age, 23,28,33,38,43,48,53,58,63,68), mcolor(black) msize(small) msymbol(circle) lcolor(black)) ///
(connected risky_share_5 age if inlist(age, 23,28,33,38,43,48,53,58,63,68), mcolor(green) msize(small) msymbol(circle) lpattern(dash) lcolor(green)) ///
if ID_age==1, ytitle("Percent") ylabel(0(10)80, angle(horizontal) nogrid) xtitle("Age") ///
legend(order(3 "Conditional" 4 "Unconditional") rows(2) position(2) ring(0) region(fcolor(white))) ///
graphregion(fcolor(white) lcolor(white))
gr export stock_share.pdf, replace
* ----------------------------------- *

** PARTICIPATION RATE
twoway (connected part_5 age if inlist(age, 23,28,33,38,43,48,53,58,63,68), ///
mcolor(black) msize(small) msymbol(circle) lcolor(black)) ///
(line mean_part age, lcolor(gs8%75) lwidth(medthin)) if ID_age==1, ///
ytitle("Percent") ylabel(0(10)80, angle(horizontal) nogrid) xtitle("Age") ///
graphregion(fcolor(white) lcolor(white)) legend(off)
gr export part_rate.pdf, replace
* ----------------------------------- *

** PARTICIPATION RATE (5 YEAR AVG.)
twoway (connected part_q1_5 age, msize(small) msymbol(circle) lpattern(dash)) ///
(connected part_q3_5 age, msize(small) msymbol(triangle) lpattern(dash)) ///
(connected part_q5_5 age, mcolor(green) msymbol(X) lpattern(dash)) ///
if ID_age==1 & inlist(age, 28,33,38,43,48,53,58), ///
ytitle("Participation Share") ylabel(, angle(horizontal) ///
nogrid) xtitle("Age") legend(order(1 "Q1" 2 "Q3" 3 "Q5") rows(1) ///
region(fcolor(white) lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export part_rate_q_5Y.pdf, replace
* ----------------------------------- *

** P75/P25 (5-YEAR MOVING AVG.)
twoway (line p7525_MA age if age>20 & age<61, lcolor(%75) lwidth(medthin) lpattern(dash)) if ID_age==1, ///
ytitle("P75/P25") ylabel(, angle(horizontal) nogrid) xtitle("Age") ///
graphregion(fcolor(white) lcolor(white))
gr export p7525_MA.pdf, replace
* ----------------------------------- *

** P90/P10 (5-YEAR MOVING AVG.)
twoway (line p9010_MA age if age>20 & age<61, lcolor(%75) lwidth(medthin) lpattern(dash)) if ID_age==1, ///
ytitle("P90/P10") ylabel(, angle(horizontal) nogrid) xtitle("Age") ///
graphregion(fcolor(white) lcolor(white))
gr export p9010_MA.pdf, replace
* ----------------------------------- *
