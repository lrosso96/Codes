**********************************
*** RA APPLICATION CODING TASK ***
**********************************

* Date: 10-10-2020
* Author: Lucas Rosso

clear all
set more off, permanently
*cd "C:\...\Coding_Task
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task"
import delimited RA_21_22.csv

* Assumption: consider only households with positive total assets 
*keep if asset_total>=0 // 7 households with asset_total<0

* Generating Wealth
g wealth_total = asset_total - debt_total
g wealth_housing = asset_housing - debt_housing
* --------------------------------------- *

save DATA.dta, replace

******************
*** QUESTION 1 ***
******************

** DESCRIPTIVE STATISTICS ** 
tabstat wealth_total wealth_housing [aw=weight], stat(median) by(education)
tabstat wealth_total wealth_housing [aw=weight], stat(median) by(race)
* --------------------------------------- *

* TRENDS IN TOTAL WEALTH BY RACE
collapse (median) wealth_total [aw=weight], by(year race)
rename wealth_total median_wealth_total
replace median_wealth_total = median_wealth_total/1000 // in thousand US$

*cd "C:\...\Coding_Task\Figures
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task\Figures"

** Figure 1: Median Wealth by Race
twoway (connected median_wealth_total year if race=="black", msymbol(smcircle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) (connected median_wealth_total year if ///
race=="Hispanic", msymbol(smsquare) mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected median_wealth_total year if race=="other", msymbol(smdiamond) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) (connected median_wealth_total year ///
if race=="white", msymbol(smtriangle) mfcolor(%65) mlcolor(%75) lcolor(%65)), ///
ytitle("Total Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "Black" 2 "Hispanic" 4 "White" 3 "Other") rows(1) ///
region(fcolor(white) lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export wealth_byrace.pdf, replace
* --------------------------------------- *

** Figures 2-5: Median Wealth by Race (Fot Appendix)
* Black
twoway (line median_wealth_total year if race=="black", lcolor(%60) lwidth(medthick) ///
lpattern(dash)), ytitle(" Total Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle(" ") ///
graphregion(fcolor(white) lcolor(white))
gr export t_wealth_black.pdf, replace

* Hispanic
twoway (line median_wealth_total year if race=="Hispanic", lcolor(%60) lwidth(medthick) ///
lpattern(dash)), ytitle(" Total Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle(" ") ///
graphregion(fcolor(white) lcolor(white))
gr export t_wealth_hispanic.pdf, replace

* White
twoway (line median_wealth_total year if race=="white", lcolor(%60) lwidth(medthick) ///
lpattern(dash)), ytitle(" Total Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle(" ") ///
graphregion(fcolor(white) lcolor(white))
gr export t_wealth_white.pdf, replace

* Other
twoway (line median_wealth_total year if race=="other", lcolor(%60) lwidth(medthick) ///
lpattern(dash)), ytitle(" Total Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle(" ") ///
graphregion(fcolor(white) lcolor(white))
gr export t_wealth_other.pdf, replace
* --------------------------------------- *

* TRENDS IN TOTAL WEALTH BY EDUCATIONAL LEVEL

*cd "C:\...\Coding_Task
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task"
use DATA, clear

collapse (median) wealth_total [aw=weight], by(year education)
rename wealth_total median_wealth_total
replace median_wealth_total = median_wealth_total/1000 // in thousand US$

*cd "C:\...\Coding_Task\Figures
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task\Figures"

** Figure 6: Median Wealth by Educational Level
twoway (connected median_wealth_total year if education=="college degree", msymbol(smcircle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) (connected median_wealth_total year if ///
education=="no college", msymbol(smsquare) mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected median_wealth_total year if education=="some college", msymbol(smtriangle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)), ytitle("Total Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "College Degree" 2 "No College" 3 "Some College") rows(1) ///
region(fcolor(white) lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export wealth_byeducation.pdf, replace
* --------------------------------------- *

******************
*** QUESTION 2 ***
******************

* Same figures as before, but for housing wealth and only for black and white households.
*cd "C:\...\Coding_Task
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task"
use DATA, clear

* Many zeros for black households
tabstat wealth_housing [aw=weight] if race=="black" | race=="white", stat(median) by(race)
count if wealth_housing ==0 & race=="black"
count if wealth_housing <0 & race=="black"

* TRENDS IN TOTAL WEALTH BY RACE
collapse (median) wealth_housing [aw=weight], by(year race)
rename wealth_housing median_wealth_housing
replace median_wealth_housing = median_wealth_housing/1000 // in thousand US$

*cd "C:\...\Coding_Task\Figures
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task\Figures"

** Figure 7: Median Wealth by Race
twoway (connected median_wealth_housing year if race=="black", msymbol(smcircle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) (connected median_wealth_housing year if ///
race=="white", msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)), ///
ytitle("Housing Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "Black" 2 "White") rows(1) region(fcolor(white) lcolor(white)) ///
bexpand) graphregion(fcolor(white) lcolor(white))
gr export housing_wealth_byrace.pdf, replace
* --------------------------------------- *

* TRENDS IN TOTAL WEALTH BY EDUCATIONAL LEVEL
*cd "C:\...\Coding_Task
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task"
use DATA, clear

collapse (median) wealth_housing [aw=weight], by(year education race)
rename wealth_housing median_wealth_housing
replace median_wealth_housing = median_wealth_housing/1000 // in thousand US$

*cd "C:\...\Coding_Task\Figures
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task\Figures"

** Figure 8: Median Wealth by Educational Level for black households
twoway (connected median_wealth_housing year if education=="college degree", msymbol(smcircle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) (connected median_wealth_housing year if ///
education=="no college", msymbol(smsquare) mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected median_wealth_housing year if education=="some college", msymbol(smtriangle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) if race=="black", ytitle("Housing Wealth (Median)") ///
ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "College Degree" 2 "No College" 3 "Some College") rows(1) ///
region(fcolor(white) lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export h_wealth_byeducation_black.pdf, replace
* --------------------------------------- *

** Figure 9: Median Wealth by Educational Level for white households
twoway (connected median_wealth_housing year if education=="college degree", msymbol(smcircle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) (connected median_wealth_housing year if ///
education=="no college", msymbol(smsquare) mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected median_wealth_housing year if education=="some college", msymbol(smtriangle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) if race=="white", ytitle("Housing Wealth (Median)") ///
ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "College Degree" 2 "No College" 3 "Some College") rows(1) ///
region(fcolor(white) lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export h_wealth_byeducation_white.pdf, replace
* --------------------------------------- *

** Figure 10: Median Wealth for College Educated Households
twoway (connected median_wealth_housing year if race=="black", msymbol(smcircle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected median_wealth_housing year if race=="white", msymbol(smcircle) ///
mfcolor(%65) mlcolor(%75) lcolor(%65)) if education=="college degree", ///
ytitle("Housing Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "Black" 2 "White") rows(1) ///
region(fcolor(white) lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export college_housing_wealth.pdf, replace
* --------------------------------------- *

******************
*** QUESTION 3 ***
******************

*cd "C:\...\Coding_Task
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task"
use DATA, clear

* Keeping homeowners age 25 or older
keep if age>24 & wealth_housing!=0

collapse (median) wealth_total wealth_housing [aw=weight], by(year race)
rename wealth_housing median_wealth_housing
rename wealth_total median_wealth_total
replace median_wealth_housing = median_wealth_housing/1000 // in thousand US$
replace median_wealth_total   = median_wealth_total/1000   // in thousand US$

* Generating Log median wealth
g l_median_wealth_housing = log(median_wealth_housing)
g l_median_wealth_total   = log(median_wealth_total)

*cd "C:\...\Coding_Task\Figures
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task\Figures"

* IN LEVEL

** Figure 11: Total Wealth for homeowners
twoway (connected median_wealth_total year if race=="black", ///
msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected median_wealth_total year if race=="white", ///
msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)), ///
xline(2007, lwidth(thin) lpattern(dash) lcolor(black)) ///
ytitle("Total Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "Black" 2 "White") region(fcolor(white) ///
lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export homeowners_twealth.pdf, replace
* --------------------------------------- *

** Figure 12: Housing Wealth for homeowners
twoway (connected median_wealth_housing year if race=="black", ///
msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected median_wealth_housing year if race=="white", ///
msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)), ///
xline(2007, lwidth(thin) lpattern(dash) lcolor(black)) ///
ytitle("Housing Wealth in 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "Black" 2 "White") region(fcolor(white) ///
lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export homeowners_hwealth.pdf, replace
* --------------------------------------- *

* IN LOGS

** Figure 13: Change in Total Wealth for homeowners
twoway (connected l_median_wealth_total year if race=="black", ///
msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected l_median_wealth_total year if race=="white", ///
msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)), ///
xline(2007, lwidth(thin) lpattern(dash) lcolor(black)) ///
ytitle("Total Wealth in Log 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "Black" 2 "White") region(fcolor(white) ///
lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export homeowners_log_twealth.pdf, replace

** Figure 13: Change in Housing Wealth for homeowners
twoway (connected l_median_wealth_housing year if race=="black", ///
msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)) ///
(connected l_median_wealth_housing year if race=="white", ///
msymbol(smcircle) mfcolor(%65) mlcolor(%75) lcolor(%65)), ///
xline(2007, lwidth(thin) lpattern(dash) lcolor(black)) ///
ytitle("Housing Wealth in Log 2016 US$") ylabel(, nogrid angle(horizontal)) xtitle("") ///
legend(order(1 "Black" 2 "White") region(fcolor(white) ///
lcolor(white)) bexpand) graphregion(fcolor(white) lcolor(white))
gr export homeowners_log_hwealth.pdf, replace
* --------------------------------------- *

* Table 2: Change in Wealth after 2007
sort race year
keep if inlist(year, 2007,2010) 
keep if race=="black" | race=="white"
drop l_median_wealth_total l_median_wealth_housing  
by race: g delta_wealth_total = median_wealth_total[_n]/median_wealth_total[_n -1] -1
by race: g delta_wealth_housing = median_wealth_housing[_n]/median_wealth_housing[_n -1] -1

*cd "C:\...\Coding_Task\Tables
cd "C:\Users\Lucas Rosso\Desktop\Lucas Rosso\Postulaciones\Coding_Task\Tables"
export excel using table2, replace
* --------------------------------------- *







