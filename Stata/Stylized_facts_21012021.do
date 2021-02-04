/****************************************************************************************************
* Title: An Heterogeneous Agent Model of Portfolio Choice and its Implications for Wealth Inequality
* Created by: Lucas Rosso
* Created on: 23/11/2020
* Purpose: Use data from the SCF (1998-2019) to document stylized facts from US households
* Last Modified on: 03/02/2021
* Last Modified by: LR
* Edits:
	[23/11/2020]: Created dofile
	[27/11/2020]: Re-defined assets and produced some stylized facts
	[15/12/2020]: Computed Top 1% Wealth Share and Gini across surveys
	[21/01/2021]: New Baseline Definition, including pension wealth (financial_wealth_RA)
	[03/02/2021]: Exported risky share average to a csv file
*****************************************************************************************************/

clear all
set more off, permanently

* Useful Notation: 
* NH: financial wealth including housing net worth
* RA: financial wealth including retirement accounts

* Required Packages
*ssc install xtline
*ssc install inequal7

* Global macros for file paths (main must be changed)
gl main "C:\Users\LR\Desktop\ME\Tesis\TESIS\Codigos\Stata"
gl data "${main}\Data"
gl figures "${main}\Figures"
gl tables "${main}\Tables"

cd "$data"
use SCF_Data
* ----------------------------------------------------------- *

** Top 1% share and gini coefficient across surveys
gen wealth_aux = financial_wealth_RA*wgt // to capture weigths 

bys year: egen top_w        = sum(wealth_aux)
by year: egen top1          = sum(wealth_aux) if FW_percentile_RA == 100
by year: egen top5          = sum(wealth_aux) if FW_percentile_RA > 95 
by year: egen top10         = sum(wealth_aux) if FW_percentile_RA > 90
by year: egen bottom50      = sum(wealth_aux) if FW_percentile_RA <= 50
by year: egen middle40      = sum(wealth_aux) if FW_percentile_RA > 50 & FW_percentile_RA < 90
by year: gen top1_share     = top1/top_w 
by year: gen top5_share     = top5/top_w 
by year: gen top10_share    = top10/top_w 
by year: gen bottom50_share = bottom50/top_w 
by year: gen middle40_share = middle40/top_w 
bys year FW_percentile_RA: gen id_year_p = _n

gen gini_w = .
forval i = 1998(3)2019 {
	inequal7 financial_wealth_RA [aw=wgt] if year == `i', f(%9.2f) returnscalars
	replace gini_w = r(gini) if year == `i'
}
* -----------------------------------------------------------------

****************************************
***             FIGURES              ***
****************************************

cd "$figures" 

/* I use modern default and computer modern font. To use modern default type:
net install scheme-modern, from("https://raw.githubusercontent.com/mdroste/stata-scheme-modern/master/")
set scheme modern

to get back to STATA's default type: 
set scheme s2color

computer modern font can be downloaded from https://www.fontsquirrel.com/fonts/computer-modern. The following link:
https://medium.com/the-stata-guide/stata-graphs-get-those-fonts-right-c38d35625142 explaint how to install the font.
Finally edit graph preferences for the specific font (I use "CMU serif") */
 
* Global macros for figure style
gl general_style ylabel(, glwidth(vthin) ) xlabel(, glwidth(vthin)) xtitle(, margin(small))
gl connected_style lwidth(thin) mcolor(%50)

* ID for figures 
bys FW_percentile_H: g ID_pctile_H = _n
bys FW_percentile_NH: g ID_pctile_NH = _n
bys FW_percentile: g ID_pctile = _n
bys FW_percentile_NH_RA: g ID_NH_RA = _n
bys FW_percentile_RA: g ID_RA = _n
bys FW_old_acc: g ID_old = _n
* ------------------------------------

***************************
*** WEALTH DISTRIBUTION ***
*************************** 

g fwgt = round(wgt) // rounding weights to the nearest integer to use them in histogram

* financial wealth
g financial_wealth_aux = financial_wealth_RA/1000
replace financial_wealth_aux = 500 if financial_wealth_aux>=500

histogram financial_wealth_aux [fw=fwgt], bin(50) $general_style ytitle("Density") ylabel(, format(%9.2fc)) ///
fcolor(blue%70) lcolor(white) xtitle("Financial Wealth (Thousands)") xlabel(, format(%9.0fc))
gr export fin_wealth_dist.pdf, replace

* risky wealth
g Risky_assets_aux = Risky_assets_RA/1000
replace Risky_assets_aux = 500 if Risky_assets_aux>=500

histogram Risky_assets_aux [fw=fwgt], bin(50) $general_style ytitle("Density") ylabel(, format(%9.2fc)) ///
fcolor(red%70) lcolor(white) xtitle("Risky Wealth (Thousands)") xlabel(, format(%9.0fc))
gr export risky_wealth_dist.pdf, replace

/*
* risky wealth below/above median income
xtile income_pctile = income [aw=wgt], nquantiles(100)  

* below
histogram Risky_assets_aux [fw=fwgt] if income_pctile<=50, bin(50) $general_style ytitle("Density") ylabel(, format(%9.2fc)) ///
fcolor(red%70) lcolor(white) xtitle("Risky Wealth (Thousands)") xlabel(, format(%9.0fc))

* above
histogram Risky_assets_aux [fw=fwgt] if income_pctile>50, bin(50) $general_style ytitle("Density") ylabel(, format(%9.2fc)) ///
fcolor(red%70) lcolor(white) xtitle("Risky Wealth (Thousands)") xlabel(, format(%9.0fc))
*/

* safe wealth
g Safe_assets_aux = Safe_assets_RA/1000
replace Safe_assets_aux = 500 if Safe_assets_aux>=500

histogram Safe_assets_aux [fw=fwgt], bin(50) $general_style ytitle("Density") ylabel(, format(%9.2fc)) ///
fcolor(blue%70) lcolor(white) xtitle("Safe Wealth (Thousands)") xlabel(, format(%9.0fc))
gr export safe_wealth_dist.pdf, replace

* Lorenz curve
/*
ssc install lorenz
lorenz e financial_wealth Risky_assets Safe_assets income
lorenz graph, over
*/
* ------------------------------------


**************************************************
*** WEALTH INEQUALITY EVOLUTION ACROSS SURVEYS ***
**************************************************
twoway (connected top1_share year, $connected_style sort(year)) (connected gini_w year, yaxis(2) $connected_style sort(year)) ///
if FW_percentile_RA==100 & id_year_p==1, $general_style ytitle(Top 1% Wealth Share) ytitle(Wealth Gini Index, axis(2)) ///
ylabel(, angle(horizontal) format(%9.2fc) axis(2)) ylabel(, format(%9.2fc)) xtitle("") ///
legend(order(1 "Top 1% Wealth Share" 2 "Wealth Gini Index") rows(2) position(11) ring(0))
gr export w_ineq_evol.pdf, replace
* ------------------------------------

******************************************************
*** RISKY SHARE ACROSS WEALTH DIST (BASELINE DEF.) ***
******************************************************
twoway (connected risky_share_RA_p FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(0(0.1)0.6, format(%9.1f))
gr export baseline_riskyshare.pdf, replace
* ------------------------------------

*** RISKY SHARE ACROSS WEALTH DIST (BASELINE DEF.) ***
/*
twoway (connected risky_share_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)) ///
(connected cond_risky_share_p FW_percentile if ID_pctile == 1 & FW_percentile>20 , $connected_style sort(FW_percentile)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(0(0.1)0.6, format(%9.1f)) legend(order(1 ///
"Unconditional" 2 "Conditional") rows(2) region(fcolor(white)) position(11) ring(0))
gr export riskyshare_cond_uncond.pdf, replace
* ------------------------------------
*/

*************************************************
*** RISKY SHARE PART. RATE ACROSS WEALTH DIST ***
*************************************************
twoway (connected part_Risky_assets_RA_p FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Part. Rate) ylabel(, format(%9.1f))
gr export baseline_partrisky.pdf, replace
* ------------------------------------

*********************************************
*** RISKY SHARE FOR DIFFERENT DEFINITIONS ***
*********************************************
twoway (connected risky_share_RA_p FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)) ///
(connected risky_share_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)) ///
(connected risky_share_NH_RA_p FW_percentile_NH_RA if ID_NH_RA == 1, $connected_style sort(FW_percentile_NH_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(, format(%9.1f)) legend(order(1 ///
"Baseline" 2 "Exc. RA" 3 "Inc. Housing") rows(3) region(fcolor(white)) position(11) ring(0))
gr export robustness_share.pdf, replace 
* ------------------------------------

********************************************
*** PART. RATE FOR DIFFERENT DEFINITIONS ***
********************************************
twoway (connected part_Risky_assets_RA_p FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)) ///
(connected part_Risky_assets_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)) ///
(connected part_Risky_assets_houseNH_RA_p FW_percentile_NH_RA if ID_NH_RA == 1, $connected_style sort(FW_percentile_NH_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Part. Rate) ylabel(, format(%9.1f)) legend(order(1 ///
"Baseline" 2 "Exc. RA" 3 "Inc. Housing") rows(3) region(fcolor(white)) position(5) ring(0))
gr export robustness_part.pdf, replace 
* ------------------------------------

*********************************
*** ASSETS ACROSS WEALTH DIST ***
*********************************
twoway (connected house_share_p FW_percentile_NH_RA, $connected_style sort(FW_percentile_NH_RA)) ///
(connected stocks_plus_p FW_percentile_NH_RA, $connected_style sort(FW_percentile_NH_RA)) ///
(connected RA_share_p FW_percentile_NH_RA if RA_share_p<=1 , $connected_style sort(FW_percentile_NH_RA)) ///
if ID_NH_RA == 1, $general_style xtitle(Wealth Percentile) ytitle(Share) legend(order(1 ///
"Housing (net worth)" 2 "Risky" 3 "Ret. Accounts") rows(3) region(fcolor(white)) position(11) ring(0)) ylabel(, format(%9.1f))
gr export assets_share.pdf, replace
* ------------------------------------

*** OLD SHARES OUT OF RETIREMENT WEALTH ***
/*
twoway (connected share_old_p FW_old_acc if ID_old == 1, $connected_style sort(FW_old_acc)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(0(0.1)0.6, format(%9.1f))
gr export old_riskyshare.pdf, replace
*/
* ------------------------------------

*************************************************
*** BACK OF THE ENVELOPE RETURN HETEROGENEITY ***
*************************************************

g r_a = 0.067
g r_b = 0.021

g r_savings = r_a*risky_share_RA + (1-risky_share_RA)*r_b
g r_savings_p = .

qui forval i = 1/100 {
    sum r_savings [aw=wgt] if FW_percentile_RA == `i' 
	replace r_savings_p = r(mean) if FW_percentile_RA == `i'
}

twoway (connected r_savings_p FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Return on Savings) ylabel(, format(%9.2f))
gr export return_het.pdf, replace
* ------------------------------------



****************************************
***             TABLES               ***
****************************************

cd "$tables" 

rename part_pension_safe part_total_pension_safe
rename part_pension_risky part_total_pension_risky

* part. rate FW
g part_financial_wealth_RA = (financial_wealth_RA > 0)

local assets Safe_assets_RA checking_accounts savings_accounts money_market_accounts CDS bonds_safe mutual_safe ///
Risky_assets_RA brokerage stocks mutual_risky bonds_risky total_pension_risky IRA_risky total_pension_safe IRA_safe ///
net_house_wealth old_accounts financial_wealth_RA

*** HOUSEHOLDS ASSETS DESCRIPTIVE STATS ***
texdoc i "desc_stat.tex", replace force

qui foreach var of varlist `assets' {
sum `var' [aw=wgt], detail 
local `var'_mean =  r(mean)
local `var'_mean: di %9.0fc ``var'_mean'

local `var'_sd =  r(sd)
local `var'_sd: di %9.0fc ``var'_sd'

local `var'_p10 =  r(p10)
local `var'_p10: di %9.0fc ``var'_p10'

local `var'_p50 =  r(p50)
local `var'_p50: di %9.0fc ``var'_p50'

local `var'_p90 =  r(p90)
local `var'_p90: di %9.0fc ``var'_p90'

sum part_`var' [aw=wgt]
local `var'_part =  r(mean)
local `var'_part: di %9.2f ``var'_part'
}

***  CREATING LATEX TABLE ***
tex \begin{tabular}{lccccccc} \hline \hline
tex  & Mean & Sd & P10 & Median & P90 & Part. Rate \\  \midrule 
* total safe
tex \textit{Total Safe Assets (S)}&`Safe_assets_RA_mean'&`Safe_assets_RA_sd' & `Safe_assets_RA_p10' & `Safe_assets_RA_p50' & `Safe_assets_RA_p90' & `Safe_assets_RA_part' \\ 
* checking accounts
tex Checking Account & `checking_accounts_mean' & `checking_accounts_sd' & `checking_accounts_p10' & `checking_accounts_p50' ///
			& `checking_accounts_p90' & `checking_accounts_part' \\
* savings accounts		
tex Savings Accounts & `savings_accounts_mean' & `savings_accounts_sd' & `savings_accounts_p10' & `savings_accounts_p50' /// 
			& `savings_accounts_p90' & `savings_accounts_part' \\
* money market accounts
tex Money Market Accounts & `money_market_accounts_mean' & `money_market_accounts_sd' & `money_market_accounts_p10' ///
			& `money_market_accounts_p50' & `money_market_accounts_p90' & `money_market_accounts_part' \\
* cerificates of deposits
tex Cerificates of Deposits & `CDS_mean' & `CDS_sd' & `CDS_p10' & `CDS_p50' & `CDS_p90' & `CDS_part' \\
* savings_bond (safe)
tex Savings bond (safe) & `bonds_safe_mean' & `bonds_safe_sd' & `bonds_safe_p10' & `bonds_safe_p50' & `bonds_safe_p90' & `bonds_safe_part' \\ 
* mutual funds (safe)
tex Mutual Funds (safe) & `mutual_safe_mean' & `mutual_safe_sd' & `mutual_safe_p10' & `mutual_safe_p50' & `mutual_safe_p90' & `mutual_safe_part' \\ 
* IRA safe
tex IRA (safe) & `IRA_safe_mean' & `IRA_safe_sd' & `IRA_safe_p10' & `IRA_safe_p50' & `IRA_safe_p90' & `IRA_safe_part' \\
* pensions safe
tex Pensions (safe) & `total_pension_safe_mean' & `total_pension_safe_sd' & `total_pension_safe_p10' & `total_pension_safe_p50' ///
	& `total_pension_safe_p90' & `total_pension_safe_part' \\ 
tex \\ 
* total risky
tex \textit{Total Risky Assets (R)} & `Risky_assets_RA_mean' & `Risky_assets_RA_sd' & `Risky_assets_RA_p10' ///
		& `Risky_assets_RA_p50' & `Risky_assets_RA_p90' & `Risky_assets_RA_part' \\ 
* brokerage
tex Brokerage & `brokerage_mean' & `brokerage_sd' & `brokerage_p10' & `brokerage_p50' & `brokerage_p90' & `brokerage_part' \\ 
* stocks
tex Stocks & `stocks_mean' & `stocks_sd' & `stocks_p10' & `stocks_p50' & `stocks_p90' & `stocks_part' \\
* mutual funds (risky)
tex Mutual Funds (risky) & `mutual_risky_mean' & `mutual_risky_sd' & `mutual_risky_p10' & `mutual_risky_p50' & `mutual_risky_p90' & `mutual_risky_part' \\ 
* savings bond (risky)
tex Savings bond (risky) & `bonds_risky_mean' & `bonds_risky_sd' & `bonds_risky_p10' & `bonds_risky_p50' ///
		& `bonds_risky_p90' & `bonds_risky_part' \\ 
* IRA risky
tex IRA (risky) & `IRA_risky_mean' & `IRA_risky_sd' & `IRA_risky_p10' & `IRA_risky_p50' & `IRA_risky_p90' & `IRA_risky_part' \\
* pensions risky
tex Pensions (risky) & `total_pension_risky_mean' & `total_pension_risky_sd' & `total_pension_risky_p10' & `total_pension_risky_p50' ///
	& `total_pension_risky_p90' & `total_pension_risky_part' \\ 
tex \\
* baseline total FW
tex \textit{Baseline Financial Wealth (S + R)} & `financial_wealth_RA_mean' & `financial_wealth_RA_sd' & `financial_wealth_RA_p10' ///
		& `financial_wealth_RA_p50' & `financial_wealth_RA_p90' & `financial_wealth_RA_part' \\
* net house wealth
tex \\
tex \textit{Net House Wealth} & `net_house_wealth_mean'& `net_house_wealth_sd' & `net_house_wealth_p10' & `net_house_wealth_p50' & `net_house_wealth_p90' & `net_house_wealth_part' \\ \bottomrule
tex \end{tabular}
texdoc close 
* --------------------------------------------------------------------------------

*** DEMOGRAPHIC STATS ***
local demographic income educ age children marital

texdoc i "demo_stat.tex", replace force

qui foreach var of varlist `demographic' {
sum `var' [aw=wgt], detail 
local `var'_mean =  r(mean)
local `var'_mean: di %9.2fc ``var'_mean'

local `var'_sd =  r(sd)
local `var'_sd: di %9.2fc ``var'_sd'

local `var'_min =  r(min)
local `var'_min: di %9.0fc ``var'_min'

local `var'_max =  r(max)
local `var'_max: di %9.0fc ``var'_max'
}

***  CREATING LATEX TABLE ***
tex \begin{tabular}{lcccc} \hline \hline
tex  & Mean & Sd & Min & Max \\  \midrule 
tex Income & `income_mean' & `income_sd' & `income_min' & `income_max' \\
tex Age & `age_mean' & `age_sd' & `age_min' & `age_max' \\
tex Children & `children_mean' & `children_sd' & `children_min' & `children_max' \\
tex \$ \mathbb{1}\{\textnormal{College Degree}\} \$ & `educ_mean' & `educ_sd' & `educ_min' & `educ_max' \\
tex \$ \mathbb{1}\{\textnormal{Married}\} \$ & `marital_mean' & `marital_sd' & `marital_min' & `marital_max' \\ \bottomrule
tex \end{tabular}
texdoc close
* --------------------------------------------------------------------------------

*** INEQUALITY STATS ***
local wealth_stats top1_share top5_share top10_share bottom50_share gini_w

texdoc i "wealth_stats.tex", replace force

qui foreach var of varlist `wealth_stats' {
sum `var'
local `var'_mean =  r(mean)
local `var'_mean: di %9.2fc ``var'_mean'
}

***  CREATING LATEX TABLE ***
tex \begin{tabular}{l*{2}{>{\centering\arraybackslash}m{3cm}}} \hline \hline
tex  Measure     & Value  & Data \\ \midrule
tex Top 10\%  &   & `top10_share_mean' \\
tex Top 5\% & & `top5_share_mean' \\
tex Top 1\% & & `top1_share_mean' \\
tex Gini Index & & `gini_w_mean' \\ \bottomrule
tex \end{tabular}
texdoc close
* --------------------------------------------------------------------------------


**************************************************************
***               Risky Share Avg. by Decile               ***
**************************************************************

xtile FW_decile_RA     = financial_wealth_RA [aw=wgt], nquantiles(10)

matrix A = J(10,1,.)
qui forv i = 1/10 {
	sum risky_share_RA [aw=wgt] if FW_decile_RA == `i'
	mat A[`i',1] = r(mean)
}
matrix list A


********************************************
***         SENSITIVITY ANALYSIS         ***      
********************************************

cd "$figures"

****************************
*** COLLEGE - NO COLLEGE ***
****************************

g risky_share_c =. // college
g risky_share_nc=. // no college
qui forval i = 1/100 {
    sum risky_share_RA [aw=wgt] if FW_percentile_RA == `i' & educ == 1
	replace risky_share_c = r(mean) if FW_percentile_RA == `i' & educ == 1
	
	sum risky_share_RA [aw=wgt] if FW_percentile_RA == `i' & educ == 0
	replace risky_share_nc = r(mean) if FW_percentile_RA == `i' & educ == 0
}

egen tag_educ = tag(FW_percentile_RA educ)
replace tag_educ = 2 if tag_educ==1 & educ==1

twoway (connected risky_share_c FW_percentile_RA if tag_educ == 2, $connected_style sort(FW_percentile_RA)) ///
(connected risky_share_nc FW_percentile_RA if tag_educ == 1, $connected_style sort(FW_percentile_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(, format(%9.1f)) ///
legend(order(1 "College" 2 "No College") rows(2) position(11) ring(0))
gr export educ_robustness.pdf, replace

******************************
***  HOMEOWNERS - RENTERS  ***
******************************

g risky_share_owner_p   =.
g risky_share_noowner_p =.
qui forval i = 1/100 {
	sum risky_share_RA [aw=wgt] if FW_percentile_RA ==`i' & owner==1
	replace risky_share_owner_p = r(mean) if FW_percentile_RA ==`i' & owner==1
	
	sum risky_share_RA [aw=wgt] if FW_percentile_RA ==`i' & owner==0
	replace risky_share_noowner_p = r(mean) if FW_percentile_RA ==`i' & owner==0
}

egen id_homeowner = tag(FW_percentile_RA owner)
replace id_homeowner = id_homeowner + 1 if owner == 1 & id_homeowner==1

sort FW_percentile_RA owner
twoway (connected risky_share_owner_p FW_percentile_RA if id_homeowner == 2, $connected_style ) ///
(connected risky_share_noowner_p FW_percentile_RA if id_homeowner == 1, $connected_style ), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(0(0.1)0.6, format(%9.1f)) ///
legend(order(1 "Homeowners" 2 "Renters") rows(2) region(fcolor(white)) position(11) ring(0))
gr export homeowner_renter_rs.pdf, replace
* ------------------------------------

**************************
*** MEDIAN RISKY SHARE ***
**************************
/*
g risky_share_m=.
qui forval i = 1/100 {
    sum risky_share_RA [aw=wgt] if FW_percentile_RA == `i', detail
	replace risky_share_m = r(p50) if FW_percentile_RA == `i'
}

twoway (connected risky_share_m FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Median Risky Share) ylabel(0(0.1)0.6, format(%9.1f)) 
*/

***************************
*** VOLATILITY IN SHARE ***
***************************

g risky_share_sd=.
qui forval i = 1/100 {
    sum risky_share_RA [aw=wgt] if FW_percentile_RA == `i', detail
	replace risky_share_sd = r(sd) if FW_percentile_RA == `i'
}

twoway (connected risky_share_sd FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share Sd.) ylabel(0(0.1)0.4, format(%9.1f))


**********************************************
*** CONTROLLING FOR AGE, EDUCATION, INCOME ***
********************************************** 

g age2 = age*age
g age3 = age2*age
g age4 = age3*age

* risky share
reg risky_share_RA i.FW_percentile_RA i.educ i.owner income i.marital i.male children age4 age3 age2 age i.year [aw=wgt]
est sto reg1

predict risky_hat, xb

g risky_hat_p =.
qui forval i = 1/100 {
    sum risky_hat [aw=wgt] if FW_percentile_RA == `i'
	replace risky_hat_p = r(mean) if FW_percentile_RA == `i'
}

twoway (connected risky_hat_p FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Predicted Risky Share) ylabel(0(0.1)0.6, format(%9.1f))
gr export predicted_rs.pdf, replace 

coefplot reg1, drop(*.educ *.owner income *.marital *.male children age4 age3 age2 age *.year _cons) mfcolor(white) ciopts(recast(rcap)) level(95) vert /// 
xlabel(4 "10" 14 "20" 24 "30" 34 "40" 44 "50" 54 "60" 64 "70" 74 "80" 84 "90" 94 "100", grid) /// 
$general_style legend(off) xtitle(Wealth Percentile) ytitle(Coefficient) ylabel(0(0.1)0.6, format(%9.1f))
gr export coefplot_riskyshare.pdf, replace

* part rate
/*
probit part_Risky_assets i.FW_percentile i.educ income age3 age2 age i.year i.male children i.marital [fw=fwgt]

predict part_hat

g part_hat_p =.
qui forval i = 1/100 {
    sum part_hat [aw=wgt] if FW_percentile == `i'
	replace part_hat_p = r(mean) if FW_percentile == `i'
}

twoway (connected part_hat_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)), ///
$general_style xtitle(Wealth Percentile) ytitle(Predicted Part. Rate) ylabel(, format(%9.1f))

