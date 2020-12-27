/****************************************************************************************************
* Title: An Heterogeneous Agent Model of Portfolio Choice and its Implications for Wealth Inequality
* Created by: Lucas Rosso
* Created on: 23/11/2020
* Purpose: Use data from the SCF (1998-2019) to document stylized facts from US households
* Last Modified on: 23/11/2020
* Last Modified by: LR
* Edits:
	[23/11/2020]: Created dofile
	[27/11/2020]: Re-defined assets and produced some stylized facts
	[15/12/2020]: Computed Top 1% Wealth Share and Gini across surveys
*****************************************************************************************************/

clear all
set more off, permanently

* Useful Notation: 
* NH: includes housing net worth
* RA: includes retirement accounts
* FW: Financial wealth

* Required Packages
*ssc install xtline
*ssc install inequal7

* Global macros for file paths (main must be changed)
gl main "C:\Users\Lucas Rosso\Desktop\ME\Tesis\TESIS\Codigos\Stata"
gl data "${main}\Data"
gl figures "${main}\Figures"
gl tables "${main}\Tables"

cd "$data"
use SCF_Data
* ----------------------------------------------------------- *

gen wealth_aux = financial_wealth*wgt // to capture weigths 

bys year: egen top_w    = sum(wealth_aux)
by year: egen top1      = sum(wealth_aux) if FW_percentile == 100
by year: gen top1_share = top1/top_w 
*by year: egen risky_1   = sum()
bys year FW_percentile: gen id_year_p = _n

gen gini_w = .
forval i = 1998(3)2019 {
	inequal7 financial_wealth [aw=wgt] if year == `i', f(%9.2f) returnscalars
	replace gini_w = r(gini) if year == `i'
}

****************************************
***             FIGURES              ***
****************************************

cd "$figures" 

* Global macros for figure style
gl general_style ylabel(, angle(horizontal) glwidth(vthin) glcolor(gs14)) graphregion(fcolor(white) lcolor(white)) xlabel(, grid)
gl connected_style lwidth(thin) mcolor(%50)

* ID for figures 
bys FW_percentile_H: g ID_pctile_H = _n
bys FW_percentile_NH: g ID_pctile_NH = _n
bys FW_percentile: g ID_pctile = _n
bys FW_percentile_NH_RA: g ID_NH_RA = _n
bys FW_percentile_RA: g ID_RA = _n
bys FW_old_acc: g ID_old = _n
* ------------------------------------

*** WEALTH INEQUALITY EVOLUTION ACROSS SURVEYS ***
twoway (connected top1_share year, $connected_style sort(year)) (connected gini_w year, yaxis(2) $connected_style sort(year)) ///
if FW_percentile==100 & id_year_p==1, $general_style ytitle(Top 1% Wealth Share) ytitle(Wealth Gini Index, axis(2)) ///
ylabel(, angle(horizontal) format(%9.2fc) axis(2)) ylabel(, format(%9.2fc)) xtitle("") ///
legend(order(1 "Top 1% Wealth Share" 2 "Wealth Gini Index") rows(2) position(11) ring(0))
gr export w_ineq_evol.pdf, replace
* ------------------------------------

*** RISKY SHARE ACROSS WEALTH DIST (BASELINE DEF.) ***
twoway (connected risky_share_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(0(0.1)0.6, format(%9.1f))
gr export baseline_riskyshare.pdf, replace
* ------------------------------------

*** RISKY SHARE ACROSS WEALTH DIST (BASELINE DEF.) ***
twoway (connected risky_share_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)) ///
(connected cond_risky_share_p FW_percentile if ID_pctile == 1 & FW_percentile>20 , $connected_style sort(FW_percentile)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(0(0.1)0.6, format(%9.1f)) legend(order(1 ///
"Unconditional" 2 "Conditional") rows(2) region(fcolor(white)) position(11) ring(0))
gr export riskyshare_cond_uncond.pdf, replace
* ------------------------------------

*** RISKY SHARE PART. RATE ACROSS WEALTH DIST ***
twoway (connected part_Risky_assets_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)), ///
$general_style xtitle(Wealth Percentile) ytitle(Part. Rate) ylabel(, format(%9.1f))
gr export baseline_partrisky.pdf, replace
* ------------------------------------

*** RISKY SHARE FOR DIFFERENT DEFINITIONS ***
twoway (connected risky_share_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)) ///
(connected risky_share_RA_p FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)) ///
(connected portfolio_houseNH_p FW_percentile_NH if ID_pctile_NH == 1, $connected_style sort(FW_percentile_NH)) ///
(connected risky_share_NH_RA_p FW_percentile_NH_RA if ID_NH_RA == 1, $connected_style sort(FW_percentile_NH_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(, format(%9.1f)) legend(order(1 ///
"Baseline" 2 "Inc. RA" 3 "Inc. Housing" 4 "Inc. Both") rows(4) region(fcolor(white)) position(11) ring(0))
gr export robustness_share.pdf, replace 
* ------------------------------------

*** PART. RATE FOR DIFFERENT DEFINITIONS ***
twoway (connected part_Risky_assets_p FW_percentile if ID_pctile == 1, $connected_style sort(FW_percentile)) ///
(connected part_Risky_assets_RA_p FW_percentile_RA if ID_RA == 1, $connected_style sort(FW_percentile_RA)) ///
(connected part_Risky_assets_houseNH_p FW_percentile_NH if ID_pctile_NH == 1, $connected_style sort(FW_percentile_NH)) ///
(connected part_Risky_assets_houseNH_RA_p FW_percentile_NH_RA if ID_NH_RA == 1, $connected_style sort(FW_percentile_NH_RA)), ///
$general_style xtitle(Wealth Percentile) ytitle(Part. Rate) ylabel(, format(%9.1f)) legend(order(1 ///
"Baseline" 2 "Inc. RA" 3 "Inc. Housing" 4 "Inc. Both") rows(4) region(fcolor(white)) position(5) ring(0))
gr export robustness_part.pdf, replace 
* ------------------------------------

*** ASSETS ACROSS WEALTH DIST ***
twoway (connected house_share_p FW_percentile_NH_RA, $connected_style sort(FW_percentile_NH_RA)) ///
(connected stocks_plus_p FW_percentile_NH_RA, $connected_style sort(FW_percentile_NH_RA)) ///
(connected RA_share_p FW_percentile_NH_RA if RA_share_p<=1 , $connected_style sort(FW_percentile_NH_RA)) ///
if ID_NH_RA == 1, $general_style xtitle(Wealth Percentile) ytitle(Share) legend(order(1 ///
"Housing (net worth)" 2 "Risky" 3 "Ret. Accounts") rows(3) region(fcolor(white)) position(11) ring(0)) ylabel(, format(%9.1f))
gr export assets_share.pdf, replace
* ------------------------------------

*** OLD SHARES OUT OF RETIREMENT WEALTH ***
twoway (connected share_old_p FW_old_acc if ID_old == 1, $connected_style sort(FW_old_acc)), ///
$general_style xtitle(Wealth Percentile) ytitle(Risky Share) ylabel(0(0.1)0.6, format(%9.1f))
gr export old_riskyshare.pdf, replace
* ------------------------------------

****************************************
***             TABLES               ***
****************************************

cd "$tables" 

local assets Safe_assets checking_accounts savings_accounts money_market_accounts CDS bonds_safe mutual_safe ///
Risky_assets brokerage stocks mutual_risky bonds_risky total_IRA total_pension ///
net_house_wealth old_accounts financial_wealth

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
tex \textit{Total Safe Assets (S)} & `Safe_assets_mean' & `Safe_assets_sd' & `Safe_assets_p10' & `Safe_assets_p50' & `Safe_assets_p90' & `Safe_assets_part' \\ 
* checking accounts
tex Checking Account & `checking_accounts_mean' & `checking_accounts_sd' & `checking_accounts_p10' & `checking_accounts_p50' ///
			& `checking_accounts_p90' & `checking_accounts_part' \\
* savings accounts		
tex Savings Accounts & `savings_accounts_mean' & `savings_accounts_sd' & `savings_accounts_p10' & `savings_accounts_p50' /// 
			& `savings_accounts_p90' & `savings_accounts_part' \\
* money market accounts
tex Money Market Accounts & `money_market_accounts_mean' & `money_market_accounts_sd' & `money_market_accountsp10' ///
			& `money_market_accounts_p50' & `money_market_accounts_p90' & `money_market_accounts_part' \\
* cerificates of deposits
tex Cerificates of Deposits & `CDS_mean' & `CDS_sd' & `CDS_p10' & `CDS_p50' & `CDS_p90' & `CDS_part' \\
* savings_bond (safe)
tex Savings bond (safe) & `bonds_safe_mean' & `bonds_safe_sd' & `bonds_safe_p10' & `bonds_safe_p50' & `bonds_safe_p90' & `bonds_safe_part' \\ 
* mutual funds (safe)
tex Mutual Funds (safe) & `mutual_safe_mean' & `mutual_safe_sd' & `mutual_safe_p10' & `mutual_safe_p50' & `mutual_safe_p90' & `mutual_safe_part' \\ 
tex \\ 
* total risky
tex \textit{Total Risky Assets (R)} & `Risky_assets_mean' & `Risky_assets_sd' & `Risky_assets_p10' ///
		& `Risky_assets_p50' & `Risky_assets_p90' & `Risky_assets_part' \\ 
* brokerage
tex Brokerage & `brokerage_mean' & `brokerage_sd' & `brokerage_p10' & `brokerage_p50' & `brokerage_p90' & `brokerage_part' \\ 
* stocks
tex Stocks & `stocks_mean' & `stocks_sd' & `stocks_p10' & `stocks_p50' & `stocks_p90' & `stocks_part' \\
* mutual funds (risky)
tex Mutual Funds (risky) & `mutual_risky_mean' & `mutual_risky_sd' & `mutual_risky_p10' & `mutual_risky_p50' & `mutual_risky_p90' & `mutual_risky_part' \\ 
* savings bond (risky)
tex Savings bond (risky) & `bonds_risky_mean' & `bonds_risky_sd' & `bonds_risky_p10' & `bonds_risky_p50' ///
		& `bonds_risky_p90' & `bonds_risky_part' \\ 
tex \\
* baseline total FW
tex \textit{Baseline Financial Wealth (S + R)} & `financial_wealth_mean' & `financial_wealth_sd' & `financial_wealth_p10' ///
		& `financial_wealth_p50' & `financial_wealth_p90' & `financial_wealth_part' \\
tex \\
* total retirement accounst 
tex \textit{Retirement Wealth} & `old_accounts_mean' & `old_accounts_sd' & `old_accounts_p10' ///
		& `old_accounts_p50' & `old_accounts_p90' & `old_accounts_part' \\
* total IRA
tex IRA & `total_IRA_mean' & `total_IRA_sd' & `total_IRA_p10' & `total_IRA_p50' & `total_IRA_p90' & `total_IRA_part' \\
* total pensions 
tex Pensions & `total_pension_mean' & `total_pension_sd' & `total_pension_p10' & `total_pension_p50' & `total_pension_p90' & `total_pension_part' \\ 
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


/*
cd "$data"

outsheet financial_wealth Safe_assets Risky_assets using wealth.csv , comma nolabel

