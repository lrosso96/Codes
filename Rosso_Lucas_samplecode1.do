///////////////////Report #3: Holidays //////////////////////
//                  Date: 26-08-2020                       //
/////////////////////////////////////////////////////////////

* Sample code
* Author: Lucas Rosso
* Code extract from RA project with Rodrigo Wagner
* This code was originally separeted into a main code and other codes for better understanding

* Required packages 
*ssc install texdoc
*ssc install reghdfe
*ssc install binscatter

clear all
set more off, permanently
gl main "C:\Users\Lucas Rosso\Dropbox\Holidays and Growth" 
cd "${main}\Data"

import delimited "holidays_table_dm_csv.csv", varnames(1)
	 
* Create new date-related variables
gen year      = real(substr( countryyear, -4,4))
gen country   = substr(countryyear,1,  strlen(countryyear) - 4)
gen month_str = substr( date, -3,3)
gen day_num   = real(substr(date,1,  strlen(date) - 4))
label define months_label3  1 "Jan"  2 "Feb" 3 "Mar" 4 "Apr"  5 "May"  6 "Jun"  7 "Jul"  8 "Aug" 9 "Sep" 10 "Oct" 11 "Nov" 12 "Dec" 
encode month_str, generate(month) label(months_label3)
gen 		td = mdy(month, day_num,year)
format      td %td
gen        		dow = dow(td)
label define 		dow_label 0 "Sunday" 1 "Monday" 2 "Tuesday" 3 "Wednesday" 4 "Thursday" 5 "Friday" 6 "Saturday"
label values 	dow dow_label 

drop day countryyear v7 v8 month_str
ren  date date_dm_str
ren  day_num day

* Country codes using Raciborski, R. (2008). "kountry: A Stata utility for merging cross-country data from multiple sources," Stata Journal, 8(3), 390-400.
gen country2 = ustrregexra(country,"-"," ") // eliminates dash and uses space between words (e.g. dominican republic instead of dominican-republic)
kountry country2 , from(other) stuck 
ren _ISO  iso3_num
kountry   iso3_num , from(iso3n) to(iso3c)
ren _ISO  iso3c
*br country iso3*
		 
order td country dow day month year date  

*Merge with type of weekend data (flagged countries were left out)
sort  iso3c 
merge iso3c using "_temp_weekedays.dta", uniqusing
keep if _m==3
drop _m // keeps only the countries that were included (i.e. non flagged)

* Analysis of weekends
*gen elegible_type = inlist(type, "Public Holiday", "Public holiday", "Regular Holiday", "National holiday", "National Holiday", "Federal Holiday", "Federal Public Holiday")
g religious_holiday = 1 if inlist(type, "National holiday, Christian", "National holiday, Hebrew", "National holiday, Orthodox")
g elegible_type = inlist(type, "Public Holiday", "Public holiday", "Regular Holiday", "National holiday", "National Holiday", ///
"National holiday, Flag day", "Federal Holiday", "Federal Public Holiday")|religious_holiday==1
gen isweekend = 1 if inlist(dow,weekend_begins,weekend_ends)
g not_weekend = 1 if dow!=weekend_begins & dow!=weekend_ends
* below is the alternative analysis, assuming all countries follow MON-FRI workweek 
*  gen isweekend = inlist(dow, 0, 6)
keep if elegible_type

* Table for Appendix A1: Types of Holidays on our Sample
tab type if iso3c != "CHN"

collapse (sum) isweekend not_weekend religious_holiday, by(year iso3c)

/*
Download GDP growt data from World Bank
clear
wbopendata , indicator(NY.GDP.MKTP.KD.ZG , NY.GDP.PCAP.KN , NY.GDP.PCAP.PP.KD ) long
*/

***** MERGING DATA ****
* Merge WDI data
gen   countrycode = iso3c
merge year countrycode using "WDI_GDPgrowth.dta", unique sort
keep if _m==3
drop _m

* Merge IMF data
mer 1:1 iso3c year using "growth_data.dta"
keep if _m==3
drop _m

* Merge UNCTAD Data
mer 1:1 iso3c year using "UNCTAD_data.dta"

* Merge UNCTAD shares DATA
mer 1:1 iso3c year using UNCTAD_shares.dta, gen(merge2)

* Merge ToT data from PCTOT
mer 1:1 iso3c year using PCTOT.dta, gen(merge3)
drop if merge3==2

* Merge Annual Leave Data (internet)
mer 1:1 iso3c year using Annual_Leave_Data.dta, gen(merge4)
drop if merge4==2

* Merge Paid Holidays Data from Botero et. al.
mer m:1 iso3c using Regulation_of_labor.dta, gen(merge5)

* Merge Labour Share estimates Karabarbounis and Neiman
mer 1:1 iso3c year using "KN Labor Share Dataset.dta", gen(merge6)

* Merge Average Hours Worked from Penn World Tables 
mer 1:1 iso3c year using extract_pwt91.dta, gen(merge7)
drop if merge7==2

* Merge importance of family ties and leisure time from World Value Survey
*mer 1:1 iso3c year using WVS.dta, gen(merge8)
*drop if merge8==2

* Merge different labor regulation indices
mer m:1 iso3c using "Country Level Labor Regulation and Related circa 2007.dta", gen(merge9)
***********************************************************

* Preparing Panel
encode iso3c , gen(id)
tsset id year

** Dropping China **
drop if iso3c == "CHN" // See Appendix in document.

* Variables for (panel) regression
g log_rgdp_lc = ln(rgdp_lc) // Local currency (IMF data)
g log_rgdp = ln(ny_gdp_mktp_kd) // Usd (WB data)
g log_rgdp_pc = ln(ny_gdp_pcap_kd) // log real gdp per capita (in usd)
g l_rgdp_pc_ppp = ln(rgdp_pc_ppp) // log real gdp per capita (ppp)
g l_rgdp_ppp = ln(ny_gdp_mktp_pp_kd) // log real gdp (ppp)
g delta_gdp_dollar = d.log_rgdp
g delta_gdp_lc = d.log_rgdp_lc
g gdp_growth_lc = rgdp_lc_g/100
g working_days = 260 - not_weekend
g l_working_days = ln(working_days)
g l_not_weekend = ln(not_weekend)
g dl_working_days = d.l_working_days

g delta_gdp = d.log_rgdp *100
g delta_days = d.working_days

* Creating Sectorial Variables (in logs)
g l_wholesale = ln(wholesale)
g l_construction = ln(construction)
g l_mining = ln(mining-manufacturing) // Mining is reported in data as mining + manufacturing
g l_industry = ln(industry)
g l_services = ln(services)
g l_agriculture = ln(agriculture)
g l_manufacturing = ln(manufacturing)

* Demand variables
g l_hh_cons = ln(hh_cons) // HH consumption expenditure
g l_cons_exp = ln(cons_exp) // Final consumption expenditure
g l_gfcf = ln(gfcf) // gross capital fixed formation
g l_gov_cons = ln(gov_cons) // General government final consumption expenditure

g holidays = isweekend + not_weekend
g not_religious = holidays - religious_holiday

* For Figure_4.do
g extra_days = working_days - 250 
g d_log_rgdp_lc = d.log_rgdp_lc *100 // For Binscatter
*

* Change in ToT (%)
g tot_change = (d.CNPI/l.CNPI)*100
g l_tot = log(CNPI)
*

encode regioncode, gen(region_id)

* For Figure_6.do
bys iso3c: egen mean_holiday_1419 = mean(holidays) if inlist(year, 2014,2015,2016,2017,2018,2019)
by iso3c: egen mean_holiday_0005 = mean(holidays) if inlist(year, 2000,2001,2002,2003,2004,2005)

by iso3c: egen H_1419 = mean(mean_holiday_1419)
by iso3c: egen H_0005 = mean(mean_holiday_0005)
g diag_line = 25
*

* For Figure_7.do
bys iso3c: g gdp_19g = ((rgdp_lc[_n] - rgdp_lc[_n - 19])/rgdp_lc[_n - 19])*100 if year==2019
by iso3c: egen mean_holiday_0003 = mean(holidays) if inlist(year, 2000,2001,2002)
by iso3c: egen mean_holiday_1719 = mean(holidays) if inlist(year, 2017,2018,2019)
by iso3c: egen H_0003 = mean(mean_holiday_0003)
by iso3c: egen H_1719 = mean(mean_holiday_1719)
g holiday_19g = (H_1719/H_0003 -1)*100 if year==2019

bys iso3c: g CAGR_gdp = ((gdp_19g/100 +1)^(1/19)-1)*100
by iso3c: g CAGR_holiday = ((holiday_19g/100 +1)^(1/19)-1)*100
*

* For Regressions
g extra_day = -1*not_weekend 
replace unemployment_rate = unemployment_rate/100
replace rgdp_lc_g = rgdp_lc_g/100
g agriculture_share_v2 = agriculture_share/100
g industry_share_v2 = industry_share/100
g l_deflator = log(deflator)
egen y_bar = mean(log_rgdp_pc)
g y_tilda = log_rgdp_pc - y_bar
*

* For Table 2: Measures of Dispersion in Total Countries *
bys iso3c: egen aux1 = max(not_weekend) 
by iso3c: egen aux2 = min(not_weekend) 
g range_notweekend = aux1-aux2 //range of non-weekend holidays by country
by iso3c: egen sd_notweekend = sd(not_weekend) 
by iso3c: egen sd_lworking = sd(l_working_days) 
by iso3c: egen aux3 = max(l_working_days) 
by iso3c: egen aux4 = min(l_working_days) 
g range_lworking = aux3-aux4 // range of log working days

* Labor Share (in %)
g labor_share = TLS*100
egen labor_share_bar = mean(TLS)
g labor_share_tilda = TLS - labor_share_bar

sort  id year

* PWT Variables
g l_capital_stock = ln(rnna) //Log Capital Stock
g d_tfp = d.rtfpna/rtfpna // Percentage Change TFP
g l_tfp = ln(rtfpna) // Log TFP
g l_avh = ln(avh) // Log Avg. Hours Worked
g l_gdp_pwt = ln(rgdpna*1000000) // Log Real GDP at constant 2011 national prices (in mil. 2011US$)
g l_employment = ln(1000000*emp) // emp: employment (in millions)
*g l_labor = labsh*l_gdp_pwt // Labor outout
g aux_avh = d.l_avh // for table with identical number of obs.


* For Figure_11.do
g d_avh = d.avh

* For Country Time Trend
encode iso3c, gen(iso3c_num)

* Labor Regulation Variables
la var index_labor7a          "Employment Laws Index: Measures the protection of labor and employment laws
la var index_col_barg13       "Labor Union Power: Measures the statutory protection and power of unions"
la var index_col_disp13       "Collective Disputes: Measures the protection of workers during collective disputes"
la var index_industrial4a     "Collective Relation Laws Index: Measures the protection of collective relations laws"
la var cost_overtimen_a       "Cost of increasing the number of hours worked"
la var index_altern12         "Existence and cost of alternatives to the standard employment contract"
la var firing_cost_3years_n   "Cost of firing 20 percent of the firm’s workers (10% are fired for redundancy and 10% without cause)"
la var index_dism2            "Worker protection granted by law or mandatory collective agreements against dismissal"
* --------------------------

** Demeaned variables **
sum index_labor7a if year==2004 //fixing a year
g index_1 = index_labor7a - r(mean)
sum index_col_barg13 if year==2004
g index_2 = index_col_barg13 - r(mean)
sum index_col_disp13 if year==2004
g index_3 = index_col_disp13 - r(mean)
sum index_industrial4a if year==2004
g index_4 = index_industrial4a - r(mean)
sum bot_04 if year==2004
g index_5 = bot_04 - r(mean)
sum el2_04 if year==2004
g index_6 = el2_04 - r(mean)
sum HP if year==2004
g index_7 = HP - r(mean)
sum flex_hiring if year==2004
g index_8 = flex_hiring - r(mean)
sum cost_overtimen_a if year==2014
g index_9 = cost_overtimen_a - r(mean)
sum index_altern12 if year==2014
g index_10 = index_altern12 - r(mean)
sum firing_cost_3years_n if year==2014
g index_11 = firing_cost_3years_n - r(mean) 
sum index_dism2 if year==2014
g index_12 = index_dism2 - r(mean)
sum EMP_LAWS if year==2014
g index_13 = EMP_LAWS - r(mean)
* ----------------------------

cd "${main}\Figures"

***************************
*** No clustering in SE ***
***************************

** Industry Share
eststo clear
eststo: reghdfe d.log_rgdp_lc c.industry_share##c.l_working_days, abs(id year) 
qui margins, dydx(l_working_days) at(industry_share=(10(10)80)) 
marginsplot, plotopts(lcolor(%75) lpattern(dash) mcolor(%50) ///
msize(medsmall) mfcolor(%30)) recastci(rarea) level(90) ciopts(fcolor(%15) lcolor(%15)) ///
ytitle("Estimated Elasticity of Working Days") ///
yline(0, lwidth(thin) lpattern(solid) lcolor(red) noextend) ///
ylabel(, angle(horizontal) nogrid) ///
xtitle("Industry Share") graphregion(fcolor(white) ///
lcolor(white)) plotregion(fcolor(white) lcolor(white)) title("")
gr export het_by_industry.pdf, replace
**

** Agriculture Share
eststo clear
eststo: reghdfe d.log_rgdp_lc c.agriculture_share##c.l_working_days, abs(id year) 
margins, dydx(l_working_days) at(agriculture_share=(10(10)80)) 
marginsplot, plotopts(lcolor(%75) lpattern(dash) mcolor(%50) ///
msize(medsmall) mfcolor(%30)) recastci(rarea) level(90) ciopts(fcolor(%15) lcolor(%15)) ///
ytitle("Estimated Elasticity of Working Days") ///
yline(0, lwidth(thin) lpattern(solid) lcolor(red) noextend) ///
ylabel(, angle(horizontal) nogrid) ///
xtitle("Agriculture Share") graphregion(fcolor(white) ///
lcolor(white)) plotregion(fcolor(white) lcolor(white)) title("")
gr export het_by_agriculture.pdf, replace
**

* Binscatter for relation between GDP growth and Extra Days
binscatter d_log_rgdp_lc extra_days, control(i.year) absorb(id) ///
ytitle(Real GDP growth (%)) ylabel(, ///
angle(horizontal) nogrid) xtitle("Extra Days") ///
mc(blue%50) legend(off) graphregion(fcolor(white) lcolor(white))
gr export figure_4.pdf, replace

cd "${main}\Tables"


*********************************************************
** Regression 1.- Growth v/s holidays,   Exc 2008-2009 **
*********************************************************

*********************
*** Days in Level ***
*********************

texdoc i "New_reg1_exc0809.tex", replace force

* (1) *
reghdfe d.log_rgdp_lc working_days                           if inlist(year,2008,2009)==0,    abs(id year) cluster(id ) 

local working_days1 = _b[working_days]
local working_days1: di %9.4f `working_days1'

qui test _b[working_days] = 0
global p_working_days1: di %12.3fc r(p)
global star_working_days1 = cond(${p_working_days1}<.01,"***",cond(${p_working_days1}<.05,"**",cond(${p_working_days1}<.1,"*","")))

matrix var_cov1 = e(V)
local se_working_days1 = sqrt(var_cov1[1,1])
local se_working_days1: di %9.4f `se_working_days1'  

local obs1 = e(N)
local obs1: di %9.0f `obs1'

local countries1 = e(N_clust) 
local countries1: di %9.0f `countries1'

local r2_a1 = e(r2_a) 
local r2_a1: di %9.2f `r2_a1'

test _b[working_days] = 1/250
local test1_1 = r(p)
local test1_1: di %9.3f `test1_1'

local individual1: di "\checkmark"
local year1: di "\checkmark"
* ------------------------------------------

* (2) *
reghdfe d.log_rgdp_lc working_days l_tot                     if inlist(year,2008,2009)==0,    abs(id year) cluster(id ) 

local working_days2 = _b[working_days]
local working_days2: di %9.4f `working_days2'

qui test _b[working_days] = 0
global p_working_days2: di %12.3fc r(p)
global star_working_days2 = cond(${p_working_days2}<.01,"***",cond(${p_working_days2}<.05,"**",cond(${p_working_days2}<.1,"*","")))

matrix var_cov2 = e(V)
local se_working_days2 = sqrt(var_cov2[1,1])
local se_working_days2: di %9.4f `se_working_days2'

local l_tot1 = _b[l_tot]
local l_tot1: di %9.2f `l_tot1'

qui test _b[l_tot] = 0
global p_l_tot1: di %12.3fc r(p)
global star_l_tot1 = cond(${p_l_tot1}<.01,"***",cond(${p_l_tot1}<.05,"**",cond(${p_l_tot1}<.1,"*","")))

local se_l_tot1 = sqrt(var_cov2[2,2])
local se_l_tot1: di %9.2f `se_l_tot1'

local obs2 = e(N)
local obs2: di %9.0f `obs2'

local countries2 = e(N_clust) 
local countries2: di %9.0f `countries2'

local r2_a2 = e(r2_a) 
local r2_a2: di %9.2f `r2_a2'

test _b[working_days] = 1/250
local test1_2 = r(p)
local test1_2: di %9.3f `test1_2'

local individual2: di "\checkmark"
local year2: di "\checkmark"
* ------------------------------------------


********************
*** Days in Logs ***
********************

* (3) *
reghdfe d.log_rgdp_lc l_working_days                           if inlist(year,2008,2009)==0,    abs(id year) cluster(id ) 

local l_working_days1 = _b[l_working_days]
local l_working_days1: di %9.2f `l_working_days1'

qui test _b[l_working_days] = 0
global p_l_working_days1: di %12.3fc r(p)
global star_l_working_days1 = cond(${p_l_working_days1}<.01,"***",cond(${p_l_working_days1}<.05,"**",cond(${p_l_working_days1}<.1,"*","")))

matrix l_var_cov1 = e(V)
local se_l_working_days1 = sqrt(l_var_cov1[1,1])
local se_l_working_days1: di %9.2f `se_l_working_days1'  

local l_obs1 = e(N)
local l_obs1: di %9.0f `l_obs1'

local l_countries1 = e(N_clust) 
local l_countries1: di %9.0f `l_countries1'

local l_r2_a1 = e(r2_a) 
local l_r2_a1: di %9.2f `l_r2_a1'

test _b[l_working_days] = 1
local l_test1_1 = r(p)
local l_test1_1: di %9.3f `l_test1_1'

test _b[l_working_days] = 0.4
local l_test2_1 = r(p)
local l_test2_1: di %9.3f `l_test2_1'

local l_individual1: di "\checkmark"
local l_year1: di "\checkmark"
* ------------------------------------------

* (4) *
reghdfe d.log_rgdp_lc l_working_days l_tot                     if inlist(year,2008,2009)==0,    abs(id year) cluster(id ) 

local l_working_days2 = _b[l_working_days]
local l_working_days2: di %9.2f `l_working_days2'

qui test _b[l_working_days] = 0
global p_l_working_days2: di %12.3fc r(p)
global star_l_working_days2 = cond(${p_l_working_days2}<.01,"***",cond(${p_l_working_days2}<.05,"**",cond(${p_l_working_days2}<.1,"*","")))


matrix l_var_cov2 = e(V)
local se_l_working_days2 = sqrt(l_var_cov2[1,1])
local se_l_working_days2: di %9.2f `se_l_working_days2'

local l_tot2 = _b[l_tot]
local l_tot2: di %9.2f `l_tot2'

local se_l_tot2 = sqrt(l_var_cov2[2,2])
local se_l_tot2: di %9.2f `se_l_tot2'

qui test _b[l_tot] = 0
global p_l_tot2: di %12.3fc r(p)
global star_l_tot2 = cond(${p_l_tot2}<.01,"***",cond(${p_l_tot2}<.05,"**",cond(${p_l_tot2}<.1,"*","")))

local l_obs2 = e(N)
local l_obs2: di %9.0f `l_obs2'

local l_countries2 = e(N_clust) 
local l_countries2: di %9.0f `l_countries1'

local l_r2_a2 = e(r2_a) 
local l_r2_a2: di %9.2f `l_r2_a2'

test _b[l_working_days] = 1
local l_test1_2 = r(p)
local l_test1_2: di %9.3f `l_test1_2'

test _b[l_working_days] = 0.4
local l_test2_2 = r(p)
local l_test2_2: di %9.3f `l_test2_2'

local l_individual2: di "\checkmark"
local l_year2: di "\checkmark"
* ------------------------------------------


***  CREATING LATEX TABLE ***
tex \begin{tabular}{l*{2}{>{\centering\arraybackslash}m{2cm}}*{1}c*{2}{>{\centering\arraybackslash}m{2cm}}} \hline \hline
tex  & \multicolumn{5}{c}{Dependent Variable: GDP Growth} \\ \cline{2-6}
tex  & \multicolumn{2}{c}{One Day Change} & & \multicolumn{2}{c}{Elasticity} \\ \cline{2-3} \cline{5-6}  
tex  & (1) & (2) & & (3) & (4) \\ \midrule 
tex Extra Days (level or log)                   & `working_days1'$star_working_days1 &  `working_days2'$star_working_days2 & & ///
                                                 `l_working_days1'$star_l_working_days1 &  `l_working_days2'$star_l_working_days2  \\
tex                                             & \RemoveSpaces{[`se_working_days1']}  & \RemoveSpaces{[`se_working_days2']} & & ///
	\RemoveSpaces{[`se_l_working_days1']}  & \RemoveSpaces{[`se_l_working_days2']}  \\ 
tex  Log Terms of Trade                         &                 &  `l_tot1'$star_l_tot1 & &                 &  `l_tot2'$star_l_tot2  \\ 
tex                                             &                 &  \RemoveSpaces{[`se_l_tot1']} & &  &  \RemoveSpaces{[`se_l_tot2']} \\ 
tex                                             &                 &          & &                 &                \\ 
tex Country Fixed Effect                        & `individual1'   & `individual2' & & `l_individual1'   & `l_individual2' \\ 
tex Year Fixed Effect                           & `year1'         & `year2'  & & `l_year1'         & `l_year2'          \\ 
tex Observations                                & `obs1'          & `obs2' & & `l_obs1'          & `l_obs2'     \\ 
tex Countries                                   & `countries1'    & `countries2' & & `l_countries1'    & `l_countries2' \\ 
tex Adj. R^{2}                                  & `r2_a1'         & `r2_a2' & & `l_r2_a1'         & `l_r2_a2'     \\ 
tex F test Extra Days = 1/250                   &  \RemoveSpaces{(`test1_1')}      & \RemoveSpaces{(`test1_2')} & &      &         \\ 
tex F test Labor Elast. = 1                     &                 &           & &  \RemoveSpaces{(`l_test1_1')}  & \RemoveSpaces{(`l_test1_2')}  \\ 
tex F test Labor Elast. = 0.4                     &            &           & &  \RemoveSpaces{(`l_test2_1')} & \RemoveSpaces{(`l_test2_2')}     \\ 
tex  \bottomrule
tex \end{tabular}
texdoc close 
* --------------------------------------------------------------------------------
