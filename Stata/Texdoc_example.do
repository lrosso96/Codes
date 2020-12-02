***************************************************
***             HOLIDAYS AND GROWTH             ***
***************************************************

cd "${main}\Tables"

* must add in latex preamble to remove space between () and standard errors 
* when importing tables from stata
*\newcommand{\RemoveSpaces}[1]{%
*  \begingroup
*  \spaceskip=2sp
*  \xspaceskip=1sp
*  #1%
*  \endgroup}

* to add math symbols use \$ otherwise stata will interpret $ as global.

*ssc install texdoc

**********************************************************
**   Regression 1.- Growth v/s holidays, Full Sample    **
**********************************************************

*********************
*** Days in Level ***
*********************

texdoc i "New_reg1_fullsample.tex", replace force

* (1) *
reghdfe d.log_rgdp_lc working_days,    abs(id year) cluster(id ) 

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
reghdfe d.log_rgdp_lc working_days l_tot,    abs(id year) cluster(id ) 

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
local countries2: di %9.0f `countries1'

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
reghdfe d.log_rgdp_lc l_working_days,    abs(id year) cluster(id ) 

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
reghdfe d.log_rgdp_lc l_working_days l_tot,    abs(id year) cluster(id ) 

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
tex  & \multicolumn{5}{c}{GDP Growth} \\ \cline{2-6}
tex  & \multicolumn{2}{c}{Extra Days in Level} & & \multicolumn{2}{c}{Extra Days in Logs} \\ \cline{2-3} \cline{5-6}  
tex  & (1) & (2) & & (3) & (4) \\ \midrule 
tex Extra Days                                  & `working_days1'$star_working_days1 &  `working_days2'$star_working_days2 & & ///
                                                 `l_working_days1'$star_l_working_days1 &  `l_working_days2'$star_l_working_days2  \\
tex                                             & \RemoveSpaces{[`se_working_days1']}  & \RemoveSpaces{[`se_working_days2']} ///
	& & \RemoveSpaces{[`se_l_working_days1']}  & \RemoveSpaces{[`se_l_working_days2']}  \\ 
tex  Log Terms of Trade                         &                 &  `l_tot1'$star_l_tot1 & &                 &  `l_tot2'$star_l_tot2  \\ 
tex                                             &                 &  \RemoveSpaces{[`se_l_tot1']} & &                 &  \RemoveSpaces{[`se_l_tot2']} \\ 
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
