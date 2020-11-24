/****************************************************************************************************
* Title: A Heterogeneous Agent Model of Portfolio Choice and its Implications for Wealth Inequality
* Created by: Lucas Rosso
* Created on: 19/11/2020
* Purpose: This dofile merges data from the SCF (1998-2019) and creates relevant variables for
			analysis and calibration.
* Last Modified on: 19/11/2020
* Last Modified by: LR
* Edits:
	[19/11/2020]: Created dofile
	[23/11/2020]: Added Surveys of 2010, 2013, 2016 and 2019 + CPI data
*****************************************************************************************************/

clear all
set more off, permanently
set maxvar 30000
gl main "C:\Users\Lucas Rosso\Desktop\ME\Tesis\TESIS\Codigos\Stata" // this must change
cd "${main}\Data\Raw"

// variables names
* 1998: https://www.federalreserve.gov/econres/files/codebk98.txt
* 2001: https://www.federalreserve.gov/econres/files/codebk2001.txt
* 2004: https://www.federalreserve.gov/econres/files/codebk2004.txt
* 2007: https://www.federalreserve.gov/econres/files/codebk2007.txt
* 2010: https://www.federalreserve.gov/econres/files/codebk2010.txt
* 2013: https://www.federalreserve.gov/econres/files/codebk2013.txt
* 2016: https://www.federalreserve.gov/econres/files/codebk2016.txt
* 2019: https://www.federalreserve.gov/econres/files/codebk2019.txt
* -------------------------------

**********************************
*********  1998 Edition  ********* 
**********************************

use p98i6

g year = 1998
gen wgt = X42001

*****************************************
*********    DEMOGRAPHICS  **************
*****************************************
  
gen male=1 if X8021==1
replace male=0 if (male==.)

gen person1_age = X8022
gen person2_age = X104
gen person3_age = X110
gen person4_age = X116
gen person5_age = X122
gen person6_age = X128
gen age= person1_age

gen educ= (X5904==1) // college degree

gen income= max(X5729,0) // (yearly)
gen children=X5910
replace children=0 if children==-1
gen marital=1 if  (X7372==1) 
replace marital=0 if marital==.
* --------------------------------------

******************************************
**************    ASSETS    ************** 
******************************************

*** Checking Accounts ***
gen checking_accounts1 = max(X3506,0) if (X3507==5)
replace  checking_accounts1=0  if ( checking_accounts1 ==.  )
gen checking_accounts2 = max(X3510,0) if (X3511==5)
replace  checking_accounts2=0  if ( checking_accounts2 ==.  )
gen checking_accounts3 = max(X3514,0) if (X3515==5)
replace  checking_accounts3=0  if ( checking_accounts3 ==.  )
gen checking_accounts4 = max(X3518,0) if (X3519==5)
replace  checking_accounts4=0  if ( checking_accounts4 ==.  )
gen checking_accounts5 = max(X3522,0) if (X3523==5)
replace  checking_accounts5=0  if ( checking_accounts5 ==.  )
gen checking_accounts6 = max(X3526,0) if (X3527==5)
replace  checking_accounts6=0  if ( checking_accounts6 ==.  )
 
 
gen checking_accounts = checking_accounts1 +  checking_accounts2  + checking_accounts3 +  checking_accounts4   ///
                        + checking_accounts5   +  checking_accounts6  
drop checking_accounts1  checking_accounts2  checking_accounts3  checking_accounts4 checking_accounts5  checking_accounts6  
gen part_checking_accounts = ( checking_accounts>0)
* --------------------------------------

*** Savings Accounts ***
gen savings_accounts = max(X3804,0) + max(X3807,0) + max(X3810,0)  + max(X3813,0)  + max(X3816,0) 
gen part_savings_accounts = ( savings_accounts>0)
* --------------------------------------

*** Money Market Accounts ***
gen MMA1 = max(X3506,0) if (X3507==1) // money market accounts reported on checkings account section (X3507==1) 
replace  MMA1=0  if ( MMA1 ==.  )
gen  MMA2 = max(X3510,0) if (X3511==1)
replace  MMA2=0  if ( MMA2 ==.  )
gen MMA3 = max(X3514,0) if (X3515==1)
replace  MMA3=0  if ( MMA3 ==.  )
gen MMA4 = max(X3518,0) if (X3519==1)
replace  MMA4=0  if ( MMA4 ==.  )
gen MMA5 = max(X3522,0) if (X3523==1)
replace  MMA5=0  if ( MMA5 ==.  )
gen MMA6 = max(X3526,0) if (X3527==1)
replace  MMA6=0  if ( MMA6 ==.  )
gen  MMA7= max(0,X3706) + max(0,X3711) + max(0,X3716) // money market accounts reported on the correct section.
gen money_market_accounts= MMA1+ MMA2 + MMA3 + MMA4 + MMA5 + MMA6 + MMA7
drop MMA1-MMA7
gen part_money_market_accounts= ( money_market_accounts>0)
* --------------------------------------

*** Cerificates of Deposits ***
gen CDS =max(0,X3721) if (X7620<=3) // X7620: 1. Joint account 2. R's account 3. Spouse's/partner's account
replace CDS=0 if (CDS==.)
gen part_CDS = ( CDS>0 )
* --------------------------------------

*** Mutual funds ***
* Stock funds  
gen  mutual_stock_funds= max(0,X3822) if ( X3821==1)
replace mutual_stock_funds =0 if  (  mutual_stock_funds ==. ) 
* tax free bonds
gen  mutual_taxfree_bfunds= max(0,X3824) if ( X3823==1)
replace mutual_taxfree_bfunds =0 if  (  mutual_taxfree_bfunds ==. ) 
* governemnt backed bonds
gen  mutual_gov_bfunds= max(0,X3826) if ( X3825==1)
replace mutual_gov_bfunds =0 if  (  mutual_gov_bfunds ==. ) 
* other bond funds
gen  mutual_other_bfunds= max(0,X3828) if ( X3827==1)
replace mutual_other_bfunds =0 if  (  mutual_other_bfunds ==. ) 
*combination funds
gen  mutual_comb_funds= max(0,X3830) if ( X3829==1)
replace mutual_comb_funds =0 if  (  mutual_comb_funds ==. ) 

gen mutual_safe = mutual_taxfree_bfunds + mutual_gov_bfunds +  mutual_other_bfunds + (mutual_comb_funds)/2
gen mutual_risky = mutual_stock_funds + (mutual_comb_funds)/2 

gen part_mutual_safe = (mutual_safe > 0)
gen part_mutual_risky = (mutual_risky > 0)
* --------------------------------------

*** Bonds ***
gen savings_bonds_safe = X3902 +  X3910 + X3908   // face value               
gen savings_bonds_risky = X3906 + X7634 + X7633   

gen part_bonds_safe = (savings_bonds_safe)
gen part_bonds_risky = (savings_bonds_risky)
* --------------------------------------

*** Life Insurance & Other assets ***
gen life_insurance = max(0,X4006 - X4010)
gen misc_assets = max(0,X4018 + X4022 + X4026 +  X4030   - X4032) 

gen part_life_insurance = (life_insurance>0)
gen part_misc_assets    = (misc_assets>0)
* --------------------------------------

*** Publicly traded Stocks ***
gen  stocks1= max(0,X3915) if ( X3913==1)
replace stocks1 =0 if  (  stocks1 ==. )
gen stocks2= max(0,X7641) if (X7192==3)
replace stocks2=0 if (stocks2==.)

gen stocks =stocks1+stocks2 
gen part_stocks = (stocks >0)
drop stocks1-stocks2
* --------------------------------------

*** Individual Retirement Accounts ***
/*
           1.  *CDs/Bank accounts; "money market"
           2.  *Stock; "mutual funds"
           3.  *Bonds/Similar assets; T-Bills; treasury notes
           4.  Combinations of 1, 2, & 3 ; "mixed"/"diversified" -- NFS
           5.  Combination of 2 & 3 above
           6.  Combination of 1 & 2 above
           11.  Universal life policy or other similar insurance product
           12.  Annuity
           13.  Commodities
           14.  Real estate/mortgages
           15.  Limited partnership/other similar investment
           16.  Brokerage account/cash management account (CMA)
           -7.  *Split/Other
*/

gen IRA_safe = (X3610 + X3620 + X3630)  if inlist(X3631, 1,3,11)     
replace IRA_safe=0 if IRA_safe==. 
replace IRA_safe=  2*( X3610 + X3620 + X3630 )/3 if X3631==4
* if ==4 then it is a combination of bank accounts, stocks and bonds so you have to give 2/3 to safe
replace IRA_safe=  ( X3610 + X3620 + X3630 ) /2 if inlist(X3631, 5,6,-7)                
* give half to safe if both and half to risky

gen IRA_risky = (X3610 + X3620 + X3630)  if inlist(X3631, 2,14,15,16)  
* stocks  or investments or brokerage accounts
replace IRA_risky=0 if IRA_risky==. 
replace IRA_risky= ( X3610 + X3620 + X3630 )/3 if X3631==4
replace IRA_risky=   ( X3610 + X3620 + X3630 )/2 if inlist(X3631, 5,6,-7) 

gen part_IRA_safe  = (IRA_safe>0)
gen part_IRA_risky = (IRA_risky>0)
* --------------------------------------

*** Brokerage Accounts ***
gen brokerage = max(0,X3930-X3932)
gen part_brokerage = (brokerage>0)
* --------------------------------------

*** Annuities ***
/*
           1.  *STOCKS; MUTUAL FUND
           2.  *BONDS/INTEREST; CDS/MONEY MARKET
           3.  REAL ESTATE
           5.  *SPLIT BETWEEN STOCKS/INTEREST; COMBINATION OF 1 & 2; MUTUAL FUNDS AND CD'S
           6.  MIXED OR DIVERSIFIED
           7.  *LIFE INSURANCE/FIXED CONTRACT; ANNUITIES
           8.  TANGIBLE ASSETS OTHER THAN REAL ESTATE
           9.  INTANGIBLE ASSETS, N.E.C.
           -7.  *OTHER
*/

gen annuity_safe = (X6820)  if inlist(X6826, 2,7)      
replace annuity_safe=0 if annuity_safe==. 
replace annuity_safe=  (X6820)/2 if inlist(X6826, 5,6,-7)        
* give half to safe if both and half to risky

gen annuity_risky = (X6820)   if inlist(X6826, 1,3)   
* stocks or real estate
replace annuity_risky=0 if annuity_risky==. 
replace annuity_risky=  (X6820) /2 if inlist(X6826, 5,6,-7) 

gen part_annuity_safe = (annuity_safe>0)
gen part_annuity_risky = (annuity_risky>0)
* --------------------------------------

*** Trusts ***
/*
           1.  *STOCKS; MUTUAL FUND
           2.  *BONDS/INTEREST/CDS/MONEY MARKET
           3.  REAL ESTATE
           5.  *COMBINATION OF 1 & 2; MUTUAL FUNDS AND CD'S
           6.  MIXED OR DIVERSIFIED
           7.  *LIFE INSURANCE/FIXED CONTRACT; ANNUITIES
           8.  TANGIBLE ASSETS OTHER THAN REAL ESTATE
           9.  INTANGIBLE ASSETS, N.E.C.
          -7.  *OTHER
*/
gen trust_safe = (X6835)  if inlist(X6841, 2,7)    
* bonds or life insurance
replace trust_safe=0 if trust_safe==. 
replace trust_safe= (X6835)/2 if inlist(X6841, 5,6,-7)       
* give half to safe if both and half to risky
gen trust_risky = (X6835)   if inlist(X6841, 1,3)   
* stocks or real estate
replace trust_risky=0 if trust_risky==. 
replace trust_risky= (X6835) /2 if inlist(X6841, 5,6,-7)

gen part_trust_safe  = (trust_safe>0)
gen part_trust_risky = (trust_risky>0)
* --------------------------------------

*** Housing ***       
gen farm_value  = X513 + X526
gen mobile_site = X604 + X614 + X623 
gen value_home  = X716  

gen house_mortgage= (X805 + X905 + X1005) 
gen other_loans= X1044

gen lines_of_Credit =  X1108 + X1119 + X1130 + X1136 
replace lines_of_Credit = 0 if  lines_of_Credit<0 

gen land_contract_lend = (X1409 + X1509 + X1609  + X1619)   
gen land_contract_borrow = X1417 + X1517  + X1617 + X1621 
gen real_estate =  X1706*(X1705/10000)  + X1806*(X1805/10000)  + X1906*( X1905/10000 ) + X2002 + X2012       
gen borrowing_housing =  X1715 *(X1705/10000) + X1815*(X1805/10000)  + X1915*(X1905/10000) + X2006 + X2016


gen house_wealth = ( farm_value + mobile_site + value_home )  +  land_contract_lend  +  real_estate
gen house_wealth2= ( farm_value + mobile_site + value_home )  +  real_estate

gen mortgages=  house_mortgage + other_loans +  land_contract_borrow + borrowing_housing
gen mortgages2= house_mortgage + other_loans + borrowing_housing

gen net_house_wealth = house_wealth-mortgages
gen net_house_wealth2 = house_wealth2-mortgages2

gen owner=0 if net_house_wealth==0
replace owner=1 if owner==.
* --------------------------------------

*** Pensions ***
*Did you tel me about this (pension) loan earlier?
* If yes I don't count it 
replace X4229=X4229 if X4230==5     
replace X4229=0 if X4229==. 
replace X4329=X4329 if X4330==5 
replace X4329=0 if X4329==. 
replace X4429=X4429 if X4430==5 
replace X4429=0 if X4429==. 
replace X4829=X4829 if X4830==5 
replace X4829=0 if X4829==. 
replace X4929=X4929 if X4930==5 
replace X4929=0 if X4929==. 
replace X5029=X5029 if X5030==5 
replace X5029=0 if X5029==. 
* Now create the pension accounts net of loans  
gen net_pension1=  max(0,X4226 - X4229)
gen net_pension2=  max(0,X4326 - X4329)
gen net_pension3=  max(0,X4426 - X4429)
gen net_pension4=  max(0,X4826 - X4829)
gen net_pension5=  max(0,X4926 - X4929)
gen net_pension6=  max(0,X5026 - X5029)

* Then break them up in safe/risky 
gen weight_1=0.5

/*
			1.  *Thrift or Savings
			2.  *401K/403B/SRA
            3.  *Profit Sharing
            4.  *Stock purchase/ESOP
            7.  Deferred compensation
            11.  IRA-SEP (not to be confused with a regular IRA)
            12.  Defined-contribution plan; TIAA-CREF (Teachers
                 Insurance and Annuity Association/College Retirement
                 Equity Fund)
            13.  Money purchase plan
            14.  Tax-deferred annuity (TDA); tax-sheltered annuity (TSA)
            17.  Other type of annuity (include ERISA plans here
                 unless otherwise specified)
            18.  Other salary reduction plan; deferred compensation plan
            24.  Other state/local government plan
            25.  Other federal government plan
            26.  Other type of account
            -7.  *Other; combination
*/

*safe pensions  (I consider 401K both risky and non-risky) 
gen k1=weight_1*net_pension1 if X4216==1  | X4216==2  | X4216==12 |  X4216==13   |  X4216==14  |  X4216==17  | X4216==-7                     
replace k1=0 if k1==. 
gen k2=weight_1*net_pension2 if X4316==1  | X4316==2  | X4316==12  |  X4316==13  |  X4316==14  |  X4316==17 | X4316 ==-7
replace k2=0 if k2==. 
gen k3=weight_1*net_pension3 if X4416==1  | X4416==2  | X4416==12  |  X4416==13  |  X4416==14  |  X4416==17  | X4416 ==-7              
replace k3=0 if k3==. 
gen k4=weight_1*net_pension4 if X4816==1  | X4816==2  | X4816==12  |  X4816==13  |  X4816==14  |  X4816==17 | X4816 ==-7               
replace k4=0 if k4==. 
gen k5=weight_1*net_pension5 if X4916==1  | X4916==2   | X4916==12 |  X4916==13  |  X4916==14  |  X4916==17 | X4916 ==-7                
replace k5=0 if k5==. 
gen k6=weight_1*net_pension6 if X5016==1  | X5016==2   | X5016==12 |  X5016==13  |  X5016==14  |  X5016==17 | X5016 ==-7              
replace k6=0 if k6==. 
gen pensions_safe =  k1 + k2 + k3 + k4 + k5 +k6          


*Risky pensions 
gen k7=net_pension1 if  X4216==3|X4216==4                
replace k7=(1-weight_1)*net_pension1 if inlist(X4216, 1,2,12,13,14,17,-7) 
replace k7=0 if k7==.   

gen k8=net_pension2 if  X4316==3|X4316==4                   
replace k8=(1-weight_1)*net_pension2 if inlist(X4316, 1,2,12,13,14,17,-7)
replace k8=0 if k8==. 
  
gen k9=net_pension3 if   X4416==3|X4416==4                  
replace k9=(1-weight_1)*net_pension3 if inlist(X4416, 1,2,12,13,14,17,-7)             
replace k9=0 if k9==. 


gen k10=net_pension4 if X4816==3|X4816==4                    
replace k10=(1-weight_1)*net_pension4 if inlist(X4816, 1,2,12,13,14,17,-7)                   
replace k10=0 if k10==. 
  
gen k11=net_pension5 if X4916==3|X4916==4 
replace k11=(1-weight_1)*net_pension5 if inlist(X4916, 1,2,12,13,14,17,-7)                                 
replace k11=0 if k11==. 
  
gen k12=net_pension6 if X5016==3|X5016==4                  
replace k12=(1-weight_1)*net_pension6 if inlist(X5016, 1,2,12,13,14,17,-7)                
replace k12=0 if k12==. 

gen pensions_risky =  k7 + k8 + k9 + k10 + k11 +k12          

gen other_pension_plans_safe =      weight_1*max(0,X5604 + X5612 + X5620 + X5628 + X5636 + X5644 )     
gen other_pension_plans_risky = (1-weight_1)*max(0,X5604 + X5612 + X5620 + X5628 + X5636 + X5644 ) 

gen total_pension_safe  = pensions_safe + other_pension_plans_safe
gen total_pension_risky = pensions_risky + other_pension_plans_risky

gen part_pension_safe  = (total_pension_safe>0)
gen part_pension_risky = (total_pension_risky>0)
* --------------------------------------


***********************************************
*************** Other Variables *************** 
***********************************************

*** Educational Debt ***
gen education_loan= max(0,X7824) + max(0,X7847) + max(0,X7870) +max(0,X7924) + max(0,X7947)  + max(0,X7970)  
gen part_education_loan = (education_loan>0)
* --------------------------------------

*** Other Consumer Loans ***
gen other1= X2723 if X6842==5 // X6842 - X6847: Is this loan one that you told me about when we talked about your business?
replace other1=0  if other1==.
gen other2= X2740 if X6843==5
replace other2=0  if other2==.
gen other3= X2823 if X6844==5
replace other3=0  if other3==.
gen other4= X2840 if X6845==5
replace other4=0  if other4==.
gen other5= X2923 if X6846==5
replace other5=0  if other5==.
gen other6= X2940 if X6847==5
replace other6=0  if other6==.

gen other_consumer_loans= other1 + other2  +  other3 +  other4 +  other5 +  other6
drop other1-other5
gen part_other_consumer_loans = (other_consumer_loans>0)

*** Credit Card Debt ***
gen credit_cards_debt= max(0,X413) + max(0,X421) + max(0,X424) +max(0,X427) + max(0,X430)  + max(X7575,0) 
gen part_credit_cards_debt = (credit_cards_debt>0)

*** Interest Rates in Consumer Loans ***
/*
gen rate_c1 = X2724 if X2724>0 // anual rate
gen rate_c2 = X2741 if X2741>0
gen rate_c3 = X2824 if X2824>0
gen rate_c4 = X2841 if X2841>0
gen rate_c5 = X2924 if X2924>0
gen rate_c6 = X2941 if X2941>0

egen rate_other_consumer_loans = rowmean(rate_c1 rate_c2 rate_c3 rate_c4 rate_c5 rate_c6) // mean anual rate, conditional on having debt.
*/
* --------------------------------------


**********************************************************
**************        Definitions       ******************
**********************************************************

*** SAFE ASSETS ***
* benchmark definition for safe assets
gen Safe_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe  ///
		       + IRA_safe + pensions_safe  + other_pension_plans_safe ) 

* safe assets net of debt (not including educational loans)
gen net_Safe_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe ///
			   + IRA_safe + pensions_safe  +  other_pension_plans_safe )  - credit_cards_debt  - other_consumer_loans 

* safe assets net of debt (including educational loans)
gen net_Safe_assets2 = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe  /// 
			   + IRA_safe + pensions_safe  +  other_pension_plans_safe )  - credit_cards_debt  - other_consumer_loans -education_loan
* --------------------------------------

*** LIQUID ASSETS ***
* benchmark definition for liquid assets
gen Liquid_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)

* liquid assets net of debt (not including educational loans)
gen net_Liquid_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)  /// 
			   - credit_cards_debt  - other_consumer_loans			   

* liquid assets net of debt (including educational loans)
gen net_Liquid_assets2 = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)  /// 
			   - credit_cards_debt  - other_consumer_loans -education_loan			   
* --------------------------------------			   

*** RISKY ASSETS ***			  
* benchmark definition for risky assets
gen Risky_assets          = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                + IRA_risky  + pensions_risky  +  other_pension_plans_risky  )

* risky assets including housing net worth				
gen Risky_assets_houseNH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                  + IRA_risky  + pensions_risky  + other_pension_plans_risky  +  net_house_wealth)

* risky assets including housing gross worth				  
gen Risky_assets_houseH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky  ///
                  + IRA_risky  + pensions_risky  + other_pension_plans_risky  +  house_wealth)
* --------------------------------------			   

*** ILIQUID ASSETS ***
* benchmark definition for illiquid assets
gen Illiquid_assets          = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                + IRA_risky  + pensions_risky  +  other_pension_plans_risky + IRA_safe + pensions_safe  + other_pension_plans_safe  )

* illiquid assets including housing net worth				
gen Illiquid_assets_houseNH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                  + IRA_risky  + pensions_risky  + other_pension_plans_risky + IRA_safe + pensions_safe  + other_pension_plans_safe  +  net_house_wealth)

* illiquid assets including housing gross worth				  
gen Illiquid_assets_houseH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky  ///
                  + IRA_risky  + pensions_risky  + other_pension_plans_risky + IRA_safe + pensions_safe  + other_pension_plans_safe  +  house_wealth)
* --------------------------------------

*** FINANCIAL WEALTH ***		  
gen financial_wealth = Safe_assets + Risky_assets // equivalent to liquid + illiquid
gen financial_wealth_houseH = Safe_assets + Risky_assets_houseH 
gen financial_wealth_houseNH = Safe_assets + Risky_assets_houseNH 
gen financial_wealth_debt = net_Safe_assets + Risky_assets 
gen financial_wealth_debt2 = net_Safe_assets2 + Risky_assets 
* --------------------------------------

*** PORTFOLIO CHOICE ***
* risky/safe
gen risky_share         = max(0,Risky_assets/financial_wealth)  
gen portfolio_houseH    = max(0,(Risky_assets_houseH)/financial_wealth_houseH)
gen portfolio_houseNH   = max(0,(Risky_assets_houseNH)/financial_wealth_houseNH)  
gen portfolio_houseHNH  = max(0,(Risky_assets_houseH)/financial_wealth_houseNH) 
gen portfolio_debt      = max(0,(Risky_assets)/financial_wealth_debt)  
gen portfolio_debt2     = max(0,(Risky_assets)/financial_wealth_debt2)  
* liquid/illiquid
gen illiquid_share         = max(0,Illiquid_assets/financial_wealth)
gen illiquid_share_houseH  = max(0,Illiquid_assets_houseH/financial_wealth_houseH)
gen illiquid_share_houseNH = max(0,Illiquid_assets_houseNH/financial_wealth_houseNH)
gen illiquid_share_debt    = max(0,Illiquid_assets/financial_wealth_debt)
gen illiquid_share_debt2   = max(0,Illiquid_assets/financial_wealth_debt2)

sum portfolio_debt, detail
gen tail=r(p99)
sum portfolio_debt2, detail
gen tail2=r(p99)
* --------------------------------------

gen Retirement_Accounts_S= IRA_safe  + pensions_safe + other_pension_plans_safe
gen Retirement_Accounts_R= IRA_risky + pensions_risky + other_pension_plans_risky
gen old_accounts=  Retirement_Accounts_S + Retirement_Accounts_R
gen share_old =   Retirement_Accounts_R /old_accounts 
 
cd "${main}\Data"
save scf98.dta, replace 

*********************************************************************************
*********************************************************************************
*********************************************************************************




**********************************
*********  2001 Edition  ********* 
**********************************

cd "${main}\Data\Raw"
use p01i6, clear

g year = 2001
gen wgt = X42001

***************************************** 
*********    DEMOGRAPHICS  ************** 
***************************************** 
  
gen male=1 if X8021==1
replace male=0 if (male==.)

gen person1_age = X8022
gen person2_age = X104
gen person3_age = X110
gen person4_age = X116
gen person5_age = X122
gen person6_age = X128
gen age= person1_age

gen     educ= (X5904==1) // college degree

gen income= max(X5729,0) // (yearly)
gen children=X5910
replace children=0 if children==-1
gen marital=1 if  (X7372==1) 
replace marital=0 if marital==.
* --------------------------------------

***************************************** 
**************  Assets          ********* 
***************************************** 
  
**** Checking Accounts
gen checking_accounts1 = max(X3506,0) if (X3507==5)
replace  checking_accounts1=0  if ( checking_accounts1 ==.  )
gen checking_accounts2 = max(X3510,0) if (X3511==5)
replace  checking_accounts2=0  if ( checking_accounts2 ==.  )
gen checking_accounts3 = max(X3514,0) if (X3515==5)
replace  checking_accounts3=0  if ( checking_accounts3 ==.  )
gen checking_accounts4 = max(X3518,0) if (X3519==5)
replace  checking_accounts4=0  if ( checking_accounts4 ==.  )
gen checking_accounts5 = max(X3522,0) if (X3523==5)
replace  checking_accounts5=0  if ( checking_accounts5 ==.  )
gen checking_accounts6 = max(X3526,0) if (X3527==5)
replace  checking_accounts6=0  if ( checking_accounts6 ==.  )
 gen checking_accounts = checking_accounts1 +  checking_accounts2  + checking_accounts3 +  checking_accounts4   ///
                        + checking_accounts5   +  checking_accounts6  
drop checking_accounts1  checking_accounts2  checking_accounts3  checking_accounts4 checking_accounts5  checking_accounts6  
gen part_checking_accounts = ( checking_accounts>0)
* --------------------------------------

*** Savings Accounts ***
gen savings_accounts = max(X3804,0) + max(X3807,0) + max(X3810,0)  + max(X3813,0)  + max(X3816,0) 
gen part_savings_accounts = (savings_accounts>0)
* --------------------------------------

*** Money Market Accounts ***
gen MMA1 = max(X3506,0) if (X3507==1)
replace  MMA1=0  if ( MMA1 ==.  )
gen  MMA2 = max(X3510,0) if (X3511==1)
replace  MMA2=0  if ( MMA2 ==.  )
gen MMA3 = max(X3514,0) if (X3515==1)
replace  MMA3=0  if ( MMA3 ==.  )
gen MMA4 = max(X3518,0) if (X3519==1)
replace  MMA4=0  if ( MMA4 ==.  )
gen MMA5 = max(X3522,0) if (X3523==1)
replace  MMA5=0  if ( MMA5 ==.  )
gen MMA6 = max(X3526,0) if (X3527==1)
replace  MMA6=0  if ( MMA6 ==.  )
gen  MMA7= max(0,X3706) + max(0,X3711) + max(0,X3716)
gen money_market_accounts= MMA1+ MMA2 + MMA3 + MMA4 + MMA5 + MMA6 + MMA7
drop MMA1-MMA7
gen part_money_market_accounts= ( money_market_accounts>0)
* --------------------------------------

*** Cerificates of Deposits ***
gen CDS =max(0,X3721) if (X7620<=3)
replace CDS=0 if (CDS==.)
gen part_CDS = ( CDS>0 )
* --------------------------------------

*** Mutual funds ***
* Stock funds  (only the first one should be risky?)
gen  mutual_stock_funds= max(0,X3822) if ( X3821==1)
replace mutual_stock_funds =0 if  (  mutual_stock_funds ==. ) 
* tax free bonds
gen  mutual_taxfree_bfunds= max(0,X3824) if ( X3823==1)
replace mutual_taxfree_bfunds =0 if  (  mutual_taxfree_bfunds ==. ) 
* governemnt backed bonds
gen  mutual_gov_bfunds= max(0,X3826) if ( X3825==1)
replace mutual_gov_bfunds =0 if  (  mutual_gov_bfunds ==. ) 
* other bond funds
gen  mutual_other_bfunds= max(0,X3828) if ( X3827==1)
replace mutual_other_bfunds =0 if  (  mutual_other_bfunds ==. ) 
*combination funds
gen  mutual_comb_funds= max(0,X3830) if ( X3829==1)
replace mutual_comb_funds =0 if  (  mutual_comb_funds ==. ) 

gen mutual_safe = mutual_taxfree_bfunds + mutual_gov_bfunds +  mutual_other_bfunds + (mutual_comb_funds)/2
gen mutual_risky = mutual_stock_funds + (mutual_comb_funds)/2

gen part_mutual_safe = (mutual_safe > 0)
gen part_mutual_risky = (mutual_risky > 0)
* --------------------------------------

*** Bonds ***
gen savings_bonds_safe = X3902 +  X3910 + X3908   // face value               
gen savings_bonds_risky = X3906 + X7634 + X7633   

gen part_bonds_safe = (savings_bonds_safe)
gen part_bonds_risky = (savings_bonds_risky)
* --------------------------------------

*** Life Insurance & Other assets ***
gen life_insurance = max(0,X4006 - X4010)
gen misc_assets = max(0,X4018 + X4022 + X4026 +  X4030   - X4032) 

gen part_life_insurance = (life_insurance>0)
gen part_misc_assets    = (misc_assets>0)
* --------------------------------------

*** Publicly traded Stocks ***
gen  stocks1= max(0,X3915) if ( X3913==1)
replace stocks1 =0 if  (  stocks1 ==. )
gen stocks2= max(0,X7641) if (X7192==3)
replace stocks2=0 if (stocks2==.)

gen stocks =stocks1+stocks2 
gen part_stocks = (stocks >0)
drop stocks1-stocks2
* --------------------------------------

*** Individual Retirement Accounts ***
/*
           1.  *CDs/Bank accounts; "money market"
           2.  *Stock; "mutual funds"
           3.  *Bonds/Similar assets; T-Bills; treasury notes
           4.  Combinations of 1, 2, & 3 ; "mixed"/"diversified" -- NFS
           5.  Combination of 2 & 3 above
           6.  Combination of 1 & 2 above
           11.  Universal life policy or other similar insurance product
           12.  Annuity
           13.  Commodities
           14.  Real estate/mortgages
           15.  Limited partnership/other similar investment
           16.  Brokerage account/cash management account (CMA)
           -7.  *Split/Other
*/

gen IRA_safe = (X3610 + X3620 + X3630)  if inlist(X3631, 1,3,11)     
replace IRA_safe=0 if IRA_safe==. 
replace IRA_safe=  2*( X3610 + X3620 + X3630 )/3 if X3631==4
* if ==4 then it is a combination of bank accounts, stocks and bonds so you have to give 2/3 to safe
replace IRA_safe=  ( X3610 + X3620 + X3630 ) /2 if inlist(X3631, 5,6,-7)                
* give half to safe if both and half to risky

gen IRA_risky = (X3610 + X3620 + X3630)  if inlist(X3631, 2,14,15,16)  
* stocks  or investments or brokerage accounts
replace IRA_risky=0 if IRA_risky==. 
replace IRA_risky= ( X3610 + X3620 + X3630 )/3 if X3631==4
replace IRA_risky=   ( X3610 + X3620 + X3630 )/2 if inlist(X3631, 5,6,-7) 

gen part_IRA_safe  = (IRA_safe>0)
gen part_IRA_risky = (IRA_risky>0)
* --------------------------------------

*** Brokerage Accounts ***
gen brokerage = max(0,X3930-X3932)
gen part_brokerage = (brokerage>0)
* --------------------------------------

*** Annuities ***
/*
           1.  *STOCKS; MUTUAL FUND
           2.  *BONDS/INTEREST; CDS/MONEY MARKET
           3.  REAL ESTATE
           5.  *SPLIT BETWEEN STOCKS/INTEREST; COMBINATION OF 1 & 2; MUTUAL FUNDS AND CD'S
           6.  MIXED OR DIVERSIFIED
           7.  *LIFE INSURANCE/FIXED CONTRACT; ANNUITIES
           8.  TANGIBLE ASSETS OTHER THAN REAL ESTATE
           9.  INTANGIBLE ASSETS, N.E.C.
           -7.  *OTHER
*/

gen annuity_safe = (X6820)  if inlist(X6826, 2,7)      
replace annuity_safe=0 if annuity_safe==. 
replace annuity_safe=  (X6820)/2 if inlist(X6826, 5,6,-7)        
* give half to safe if both and half to risky

gen annuity_risky = (X6820)   if inlist(X6826, 1,3)   
* stocks or real estate
replace annuity_risky=0 if annuity_risky==. 
replace annuity_risky=  (X6820) /2 if inlist(X6826, 5,6,-7) 

gen part_annuity_safe = (annuity_safe>0)
gen part_annuity_risky = (annuity_risky>0)
* --------------------------------------

*** Trusts ***
/*
           1.  *STOCKS; MUTUAL FUND
           2.  *BONDS/INTEREST/CDS/MONEY MARKET
           3.  REAL ESTATE
           5.  *COMBINATION OF 1 & 2; MUTUAL FUNDS AND CD'S
           6.  MIXED OR DIVERSIFIED
           7.  *LIFE INSURANCE/FIXED CONTRACT; ANNUITIES
           8.  TANGIBLE ASSETS OTHER THAN REAL ESTATE
           9.  INTANGIBLE ASSETS, N.E.C.
          -7.  *OTHER
*/
gen trust_safe = (X6835)  if inlist(X6841, 2,7)    
* bonds or life insurance
replace trust_safe=0 if trust_safe==. 
replace trust_safe= (X6835)/2 if inlist(X6841, 5,6,-7)       
* give half to safe if both and half to risky
gen trust_risky = (X6835)   if inlist(X6841, 1,3)   
* stocks or real estate
replace trust_risky=0 if trust_risky==. 
replace trust_risky= (X6835) /2 if inlist(X6841, 5,6,-7)

gen part_trust_safe  = (trust_safe>0)
gen part_trust_risky = (trust_risky>0)
* --------------------------------------

*** Housing ***       
gen farm_value  = X513 + X526
gen mobile_site = X604 + X614 + X623 
gen value_home  = X716  

gen house_mortgage= (X805 + X905 + X1005) 
gen other_loans= X1044

gen lines_of_Credit =  X1108 + X1119 + X1130 + X1136 
replace lines_of_Credit = 0 if  lines_of_Credit<0 

gen land_contract_lend = (X1409 + X1509 + X1609  + X1619)   
gen land_contract_borrow = X1417 + X1517  + X1617 + X1621 
gen real_estate =  X1706*(X1705/10000)  + X1806*(X1805/10000)  + X1906*( X1905/10000 ) + X2002 + X2012       
gen borrowing_housing =  X1715 *(X1705/10000) + X1815*(X1805/10000)  + X1915*(X1905/10000) + X2006 + X2016


gen house_wealth = ( farm_value + mobile_site + value_home )  +  land_contract_lend  +  real_estate
gen house_wealth2= ( farm_value + mobile_site + value_home )  +  real_estate

gen mortgages=  house_mortgage + other_loans +  land_contract_borrow + borrowing_housing
gen mortgages2= house_mortgage + other_loans + borrowing_housing

gen net_house_wealth = house_wealth-mortgages
gen net_house_wealth2 = house_wealth2-mortgages2

gen owner=0 if net_house_wealth==0
replace owner=1 if owner==.
* --------------------------------------

*** Pensions ***
*Did you tel me about this (pension) loan earlier?
* If yes I don't count it 
replace X4229=X4229 if X4230==5     
replace X4229=0 if X4229==. 
replace X4329=X4329 if X4330==5 
replace X4329=0 if X4329==. 
replace X4429=X4429 if X4430==5 
replace X4429=0 if X4429==. 
replace X4829=X4829 if X4830==5 
replace X4829=0 if X4829==. 
replace X4929=X4929 if X4930==5 
replace X4929=0 if X4929==. 
replace X5029=X5029 if X5030==5 
replace X5029=0 if X5029==. 
* Now create the pension accounts net of loans  
gen net_pension1=  max(0,X4226 - X4229)
gen net_pension2=  max(0,X4326 - X4329)
gen net_pension3=  max(0,X4426 - X4429)
gen net_pension4=  max(0,X4826 - X4829)
gen net_pension5=  max(0,X4926 - X4929)
gen net_pension6=  max(0,X5026 - X5029)
* Then break them up in safe/risky 
gen weight_1=0.5

/*
			1.  *Thrift or Savings
			2.  *401K/403B/SRA
            3.  *Profit Sharing
            4.  *Stock purchase/ESOP
            7.  Deferred compensation
            11.  IRA-SEP; IRA Simple (not to be confused with a regular IRA)
            12.  Defined-contribution plan; TIAA-CREF (Teachers
                 Insurance and Annuity Association/College Retirement
                 Equity Fund)
            13.  Money purchase plan
            14.  Tax-deferred annuity (TDA); tax-sheltered annuity (TSA)
            17.  Other type of annuity (include ERISA plans here
                 unless otherwise specified)
            18.  Other salary reduction plan; deferred compensation plan
            24.  Other state/local government plan
            25.  Other federal government plan
            26.  Other type of account
            30.  Cash balance plan
            -7.  *Other; combination
*/

*safe pensions  (I consider 401K both risky and non-risky) 
gen k1=weight_1*net_pension1 if X4216==1  | X4216==2  | X4216==12 |  X4216==13   |  X4216==14  |  X4216==17  | X4216==-7                     
replace k1=0 if k1==. 
gen k2=weight_1*net_pension2 if X4316==1  | X4316==2  | X4316==12  |  X4316==13  |  X4316==14  |  X4316==17 | X4316 ==-7
replace k2=0 if k2==. 
gen k3=weight_1*net_pension3 if X4416==1  | X4416==2  | X4416==12  |  X4416==13  |  X4416==14  |  X4416==17  | X4416 ==-7              
replace k3=0 if k3==. 
gen k4=weight_1*net_pension4 if X4816==1  | X4816==2  | X4816==12  |  X4816==13  |  X4816==14  |  X4816==17 | X4816 ==-7               
replace k4=0 if k4==. 
gen k5=weight_1*net_pension5 if X4916==1  | X4916==2   | X4916==12 |  X4916==13  |  X4916==14  |  X4916==17 | X4916 ==-7                
replace k5=0 if k5==. 
gen k6=weight_1*net_pension6 if X5016==1  | X5016==2   | X5016==12 |  X5016==13  |  X5016==14  |  X5016==17 | X5016 ==-7              
replace k6=0 if k6==. 
gen pensions_safe =  k1 + k2 + k3 + k4 + k5 + k6    

*Risky pensions 
gen k7=net_pension1 if  X4216==3|X4216==4                
replace k7=(1-weight_1)*net_pension1 if inlist(X4216, 1,2,12,13,14,17,-7)
replace k7=0 if k7==. 
 
gen k8=net_pension2 if  X4316==3|X4316==4                   
replace k8=(1-weight_1)*net_pension2 if inlist(X4316, 1,2,12,13,14,17,-7)                 
replace k8=0 if k8==. 
  
gen k9=net_pension3 if   X4416==3|X4416==4                  
replace k9=(1-weight_1)*net_pension3 if inlist(X4416, 1,2,12,13,14,17,-7)                
replace k9=0 if k9==. 


gen k10=net_pension4 if X4816==3|X4816==4                    
replace k10=(1-weight_1)*net_pension4 if inlist(X4816, 1,2,12,13,14,17,-7)                  
replace k10=0 if k10==. 
  
gen k11=net_pension5 if X4916==3|X4916==4 
replace k11=(1-weight_1)*net_pension5 if inlist(X4916, 1,2,12,13,14,17,-7)                                 
replace k11=0 if k11==. 
  
gen k12=net_pension6 if X5016==3|X5016==4                  
replace k12=(1-weight_1)*net_pension6 if inlist(X5016, 1,2,12,13,14,17,-7)                 
replace k12=0 if k12==. 

gen pensions_risky =  k7 + k8 + k9 + k10 + k11 + k12 

*** Pensions Future
gen pensions_safe1 = X5604  if X6491 ==2    
replace pensions_safe1=0 if  pensions_safe1==. 
replace pensions_safe1=  ( X5604 ) /2 if X6491==3 |  X6491==-7
gen pensions_safe2 = X5612  if X6492 ==2    
replace pensions_safe2=0 if  pensions_safe2==. 
replace pensions_safe2=  ( X5612 ) /2 if X6492==3 |  X6492==-7
gen pensions_safe3 = X5620  if X6493 ==2    
replace pensions_safe3=0 if  pensions_safe3==. 
replace pensions_safe3=  ( X5620 ) /2 if X6493==3 |  X6493==-7
gen pensions_safe4 = X5628  if X6494 ==2    
replace pensions_safe4=0 if  pensions_safe4==. 
replace pensions_safe4=  ( X5628 ) /2 if X6494==3 |  X6494==-7
gen pensions_safe5 = X5636  if X6495 ==2    
replace pensions_safe5=0 if  pensions_safe5==. 
replace pensions_safe5=  ( X5636 ) /2 if X6495==3 |  X6495==-7
gen pensions_safe6 = X5644  if X6496 ==2    
replace pensions_safe6=0 if  pensions_safe6==. 
replace pensions_safe6=  ( X5644 ) /2 if X6496==3 |  X6496==-7

gen pensions_safeF= pensions_safe1 + pensions_safe2 + pensions_safe3 + pensions_safe4 + pensions_safe5 + pensions_safe6 

gen pensions_risky1 = X5604  if X6491 ==1 |  X6491 ==4 | X6491 ==5   
replace pensions_risky1=0 if  pensions_risky1==. 
replace pensions_risky1=  ( X5604 ) /2 if X6491==3 |  X6491==-7
gen pensions_risky2 = X5612  if X6492 ==1 |  X6492 ==4 | X6492 ==5   
replace pensions_risky2=0 if  pensions_risky2==. 
replace pensions_risky2=  ( X5612 ) /2 if X6492==3 |  X6492==-7
gen pensions_risky3 = X5620  if X6493 ==1 |  X6493 ==4 | X6493 ==5   
replace pensions_risky3=0 if  pensions_risky3==. 
replace pensions_risky3=  ( X5620 ) /2 if X6493==3 |  X6493==-7
gen pensions_risky4 = X5628  if X6494 ==1 |  X6494 ==4 | X6494 ==5   
replace pensions_risky4=0 if  pensions_risky4==. 
replace pensions_risky4=  ( X5628 ) /2 if X6494==3 |  X6494==-7
gen pensions_risky5 = X5636  if X6495 ==1 |  X6495 ==4 | X6495 ==5   
replace pensions_risky5=0 if  pensions_risky5==. 
replace pensions_risky5=  ( X5636 ) /2 if X6495==3 |  X6495==-7
gen pensions_risky6 = X5644  if X6496 ==1 |  X6496 ==4 | X6496 ==5   
replace pensions_risky6=0 if  pensions_risky6==. 
replace pensions_risky6=  ( X5644 ) /2 if X6496==3 |  X6496==-7

gen pensions_riskyF= pensions_risky1 + pensions_risky2 + pensions_risky3 + pensions_risky4 + pensions_risky5 + pensions_risky6 

gen total_pension_safe  = pensions_safe + pensions_safeF
gen total_pension_risky = pensions_risky + pensions_riskyF

gen part_pension_safe  = (total_pension_safe>0)
gen part_pension_risky = (total_pension_risky>0)
* --------------------------------------

***********************************************
*************** Other Variables *************** 
***********************************************

*** Educational Debt ***
gen education_loan= max(0,X7824) + max(0,X7847) + max(0,X7870) +max(0,X7924) + max(0,X7947)  + max(0,X7970)  
gen part_education_loan = (education_loan>0)
* --------------------------------------

*** Other Consumer Loans ***
gen other1= X2723 if X6842==5 // X6842 - X6847: Is this loan one that you told me about when we talked about your business?
replace other1=0  if other1==.
gen other2= X2740 if X6843==5
replace other2=0  if other2==.
gen other3= X2823 if X6844==5
replace other3=0  if other3==.
gen other4= X2840 if X6845==5
replace other4=0  if other4==.
gen other5= X2923 if X6846==5
replace other5=0  if other5==.
gen other6= X2940 if X6847==5
replace other6=0  if other6==.

gen other_consumer_loans= other1 + other2  +  other3 +  other4 +  other5 +  other6
drop other1-other5
gen part_other_consumer_loans = (other_consumer_loans>0)

*** Credit Card Debt ***
gen credit_cards_debt= max(0,X413) + max(0,X421) + max(0,X424) +max(0,X427) + max(0,X430)  + max(X7575,0) 
gen part_credit_cards_debt = (credit_cards_debt>0)

*** Interest Rates in Consumer Loans ***
/*
gen rate_c1 = X2724 if X2724>0 // anual rate
gen rate_c2 = X2741 if X2741>0
gen rate_c3 = X2824 if X2824>0
gen rate_c4 = X2841 if X2841>0
gen rate_c5 = X2924 if X2924>0
gen rate_c6 = X2941 if X2941>0

egen rate_other_consumer_loans = rowmean(rate_c1 rate_c2 rate_c3 rate_c4 rate_c5 rate_c6) // mean anual rate, conditional on having debt.
*/
* --------------------------------------


**********************************************************
**************        Definitions       ******************
**********************************************************

*** SAFE ASSETS ***
* benchmark definition for safe assets
gen Safe_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe  ///
		       + IRA_safe + total_pension_safe ) 

* safe assets net of debt (not including educational loans)
gen net_Safe_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe ///
			   + IRA_safe + total_pension_safe )  - credit_cards_debt  - other_consumer_loans 

* safe assets net of debt (including educational loans)
gen net_Safe_assets2 = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe  /// 
			   + IRA_safe + total_pension_safe )  - credit_cards_debt  - other_consumer_loans -education_loan
* --------------------------------------

*** LIQUID ASSETS ***
* benchmark definition for liquid assets
gen Liquid_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)

* liquid assets net of debt (not including educational loans)
gen net_Liquid_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)  /// 
			   - credit_cards_debt  - other_consumer_loans			   

* liquid assets net of debt (including educational loans)
gen net_Liquid_assets2 = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)  /// 
			   - credit_cards_debt  - other_consumer_loans -education_loan			   
* --------------------------------------			   

*** RISKY ASSETS ***			  
* benchmark definition for risky assets
gen Risky_assets          = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                + IRA_risky  + total_pension_risky  )

* risky assets including housing net worth				
gen Risky_assets_houseNH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                  + IRA_risky  + total_pension_risky  +  net_house_wealth)

* risky assets including housing gross worth				  
gen Risky_assets_houseH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky  ///
                  + IRA_risky  + total_pension_risky  +  house_wealth)
* --------------------------------------			   

*** ILIQUID ASSETS ***
* benchmark definition for illiquid assets
gen Illiquid_assets          = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  )

* illiquid assets including housing net worth				
gen Illiquid_assets_houseNH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                  + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  +  net_house_wealth)

* illiquid assets including housing gross worth				  
gen Illiquid_assets_houseH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky  ///
                  + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  +  house_wealth)
* --------------------------------------

*** FINANCIAL WEALTH ***		  
gen financial_wealth = Safe_assets + Risky_assets // equivalent to liquid + illiquid
gen financial_wealth_houseH = Safe_assets + Risky_assets_houseH 
gen financial_wealth_houseNH = Safe_assets + Risky_assets_houseNH 
gen financial_wealth_debt = net_Safe_assets + Risky_assets 
gen financial_wealth_debt2 = net_Safe_assets2 + Risky_assets 
* --------------------------------------

*** PORTFOLIO CHOICE ***
* risky/safe
gen risky_share         = max(0,Risky_assets/financial_wealth)  
gen portfolio_houseH    = max(0,(Risky_assets_houseH)/financial_wealth_houseH)
gen portfolio_houseNH   = max(0,(Risky_assets_houseNH)/financial_wealth_houseNH)  
gen portfolio_houseHNH  = max(0,(Risky_assets_houseH)/financial_wealth_houseNH) 
gen portfolio_debt      = max(0,(Risky_assets)/financial_wealth_debt)  
gen portfolio_debt2     = max(0,(Risky_assets)/financial_wealth_debt2)  
* liquid/illiquid
gen illiquid_share         = max(0,Illiquid_assets/financial_wealth)
gen illiquid_share_houseH  = max(0,Illiquid_assets_houseH/financial_wealth_houseH)
gen illiquid_share_houseNH = max(0,Illiquid_assets_houseNH/financial_wealth_houseNH)
gen illiquid_share_debt    = max(0,Illiquid_assets/financial_wealth_debt)
gen illiquid_share_debt2   = max(0,Illiquid_assets/financial_wealth_debt2)

sum portfolio_debt, detail
gen tail=r(p99)
sum portfolio_debt2, detail
gen tail2=r(p99)
* --------------------------------------

gen Retirement_Accounts_S= IRA_safe  + total_pension_safe
gen Retirement_Accounts_R= IRA_risky + total_pension_risky
gen old_accounts=  Retirement_Accounts_S + Retirement_Accounts_R
gen share_old =   Retirement_Accounts_R /old_accounts 
 
cd "${main}\Data"
save scf01.dta, replace 


*********************************************************************************
*********************************************************************************
*********************************************************************************




**********************************
*********  2004 - 2016   ********* 
**********************************

cd "${main}\Data\Raw"

foreach var in 04 07 10 13 16 {
	use p`var'i6.dta, clear
	
gen year=20`var' 
gen wgt = X42001

	
***************************************** 
*********    DEMOGRAPHICS  ************** 
***************************************** 
  
gen male=1 if X8021==1
replace male=0 if (male==.)

gen person1_age = X8022
gen person2_age = X104
gen person3_age = X110
gen person4_age = X116
gen person5_age = X122
gen person6_age = X128
gen age= person1_age

if `var' == 04 | `var' == 07 | `var' == 10 | `var' == 13 {
	gen educ= (X5904==1) // college degree
}
else {
	gen educ = (X5931 > 9) 
}
gen income= max(X5729,0)
gen children=X5910
replace children=0 if children==-1
gen marital=1 if  (X7372==1) 
replace marital=0 if marital==.
	
	
******************************************
**************    ASSETS    ************** 
******************************************

*** Checking Accounts ***
gen checking_accounts1 = max(X3506,0) if (X3507==5)
replace  checking_accounts1=0  if ( checking_accounts1 ==.  )
gen checking_accounts2 = max(X3510,0) if (X3511==5)
replace  checking_accounts2=0  if ( checking_accounts2 ==.  )
gen checking_accounts3 = max(X3514,0) if (X3515==5)
replace  checking_accounts3=0  if ( checking_accounts3 ==.  )
gen checking_accounts4 = max(X3518,0) if (X3519==5)
replace  checking_accounts4=0  if ( checking_accounts4 ==.  )
gen checking_accounts5 = max(X3522,0) if (X3523==5)
replace  checking_accounts5=0  if ( checking_accounts5 ==.  )
gen checking_accounts6 = max(X3526,0) if (X3527==5)
replace  checking_accounts6=0  if ( checking_accounts6 ==.  )
 
 
gen checking_accounts = checking_accounts1 +  checking_accounts2  + checking_accounts3 +  checking_accounts4   ///
                        + checking_accounts5   +  checking_accounts6  
drop checking_accounts1  checking_accounts2  checking_accounts3  checking_accounts4 checking_accounts5  checking_accounts6  
gen part_checking_accounts = ( checking_accounts>0)
* --------------------------------------	
	
*** Savings Accounts ***
gen savings_accounts = max(X3730,0) + max(X3736,0) + max(X3742,0)  + max(X3748,0)  + max(X3754,0) + max(X3760,0)
gen part_savings_accounts = ( savings_accounts>0)

*** Money Market Accounts ***
gen MMA1 = max(X3506,0) if (X3507==1)
replace  MMA1=0  if ( MMA1 ==.  )
gen  MMA2 = max(X3510,0) if (X3511==1)
replace  MMA2=0  if ( MMA2 ==.  )
gen MMA3 = max(X3514,0) if (X3515==1)
replace  MMA3=0  if ( MMA3 ==.  )
gen MMA4 = max(X3518,0) if (X3519==1)
replace  MMA4=0  if ( MMA4 ==.  )
gen MMA5 = max(X3522,0) if (X3523==1)
replace  MMA5=0  if ( MMA5 ==.  )
gen MMA6 = max(X3526,0) if (X3527==1)
replace  MMA6=0  if ( MMA6 ==.  )
gen money_market_accounts= MMA1+ MMA2 + MMA3 + MMA4 + MMA5 + MMA6 
drop MMA1-MMA6
gen part_money_market_accounts= ( money_market_accounts>0)

*** Cerificates of Deposits ***
gen CDS =max(0,X3721) if (X7620<=3)
replace CDS=0 if (CDS==.)
gen part_CDS = ( CDS>0 )	
	
*** Mutual funds ***
* Stock funds  
gen  mutual_stock_funds= max(0,X3822) if ( X3821==1)
replace mutual_stock_funds =0 if  (  mutual_stock_funds ==. ) 
* tax free bonds
gen  mutual_taxfree_bfunds= max(0,X3824) if ( X3823==1)
replace mutual_taxfree_bfunds =0 if  (  mutual_taxfree_bfunds ==. ) 
* governemnt backed bonds
gen  mutual_gov_bfunds= max(0,X3826) if ( X3825==1)
replace mutual_gov_bfunds =0 if  (  mutual_gov_bfunds ==. ) 
* other bond funds
gen  mutual_other_bfunds= max(0,X3828) if ( X3827==1)
replace mutual_other_bfunds =0 if  (  mutual_other_bfunds ==. ) 
*combination funds
gen  mutual_comb_funds= max(0,X3830) if ( X3829==1)
replace mutual_comb_funds =0 if  (  mutual_comb_funds ==. ) 

gen mutual_safe = mutual_taxfree_bfunds + mutual_gov_bfunds +  mutual_other_bfunds + (mutual_comb_funds)/2
gen mutual_risky = mutual_stock_funds + (mutual_comb_funds)/2 

gen part_mutual_safe = (mutual_safe > 0)
gen part_mutual_risky = (mutual_risky > 0)
* --------------------------------------	

*** Bonds ***
gen savings_bonds_safe = X3902 +  X3910 + X3908   // face value               
gen savings_bonds_risky = X3906 + X7634 + X7633   

gen part_bonds_safe = (savings_bonds_safe)
gen part_bonds_risky = (savings_bonds_risky)
* --------------------------------------

*** Life Insurance & Other assets ***
gen life_insurance = max(0,X4006 - X4010)
gen misc_assets = max(0,X4018 + X4022 + X4026 +  X4030   - X4032) 

gen part_life_insurance = (life_insurance>0)
gen part_misc_assets    = (misc_assets>0)
* --------------------------------------

*** Publicly traded Stocks ***
gen  stocks1= max(0,X3915) if ( X3913==1)
replace stocks1 =0 if  (  stocks1 ==. )
gen stocks2= max(0,X7641) if (X7192==3)
replace stocks2=0 if (stocks2==.)

gen stocks =stocks1+stocks2 
gen part_stocks = (stocks >0)
drop stocks1-stocks2
* --------------------------------------

*** Individual Retirement Accounts ***
gen ira1 = X6551 + X6552 + X6553 + X6554
gen ira2 = X6559 + X6560 + X6561 + X6562
gen ira3 = X6567 + X6568 + X6569 + X6570
gen IRA_safe1 = ira1  if X6555 ==2     
replace IRA_safe1=0 if IRA_safe1==. 
replace IRA_safe1=  ( ira1 ) /2 if X6555==3 |  X6555==-7
gen IRA_safe2 = ira2  if X6563 ==2     
replace IRA_safe2=0 if IRA_safe2==. 
replace IRA_safe2=  ( ira2 ) /2 if X6563==3 |  X6563==-7
gen IRA_safe3 = ira3  if X6571 ==2     
replace IRA_safe3=0 if IRA_safe3==. 
replace IRA_safe3=  ( ira3 ) /2 if X6571==3 |  X6571==-7

gen IRA_safe = IRA_safe1 + IRA_safe2 + IRA_safe3

gen IRA_risky1 = ira1  if X6555 ==1 |   X6555 ==4 | X6555 ==5   
replace IRA_risky1=0 if IRA_risky1==. 
replace IRA_risky1=  ( ira1 ) /2 if X6555==3 |  X6555==-7
gen IRA_risky2 = ira2  if X6563 ==1 |   X6563 ==4 | X6563 ==5   
replace IRA_risky2=0 if IRA_risky2==. 
replace IRA_risky2=  ( ira2 ) /2 if X6563==3 |  X6563==-7
gen IRA_risky3 = ira3  if X6571 ==1 |   X6571 ==4 | X6571 ==5   
replace IRA_risky3=0 if IRA_risky3==. 
replace IRA_risky3=  ( ira3 ) /2 if X6571==3 |  X6571==-7

gen IRA_risky = IRA_risky1 + IRA_risky2 + IRA_risky3

gen part_IRA_safe  = (IRA_safe>0)
gen part_IRA_risky = (IRA_risky>0)
* --------------------------------------

*** Brokerage Accounts ***
gen brokerage = max(0,X3930-X3932)
gen part_brokerage = (brokerage>0)
* --------------------------------------

*** Annuities ***
/*
            1.  *ALL IN STOCKS
            2.  *ALL IN INTEREST EARNING ASSETS
            3.  *SPLIT
            4.  Real estate
            5.  Hedge fund
            6.  Annuities
            8.  Mineral rights
           -7.  *OTHER 
*/
gen annuity_safe = (X6577)  if X6581 ==2      
replace annuity_safe=0 if annuity_safe==. 
replace annuity_safe=  (X6577)/2 if X6581==3| X6581==6 |X6581==-7      
* give half to safe if both and half to risky
gen annuity_risky = (X6577)   if X6581 ==1  |X6581==4 | X6581==5
* stocks or real estate
replace annuity_risky=0 if annuity_risky==. 
replace annuity_risky=  (X6577) /2 if X6581==3| X6581== 6 |X6581==-7
* --------------------------------------

*** Trust ***
/*
            1.  *ALL IN STOCKS
            2.  *ALL IN INTEREST EARNING ASSETS
            3.  *SPLIT
            4.  Real estate
            5.  Hedge fund
            6.  Annuities
            8.  Mineral rights
           -7.  *OTHER 
*/
gen trust_safe = (X6587)  if X6591 ==2    
* bonds or life insurance
replace trust_safe=0 if trust_safe==. 
replace trust_safe= (X6587)/2 if X6591==3| X6591==6 | X6591==-7      
* give half to safe if both and half to risky
gen trust_risky = (X6587)   if X6591 ==1 |X6591==4  |X6591==5
* stocks or real estate
replace trust_risky=0 if trust_risky==. 
replace trust_risky= (X6587) /2 if X6591==3| X6591== 6 | X6591==-7
* --------------------------------------

*** Housing ***         
gen farm_value= X513 + X526
gen mobile_site= X604 + X614 + X623 
gen value_home =  X716  
gen house_mortgage= (X805 + X905 + X1005) 
gen other_loans= X1044
gen lines_of_Credit =  X1108 + X1119 + X1130 + X1136 
replace lines_of_Credit = 0 if  lines_of_Credit<0 

if `var' == 04 | `var' == 07 {  
gen land_contract_lend = (X1409 + X1509 + X1609  + X1619)   
gen land_contract_borrow = X1417 + X1517  + X1617 + X1621 
gen real_estate =  X1706*(X1705/10000)  + X1806*(X1805/10000)  + X1906*( X1905/10000 ) + X2002 + X2012       
gen borrowing_housing =  X1715 *(X1705/10000) + X1815*(X1805/10000)  + X1915*(X1905/10000) + X2006 + X2016
}
else if `var' == 10 {
	gen land_contract_lend = (X1409 + X1509   + X1619)   
	gen land_contract_borrow = X1417 + X1517 + X1621 
	gen real_estate =  X1706*(X1705/10000)  + X1806*(X1805/10000) + X2002 + X2012       
	gen borrowing_housing =  X1715 *(X1705/10000) + X1815*(X1805/10000) + X2006 + X2016
}
else if `var' == 13 | `var' == 16  {
	gen land_contract_lend = (X1310 + X1328 + X1339)
	gen land_contract_borrow = X1318 + X1337 + X1342 
	gen real_estate =  X1706*(X1705/10000)  + X1806*(X1805/10000) + X2002 + X2012       
	gen borrowing_housing =  X1715 *(X1705/10000) + X1815*(X1805/10000) + X2006 + X2016
}

gen house_wealth = ( farm_value + mobile_site + value_home )  +  land_contract_lend  +  real_estate
gen house_wealth2= ( farm_value + mobile_site + value_home )  +  real_estate
gen mortgages=  house_mortgage + other_loans +  land_contract_borrow + borrowing_housing
gen mortgages2= house_mortgage + other_loans + borrowing_housing
gen net_house_wealth = house_wealth-mortgages
gen net_house_wealth2 = house_wealth2-mortgages2
gen owner=0 if net_house_wealth==0
replace owner=1 if owner==.
* --------------------------------------

*** Pensions ***
gen pensions_safe1 = X11032  if X11036 ==2    
replace pensions_safe1=0 if  pensions_safe1==. 
replace pensions_safe1=  ( X11032 ) /2 if X11036==3 |  X11036==-7

gen pensions_safe2 = X11132  if X11136 ==2    
replace pensions_safe2=0 if  pensions_safe2==. 
replace pensions_safe2=  ( X11132 ) /2 if X11136==3 |  X11136==-7

if `var' == 04 | `var' == 07 {  
gen pensions_safe3 = X11232  if X11236 ==2    
replace pensions_safe3=0 if  pensions_safe3==. 
replace pensions_safe3=  ( X11232 ) /2 if X11236==3 |  X11236==-7

gen pensions_safe6 = X11532  if X11536 ==2    
replace pensions_safe6=0 if  pensions_safe6==. 
replace pensions_safe6=  ( X11532 ) /2 if X11536==3 |  X11536==-7
}

gen pensions_safe4 = X11332  if X11336 ==2    
replace pensions_safe4=0 if  pensions_safe4==. 
replace pensions_safe4=  ( X11332 ) /2 if X11336==3 |  X11336==-7

gen pensions_safe5 = X11432  if X11436 ==2    
replace pensions_safe5=0 if  pensions_safe5==. 
replace pensions_safe5=  ( X11432 ) /2 if X11436==3 |  X11436==-7


if `var' == 04 | `var' == 07 {
gen pensions_safe= pensions_safe1 + pensions_safe2 + pensions_safe3 + pensions_safe4 + pensions_safe5 + pensions_safe6 

drop pensions_safe1-pensions_safe5
}
else {
	gen pensions_safe= pensions_safe1 + pensions_safe2 + pensions_safe4 + pensions_safe5  

	drop pensions_safe1-pensions_safe5
}

gen pensions_risky1 = X11032  if X11036 ==1  | X11036 ==4 | X11036 ==5     
replace pensions_risky1=0 if  pensions_risky1==. 
replace pensions_risky1=  ( X11032 ) /2 if X11036==3 |  X11036==-7

gen pensions_risky2 = X11132  if X11136 ==1  | X11136 ==4 | X11136 ==5     
replace pensions_risky2=0 if  pensions_risky2==. 
replace pensions_risky2=  ( X11132 ) /2 if X11136==3 |  X11136==-7

if `var' == 04 | `var' == 07 {
gen pensions_risky3 = X11232  if X11236 ==1  | X11236 ==4 | X11236 ==5     
replace pensions_risky3=0 if  pensions_risky3==. 
replace pensions_risky3=  ( X11232 ) /2 if X11236==3 |  X11236==-7

gen pensions_risky6 = X11532  if X11536 ==1  | X11536 ==4 | X11536 ==5     
replace pensions_risky6=0 if  pensions_risky6==. 
replace pensions_risky6=  ( X11532 ) /2 if X11536==3 |  X11536==-7
}

gen pensions_risky4 = X11332  if X11336 ==1  | X11336 ==4 | X11336 ==5     
replace pensions_risky4=0 if  pensions_risky4==. 
replace pensions_risky4=  ( X11332 ) /2 if X11336==3 |  X11336==-7

gen pensions_risky5 = X11432  if X11436 ==1  | X11436 ==4 | X11436 ==5     
replace pensions_risky5=0 if  pensions_risky5==. 
replace pensions_risky5=  ( X11432 ) /2 if X11436==3 |  X11436==-7

if `var' == 04 | `var' == 07 {
gen pensions_risky=  pensions_risky1 +  pensions_risky2+  pensions_risky3 +  pensions_risky4 +  pensions_risky5 +  pensions_risky6

drop pensions_risky1-pensions_risky5
}
else {
	gen pensions_risky=  pensions_risky1 +  pensions_risky2 +  pensions_risky4 +  pensions_risky5 

	drop pensions_risky1-pensions_risky5
}
* --------------------------------------

*** Pensions Future ***
gen pensions_safe1 = X5604  if X6962 ==2    
replace pensions_safe1=0 if  pensions_safe1==. 
replace pensions_safe1=  ( X5604 ) /2 if X6962==3 |  X6962==-7
gen pensions_safe2 = X5612  if X6968 ==2    
replace pensions_safe2=0 if  pensions_safe2==. 
replace pensions_safe2=  ( X5612 ) /2 if X6968==3 |  X6968==-7
gen pensions_safe3 = X5620  if X6974 ==2    
replace pensions_safe3=0 if  pensions_safe3==. 
replace pensions_safe3=  ( X5620 ) /2 if X6974==3 |  X6974==-7
gen pensions_safe4 = X5628  if X6980 ==2    
replace pensions_safe4=0 if  pensions_safe4==. 
replace pensions_safe4=  ( X5628 ) /2 if X6980==3 |  X6980==-7

if `var' == 04 | `var' == 07 {
gen pensions_safe5 = X5636  if X6986 ==2    
replace pensions_safe5=0 if  pensions_safe5==. 
replace pensions_safe5=  ( X5636 ) /2 if X6986==3 |  X6986==-7
gen pensions_safe6 = X5644  if X6992 ==2    
replace pensions_safe6=0 if  pensions_safe6==. 
replace pensions_safe6=  ( X5644 ) /2 if X6992==3 |  X6992==-7

gen pensions_safeF= pensions_safe1 + pensions_safe2 + pensions_safe3 + pensions_safe4 + pensions_safe5 + pensions_safe6 
}
else { 
gen pensions_safeF= pensions_safe1 + pensions_safe2 + pensions_safe3 + pensions_safe4 
}

gen pensions_risky1 = X5604  if X6962 ==1 |  X6962 ==4 | X6962 ==5   
replace pensions_risky1=0 if  pensions_risky1==. 
replace pensions_risky1=  ( X5604 ) /2 if X6962==3 |  X6962==-7
gen pensions_risky2 = X5612  if X6968 ==1 |  X6968 ==4 | X6968 ==5   
replace pensions_risky2=0 if  pensions_risky2==. 
replace pensions_risky2=  ( X5612 ) /2 if X6968==3 |  X6968==-7
gen pensions_risky3 = X5620  if X6974 ==1 |  X6974 ==4 | X6974 ==5   
replace pensions_risky3=0 if  pensions_risky3==. 
replace pensions_risky3=  ( X5620 ) /2 if X6974==3 |  X6974==-7
gen pensions_risky4 = X5628  if X6980 ==1 |  X6980 ==4 | X6980 ==5   
replace pensions_risky4=0 if  pensions_risky4==. 
replace pensions_risky4=  ( X5628 ) /2 if X6980==3 |  X6980==-7

if `var' == 04 | `var' == 07 {
gen pensions_risky5 = X5636  if X6986 ==1 |  X6986 ==4 | X6986 ==5   
replace pensions_risky5=0 if  pensions_risky5==. 
replace pensions_risky5=  ( X5636 ) /2 if X6986==3 |  X6986==-7
gen pensions_risky6 = X5644  if X6992 ==1 |  X6992 ==4 | X6992 ==5   
replace pensions_risky6=0 if  pensions_risky6==. 
replace pensions_risky6=  ( X5644 ) /2 if X6992==3 |  X6992==-7

gen pensions_riskyF= pensions_risky1 + pensions_risky2 + pensions_risky3 + pensions_risky4 + pensions_risky5 + pensions_risky6 
}
else {
	gen pensions_riskyF= pensions_risky1 + pensions_risky2 + pensions_risky3 + pensions_risky4 
}

* Total Pensions 
gen total_pension_safe  = pensions_safe + pensions_safeF
gen total_pension_risky = pensions_risky + pensions_riskyF

gen part_pension_safe  = (total_pension_safe>0)
gen part_pension_risky = (total_pension_risky>0)
* --------------------------------------

***********************************************
*************** Other Variables *************** 
***********************************************

*** Educational Debt ***
gen education_loan= max(0,X7824) + max(0,X7847) + max(0,X7870) +max(0,X7924) + max(0,X7947)  + max(0,X7970)  
gen part_education_loan = (education_loan>0)
* --------------------------------------

*** Other Consumer Loans ***
gen other1= X2723 if X6842==5 // X6842 - X6847: Is this loan one that you told me about when we talked about your business?
replace other1=0  if other1==.
gen other2= X2740 if X6843==5
replace other2=0  if other2==.
gen other3= X2823 if X6844==5
replace other3=0  if other3==.
gen other4= X2840 if X6845==5
replace other4=0  if other4==.
gen other5= X2923 if X6846==5
replace other5=0  if other5==.
gen other6= X2940 if X6847==5
replace other6=0  if other6==.

gen other_consumer_loans= other1 + other2  +  other3 +  other4 +  other5 +  other6
drop other1-other5
gen part_other_consumer_loans = (other_consumer_loans>0)

*** Credit Card Debt ***
if `var' == 04 | `var' == 07 {
gen credit_cards_debt= max(0,X413) + max(0,X421) + max(0,X424) +max(0,X427) + max(0,X430)  + max(X7575,0) 
}
else if `var' == 10 | `var' == 13 {
	gen credit_cards_debt= max(0,X413) + max(0,X421) + max(0,X427) + max(0,X430)  + max(X7575,0) 
}
else if `var' == 16 {
		gen credit_cards_debt= max(0,X413) + max(0,X421) + max(0,X427)  + max(X7575,0) 
}

gen part_credit_cards_debt = (credit_cards_debt>0)

*** Interest Rates in Consumer Loans ***
/*
gen rate_c1 = X2724 if X2724>0 // anual rate
gen rate_c2 = X2741 if X2741>0
gen rate_c3 = X2824 if X2824>0
gen rate_c4 = X2841 if X2841>0
gen rate_c5 = X2924 if X2924>0
gen rate_c6 = X2941 if X2941>0

egen rate_other_consumer_loans = rowmean(rate_c1 rate_c2 rate_c3 rate_c4 rate_c5 rate_c6) // mean anual rate, conditional on having debt.
*/
* --------------------------------------

**********************************************************
**************        Definitions       ******************
**********************************************************

*** SAFE ASSETS ***
* benchmark definition for safe assets
gen Safe_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe  ///
		       + IRA_safe + total_pension_safe ) 

* safe assets net of debt (not including educational loans)
gen net_Safe_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe ///
			   + IRA_safe + total_pension_safe )  - credit_cards_debt  - other_consumer_loans 

* safe assets net of debt (including educational loans)
gen net_Safe_assets2 = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe  /// 
			   + IRA_safe + total_pension_safe )  - credit_cards_debt  - other_consumer_loans -education_loan
* --------------------------------------

*** LIQUID ASSETS ***
* benchmark definition for liquid assets
gen Liquid_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)

* liquid assets net of debt (not including educational loans)
gen net_Liquid_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)  /// 
			   - credit_cards_debt  - other_consumer_loans			   

* liquid assets net of debt (including educational loans)
gen net_Liquid_assets2 = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)  /// 
			   - credit_cards_debt  - other_consumer_loans -education_loan			   
* --------------------------------------			   

*** RISKY ASSETS ***			  
* benchmark definition for risky assets
gen Risky_assets          = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                + IRA_risky  + total_pension_risky  )

* risky assets including housing net worth				
gen Risky_assets_houseNH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                  + IRA_risky  + total_pension_risky  +  net_house_wealth)

* risky assets including housing gross worth				  
gen Risky_assets_houseH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky  ///
                  + IRA_risky  + total_pension_risky  +  house_wealth)
* --------------------------------------			   

*** ILIQUID ASSETS ***
* benchmark definition for illiquid assets
gen Illiquid_assets          = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  )

* illiquid assets including housing net worth				
gen Illiquid_assets_houseNH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                  + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  +  net_house_wealth)

* illiquid assets including housing gross worth				  
gen Illiquid_assets_houseH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky  ///
                  + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  +  house_wealth)
* --------------------------------------

*** FINANCIAL WEALTH ***		  
gen financial_wealth = Safe_assets + Risky_assets // equivalent to liquid + illiquid
gen financial_wealth_houseH = Safe_assets + Risky_assets_houseH 
gen financial_wealth_houseNH = Safe_assets + Risky_assets_houseNH 
gen financial_wealth_debt = net_Safe_assets + Risky_assets 
gen financial_wealth_debt2 = net_Safe_assets2 + Risky_assets 
* --------------------------------------

*** PORTFOLIO CHOICE ***
* risky/safe
gen risky_share         = max(0,Risky_assets/financial_wealth)  
gen portfolio_houseH    = max(0,(Risky_assets_houseH)/financial_wealth_houseH)
gen portfolio_houseNH   = max(0,(Risky_assets_houseNH)/financial_wealth_houseNH)  
gen portfolio_houseHNH  = max(0,(Risky_assets_houseH)/financial_wealth_houseNH) 
gen portfolio_debt      = max(0,(Risky_assets)/financial_wealth_debt)  
gen portfolio_debt2     = max(0,(Risky_assets)/financial_wealth_debt2)  
* liquid/illiquid
gen illiquid_share         = max(0,Illiquid_assets/financial_wealth)
gen illiquid_share_houseH  = max(0,Illiquid_assets_houseH/financial_wealth_houseH)
gen illiquid_share_houseNH = max(0,Illiquid_assets_houseNH/financial_wealth_houseNH)
gen illiquid_share_debt    = max(0,Illiquid_assets/financial_wealth_debt)
gen illiquid_share_debt2   = max(0,Illiquid_assets/financial_wealth_debt2)

sum portfolio_debt, detail
gen tail=r(p99)
sum portfolio_debt2, detail
gen tail2=r(p99)
* --------------------------------------

gen Retirement_Accounts_S= IRA_safe  + total_pension_safe
gen Retirement_Accounts_R= IRA_risky + total_pension_risky
gen old_accounts=  Retirement_Accounts_S + Retirement_Accounts_R
gen share_old =   Retirement_Accounts_R /old_accounts 
 
cd "${main}\Data"
save scf`var'.dta, replace 
cd "${main}\Data\Raw"
}



**********************************
*********  2019 Edition   ******** 
**********************************

* equal to 2016, but variables are all with "x" instead of "X".

cd "${main}\Data\Raw"

use p19i6.dta, clear
	
gen year=2019 
gen wgt = x42001

	
***************************************** 
*********    DEMOGRAPHICS  ************** 
***************************************** 
  
gen male=1 if x8021==1
replace male=0 if (male==.)

gen person1_age = x8022
gen person2_age = x104
gen person3_age = x110
gen person4_age = x116
gen person5_age = x122
gen person6_age = x128
gen age= person1_age

else {
	gen educ = (x5931 > 9) 
}
gen income= max(x5729,0)
gen children=x5910
replace children=0 if children==-1
gen marital=1 if  (x7372==1) 
replace marital=0 if marital==.
	
	
******************************************
**************    ASSETS    ************** 
******************************************

*** Checking Accounts ***
gen checking_accounts1 = max(x3506,0) if (x3507==5)
replace  checking_accounts1=0  if ( checking_accounts1 ==.  )
gen checking_accounts2 = max(x3510,0) if (x3511==5)
replace  checking_accounts2=0  if ( checking_accounts2 ==.  )
gen checking_accounts3 = max(x3514,0) if (x3515==5)
replace  checking_accounts3=0  if ( checking_accounts3 ==.  )
gen checking_accounts4 = max(x3518,0) if (x3519==5)
replace  checking_accounts4=0  if ( checking_accounts4 ==.  )
gen checking_accounts5 = max(x3522,0) if (x3523==5)
replace  checking_accounts5=0  if ( checking_accounts5 ==.  )
gen checking_accounts6 = max(x3526,0) if (x3527==5)
replace  checking_accounts6=0  if ( checking_accounts6 ==.  )
 
 
gen checking_accounts = checking_accounts1 +  checking_accounts2  + checking_accounts3 +  checking_accounts4   ///
                        + checking_accounts5   +  checking_accounts6  
drop checking_accounts1  checking_accounts2  checking_accounts3  checking_accounts4 checking_accounts5  checking_accounts6  
gen part_checking_accounts = ( checking_accounts>0)
* --------------------------------------	
	
*** Savings Accounts ***
gen savings_accounts = max(x3730,0) + max(x3736,0) + max(x3742,0)  + max(x3748,0)  + max(x3754,0) + max(x3760,0)
gen part_savings_accounts = ( savings_accounts>0)

*** Money Market Accounts ***
gen MMA1 = max(x3506,0) if (x3507==1)
replace  MMA1=0  if ( MMA1 ==.  )
gen  MMA2 = max(x3510,0) if (x3511==1)
replace  MMA2=0  if ( MMA2 ==.  )
gen MMA3 = max(x3514,0) if (x3515==1)
replace  MMA3=0  if ( MMA3 ==.  )
gen MMA4 = max(x3518,0) if (x3519==1)
replace  MMA4=0  if ( MMA4 ==.  )
gen MMA5 = max(x3522,0) if (x3523==1)
replace  MMA5=0  if ( MMA5 ==.  )
gen MMA6 = max(x3526,0) if (x3527==1)
replace  MMA6=0  if ( MMA6 ==.  )
gen money_market_accounts= MMA1+ MMA2 + MMA3 + MMA4 + MMA5 + MMA6 
drop MMA1-MMA6
gen part_money_market_accounts= ( money_market_accounts>0)

*** Cerificates of Deposits ***
gen CDS =max(0,x3721) if (x7620<=3)
replace CDS=0 if (CDS==.)
gen part_CDS = ( CDS>0 )	
	
*** Mutual funds ***
* Stock funds  
gen  mutual_stock_funds= max(0,x3822) if ( x3821==1)
replace mutual_stock_funds =0 if  (  mutual_stock_funds ==. ) 
* tax free bonds
gen  mutual_taxfree_bfunds= max(0,x3824) if ( x3823==1)
replace mutual_taxfree_bfunds =0 if  (  mutual_taxfree_bfunds ==. ) 
* governemnt backed bonds
gen  mutual_gov_bfunds= max(0,x3826) if ( x3825==1)
replace mutual_gov_bfunds =0 if  (  mutual_gov_bfunds ==. ) 
* other bond funds
gen  mutual_other_bfunds= max(0,x3828) if ( x3827==1)
replace mutual_other_bfunds =0 if  (  mutual_other_bfunds ==. ) 
*combination funds
gen  mutual_comb_funds= max(0,x3830) if ( x3829==1)
replace mutual_comb_funds =0 if  (  mutual_comb_funds ==. ) 

gen mutual_safe = mutual_taxfree_bfunds + mutual_gov_bfunds +  mutual_other_bfunds + (mutual_comb_funds)/2
gen mutual_risky = mutual_stock_funds + (mutual_comb_funds)/2 

gen part_mutual_safe = (mutual_safe > 0)
gen part_mutual_risky = (mutual_risky > 0)
* --------------------------------------	

*** Bonds ***
gen savings_bonds_safe = x3902 +  x3910 + x3908   // face value               
gen savings_bonds_risky = x3906 + x7634 + x7633   

gen part_bonds_safe = (savings_bonds_safe)
gen part_bonds_risky = (savings_bonds_risky)
* --------------------------------------

*** Life Insurance & Other assets ***
gen life_insurance = max(0,x4006 - x4010)
gen misc_assets = max(0,x4018 + x4022 + x4026 +  x4030   - x4032) 

gen part_life_insurance = (life_insurance>0)
gen part_misc_assets    = (misc_assets>0)
* --------------------------------------

*** Publicly traded Stocks ***
gen  stocks1= max(0,x3915) if ( x3913==1)
replace stocks1 =0 if  (  stocks1 ==. )
gen stocks2= max(0,x7641) if (x7192==3)
replace stocks2=0 if (stocks2==.)

gen stocks =stocks1+stocks2 
gen part_stocks = (stocks >0)
drop stocks1-stocks2
* --------------------------------------

*** Individual Retirement Accounts ***
gen ira1 = x6551 + x6552 + x6553 + x6554
gen ira2 = x6559 + x6560 + x6561 + x6562
gen ira3 = x6567 + x6568 + x6569 + x6570
gen IRA_safe1 = ira1  if x6555 ==2     
replace IRA_safe1=0 if IRA_safe1==. 
replace IRA_safe1=  ( ira1 ) /2 if x6555==3 |  x6555==-7
gen IRA_safe2 = ira2  if x6563 ==2     
replace IRA_safe2=0 if IRA_safe2==. 
replace IRA_safe2=  ( ira2 ) /2 if x6563==3 |  x6563==-7
gen IRA_safe3 = ira3  if x6571 ==2     
replace IRA_safe3=0 if IRA_safe3==. 
replace IRA_safe3=  ( ira3 ) /2 if x6571==3 |  x6571==-7

gen IRA_safe = IRA_safe1 + IRA_safe2 + IRA_safe3

gen IRA_risky1 = ira1  if x6555 ==1 |   x6555 ==4 | x6555 ==5   
replace IRA_risky1=0 if IRA_risky1==. 
replace IRA_risky1=  ( ira1 ) /2 if x6555==3 |  x6555==-7
gen IRA_risky2 = ira2  if x6563 ==1 |   x6563 ==4 | x6563 ==5   
replace IRA_risky2=0 if IRA_risky2==. 
replace IRA_risky2=  ( ira2 ) /2 if x6563==3 |  x6563==-7
gen IRA_risky3 = ira3  if x6571 ==1 |   x6571 ==4 | x6571 ==5   
replace IRA_risky3=0 if IRA_risky3==. 
replace IRA_risky3=  ( ira3 ) /2 if x6571==3 |  x6571==-7

gen IRA_risky = IRA_risky1 + IRA_risky2 + IRA_risky3

gen part_IRA_safe  = (IRA_safe>0)
gen part_IRA_risky = (IRA_risky>0)
* --------------------------------------

*** Brokerage Accounts ***
gen brokerage = max(0,x3930-x3932)
gen part_brokerage = (brokerage>0)
* --------------------------------------

*** Annuities ***
/*
            1.  *ALL IN STOCKS
            2.  *ALL IN INTEREST EARNING ASSETS
            3.  *SPLIT
            4.  Real estate
            5.  Hedge fund
            6.  Annuities
            8.  Mineral rights
           -7.  *OTHER 
*/
gen annuity_safe = (x6577)  if x6581 ==2      
replace annuity_safe=0 if annuity_safe==. 
replace annuity_safe=  (x6577)/2 if x6581==3| x6581==6 |x6581==-7      
* give half to safe if both and half to risky
gen annuity_risky = (x6577)   if x6581 ==1  |x6581==4 | x6581==5
* stocks or real estate
replace annuity_risky=0 if annuity_risky==. 
replace annuity_risky=  (x6577) /2 if x6581==3| x6581== 6 |x6581==-7
* --------------------------------------

*** Trust ***
/*
            1.  *ALL IN STOCKS
            2.  *ALL IN INTEREST EARNING ASSETS
            3.  *SPLIT
            4.  Real estate
            5.  Hedge fund
            6.  Annuities
            8.  Mineral rights
           -7.  *OTHER 
*/
gen trust_safe = (x6587)  if x6591 ==2    
* bonds or life insurance
replace trust_safe=0 if trust_safe==. 
replace trust_safe= (x6587)/2 if x6591==3| x6591==6 | x6591==-7      
* give half to safe if both and half to risky
gen trust_risky = (x6587)   if x6591 ==1 |x6591==4  |x6591==5
* stocks or real estate
replace trust_risky=0 if trust_risky==. 
replace trust_risky= (x6587) /2 if x6591==3| x6591== 6 | x6591==-7
* --------------------------------------

*** Housing ***         
gen farm_value= x513 + x526
gen mobile_site= x604 + x614 + x623 
gen value_home =  x716  
gen house_mortgage= (x805 + x905 + x1005) 
gen other_loans= x1044
gen lines_of_Credit =  x1108 + x1119 + x1130 + x1136 
replace lines_of_Credit = 0 if  lines_of_Credit<0 

gen land_contract_lend = (x1310 + x1328 + x1339)
gen land_contract_borrow = x1318 + x1337 + x1342 
gen real_estate =  x1706*(x1705/10000)  + x1806*(x1805/10000) + x2002 + x2012       
gen borrowing_housing =  x1715 *(x1705/10000) + x1815*(x1805/10000) + x2006 + x2016

gen house_wealth = ( farm_value + mobile_site + value_home )  +  land_contract_lend  +  real_estate
gen house_wealth2= ( farm_value + mobile_site + value_home )  +  real_estate
gen mortgages=  house_mortgage + other_loans +  land_contract_borrow + borrowing_housing
gen mortgages2= house_mortgage + other_loans + borrowing_housing
gen net_house_wealth = house_wealth-mortgages
gen net_house_wealth2 = house_wealth2-mortgages2
gen owner=0 if net_house_wealth==0
replace owner=1 if owner==.
* --------------------------------------

*** Pensions ***
gen pensions_safe1 = x11032  if x11036 ==2    
replace pensions_safe1=0 if  pensions_safe1==. 
replace pensions_safe1=  ( x11032 ) /2 if x11036==3 |  x11036==-7

gen pensions_safe2 = x11132  if x11136 ==2    
replace pensions_safe2=0 if  pensions_safe2==. 
replace pensions_safe2=  ( x11132 ) /2 if x11136==3 |  x11136==-7

gen pensions_safe4 = x11332  if x11336 ==2    
replace pensions_safe4=0 if  pensions_safe4==. 
replace pensions_safe4=  ( x11332 ) /2 if x11336==3 |  x11336==-7

gen pensions_safe5 = x11432  if x11436 ==2    
replace pensions_safe5=0 if  pensions_safe5==. 
replace pensions_safe5=  ( x11432 ) /2 if x11436==3 |  x11436==-7

gen pensions_safe= pensions_safe1 + pensions_safe2 + pensions_safe4 + pensions_safe5  
drop pensions_safe1-pensions_safe5


gen pensions_risky1 = x11032  if x11036 ==1  | x11036 ==4 | x11036 ==5     
replace pensions_risky1=0 if  pensions_risky1==. 
replace pensions_risky1=  ( x11032 ) /2 if x11036==3 |  x11036==-7

gen pensions_risky2 = x11132  if x11136 ==1  | x11136 ==4 | x11136 ==5     
replace pensions_risky2=0 if  pensions_risky2==. 
replace pensions_risky2=  ( x11132 ) /2 if x11136==3 |  x11136==-7

gen pensions_risky4 = x11332  if x11336 ==1  | x11336 ==4 | x11336 ==5     
replace pensions_risky4=0 if  pensions_risky4==. 
replace pensions_risky4=  ( x11332 ) /2 if x11336==3 |  x11336==-7

gen pensions_risky5 = x11432  if x11436 ==1  | x11436 ==4 | x11436 ==5     
replace pensions_risky5=0 if  pensions_risky5==. 
replace pensions_risky5=  ( x11432 ) /2 if x11436==3 |  x11436==-7

gen pensions_risky=  pensions_risky1 +  pensions_risky2 +  pensions_risky4 +  pensions_risky5 

drop pensions_risky1-pensions_risky5
* --------------------------------------

*** Pensions Future ***
gen pensions_safe1 = x5604  if x6962 ==2    
replace pensions_safe1=0 if  pensions_safe1==. 
replace pensions_safe1=  ( x5604 ) /2 if x6962==3 |  x6962==-7
gen pensions_safe2 = x5612  if x6968 ==2    
replace pensions_safe2=0 if  pensions_safe2==. 
replace pensions_safe2=  ( x5612 ) /2 if x6968==3 |  x6968==-7
gen pensions_safe3 = x5620  if x6974 ==2    
replace pensions_safe3=0 if  pensions_safe3==. 
replace pensions_safe3=  ( x5620 ) /2 if x6974==3 |  x6974==-7
gen pensions_safe4 = x5628  if x6980 ==2    
replace pensions_safe4=0 if  pensions_safe4==. 
replace pensions_safe4=  ( x5628 ) /2 if x6980==3 |  x6980==-7
 
gen pensions_safeF= pensions_safe1 + pensions_safe2 + pensions_safe3 + pensions_safe4 


gen pensions_risky1 = x5604  if x6962 ==1 |  x6962 ==4 | x6962 ==5   
replace pensions_risky1=0 if  pensions_risky1==. 
replace pensions_risky1=  ( x5604 ) /2 if x6962==3 |  x6962==-7
gen pensions_risky2 = x5612  if x6968 ==1 |  x6968 ==4 | x6968 ==5   
replace pensions_risky2=0 if  pensions_risky2==. 
replace pensions_risky2=  ( x5612 ) /2 if x6968==3 |  x6968==-7
gen pensions_risky3 = x5620  if x6974 ==1 |  x6974 ==4 | x6974 ==5   
replace pensions_risky3=0 if  pensions_risky3==. 
replace pensions_risky3=  ( x5620 ) /2 if x6974==3 |  x6974==-7
gen pensions_risky4 = x5628  if x6980 ==1 |  x6980 ==4 | x6980 ==5   
replace pensions_risky4=0 if  pensions_risky4==. 
replace pensions_risky4=  ( x5628 ) /2 if x6980==3 |  x6980==-7

gen pensions_riskyF= pensions_risky1 + pensions_risky2 + pensions_risky3 + pensions_risky4 

* Total Pensions 
gen total_pension_safe  = pensions_safe + pensions_safeF
gen total_pension_risky = pensions_risky + pensions_riskyF

gen part_pension_safe  = (total_pension_safe>0)
gen part_pension_risky = (total_pension_risky>0)
* --------------------------------------

***********************************************
*************** Other Variables *************** 
***********************************************

*** Educational Debt ***
gen education_loan= max(0,x7824) + max(0,x7847) + max(0,x7870) +max(0,x7924) + max(0,x7947)  + max(0,x7970)  
gen part_education_loan = (education_loan>0)
* --------------------------------------

*** Other Consumer Loans ***
gen other1= x2723 if x6842==5 // X6842 - X6847: Is this loan one that you told me about when we talked about your business?
replace other1=0  if other1==.
gen other2= x2740 if x6843==5
replace other2=0  if other2==.
gen other3= x2823 if x6844==5
replace other3=0  if other3==.
gen other4= x2840 if x6845==5
replace other4=0  if other4==.
gen other5= x2923 if x6846==5
replace other5=0  if other5==.
gen other6= x2940 if x6847==5
replace other6=0  if other6==.

gen other_consumer_loans= other1 + other2  +  other3 +  other4 +  other5 +  other6
drop other1-other5
gen part_other_consumer_loans = (other_consumer_loans>0)

*** Credit Card Debt ***
gen credit_cards_debt= max(0,x413) + max(0,x421) + max(0,x427)  + max(x7575,0) 
gen part_credit_cards_debt = (credit_cards_debt>0)

*** Interest Rates in Consumer Loans ***
/*
gen rate_c1 = x2724 if x2724>0 // anual rate
gen rate_c2 = x2741 if x2741>0
gen rate_c3 = x2824 if x2824>0
gen rate_c4 = x2841 if x2841>0
gen rate_c5 = x2924 if x2924>0
gen rate_c6 = x2941 if x2941>0

egen rate_other_consumer_loans = rowmean(rate_c1 rate_c2 rate_c3 rate_c4 rate_c5 rate_c6) // mean anual rate, conditional on having debt.
*/
* --------------------------------------

**********************************************************
**************        Definitions       ******************
**********************************************************

*** SAFE ASSETS ***
* benchmark definition for safe assets
gen Safe_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe  ///
		       + IRA_safe + total_pension_safe ) 

* safe assets net of debt (not including educational loans)
gen net_Safe_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe ///
			   + IRA_safe + total_pension_safe )  - credit_cards_debt  - other_consumer_loans 

* safe assets net of debt (including educational loans)
gen net_Safe_assets2 = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe +  life_insurance  +   misc_assets +  annuity_safe +  trust_safe  /// 
			   + IRA_safe + total_pension_safe )  - credit_cards_debt  - other_consumer_loans -education_loan
* --------------------------------------

*** LIQUID ASSETS ***
* benchmark definition for liquid assets
gen Liquid_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS ///
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)

* liquid assets net of debt (not including educational loans)
gen net_Liquid_assets = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)  /// 
			   - credit_cards_debt  - other_consumer_loans			   

* liquid assets net of debt (including educational loans)
gen net_Liquid_assets2 = max(0,checking_accounts + savings_accounts + money_market_accounts  + CDS /// 
               + savings_bonds_safe +  mutual_safe  +   misc_assets +  annuity_safe +  trust_safe)  /// 
			   - credit_cards_debt  - other_consumer_loans -education_loan			   
* --------------------------------------			   

*** RISKY ASSETS ***			  
* benchmark definition for risky assets
gen Risky_assets          = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                + IRA_risky  + total_pension_risky  )

* risky assets including housing net worth				
gen Risky_assets_houseNH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                  + IRA_risky  + total_pension_risky  +  net_house_wealth)

* risky assets including housing gross worth				  
gen Risky_assets_houseH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky  ///
                  + IRA_risky  + total_pension_risky  +  house_wealth)
* --------------------------------------			   

*** ILIQUID ASSETS ***
* benchmark definition for illiquid assets
gen Illiquid_assets          = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  )

* illiquid assets including housing net worth				
gen Illiquid_assets_houseNH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky ///
                  + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  +  net_house_wealth)

* illiquid assets including housing gross worth				  
gen Illiquid_assets_houseH = max(0, brokerage  +   stocks  + mutual_risky + annuity_risky + savings_bonds_risky +  trust_risky  ///
                  + IRA_risky  + total_pension_risky + IRA_safe + total_pension_safe  +  house_wealth)
* --------------------------------------

*** FINANCIAL WEALTH ***		  
gen financial_wealth = Safe_assets + Risky_assets // equivalent to liquid + illiquid
gen financial_wealth_houseH = Safe_assets + Risky_assets_houseH 
gen financial_wealth_houseNH = Safe_assets + Risky_assets_houseNH 
gen financial_wealth_debt = net_Safe_assets + Risky_assets 
gen financial_wealth_debt2 = net_Safe_assets2 + Risky_assets 
* --------------------------------------

*** PORTFOLIO CHOICE ***
* risky/safe
gen risky_share         = max(0,Risky_assets/financial_wealth)  
gen portfolio_houseH    = max(0,(Risky_assets_houseH)/financial_wealth_houseH)
gen portfolio_houseNH   = max(0,(Risky_assets_houseNH)/financial_wealth_houseNH)  
gen portfolio_houseHNH  = max(0,(Risky_assets_houseH)/financial_wealth_houseNH) 
gen portfolio_debt      = max(0,(Risky_assets)/financial_wealth_debt)  
gen portfolio_debt2     = max(0,(Risky_assets)/financial_wealth_debt2)  
* liquid/illiquid
gen illiquid_share         = max(0,Illiquid_assets/financial_wealth)
gen illiquid_share_houseH  = max(0,Illiquid_assets_houseH/financial_wealth_houseH)
gen illiquid_share_houseNH = max(0,Illiquid_assets_houseNH/financial_wealth_houseNH)
gen illiquid_share_debt    = max(0,Illiquid_assets/financial_wealth_debt)
gen illiquid_share_debt2   = max(0,Illiquid_assets/financial_wealth_debt2)

sum portfolio_debt, detail
gen tail=r(p99)
sum portfolio_debt2, detail
gen tail2=r(p99)
* --------------------------------------

gen Retirement_Accounts_S= IRA_safe  + total_pension_safe
gen Retirement_Accounts_R= IRA_risky + total_pension_risky
gen old_accounts=  Retirement_Accounts_S + Retirement_Accounts_R
gen share_old =   Retirement_Accounts_R /old_accounts 
 
cd "${main}\Data"
save scf19.dta, replace 

*********************************************************************************
*********************************************************************************
*********************************************************************************

****************
*** CPI DATA ***
****************

wbopendata, country(usa) year(1998:2019) indicator(FP.CPI.TOTL) long clear
rename fp_cpi_totl cpi
replace cpi = cpi/100 
keep if inlist(year, 1998,2001,2004,2007,2010,2013,2016,2019)
save cpi.dta, replace 


*********************************
*** MERGING ALL SURVEYS + CPI ***
*********************************

use scf98, clear

foreach var in 01 04 07 10 13 16 19 {
	append using scf`var'.dta
}

mer m:1 year using cpi.dta

************************
*** ADJUSTING BY CPI ***
************************

keep year wgt cpi share_old old_accounts Retirement_Accounts_R Retirement_Accounts_S portfolio_debt portfolio_debt2 illiquid_share_debt2 ///
illiquid_share_debt illiquid_share_houseNH illiquid_share_houseH illiquid_share portfolio_houseHNH portfolio_houseNH ///
portfolio_houseH risky_share financial_wealth_debt2 financial_wealth_debt financial_wealth_houseNH financial_wealth_houseH financial_wealth ///
Illiquid_assets_houseH Illiquid_assets_houseNH Illiquid_assets Risky_assets_houseH Risky_assets_houseNH Risky_assets net_Liquid_assets2 ///
net_Liquid_assets Liquid_assets net_Safe_assets2 net_Safe_assets Safe_assets checking_accounts savings_accounts money_market_accounts ///
CDS savings_bonds_safe mutual_safe life_insurance misc_assets annuity_safe trust_safe IRA_safe total_pension_safe credit_cards_debt /// 
other_consumer_loans education_loan brokerage stocks mutual_risky annuity_risky savings_bonds_risky trust_risky ///
IRA_risky total_pension_risky net_house_wealth house_wealth income educ male children marital

local deflator cpi
local variables_to_deflate old_accounts Retirement_Accounts_R Retirement_Accounts_S ///
financial_wealth_debt2 financial_wealth_debt financial_wealth_houseNH financial_wealth_houseH financial_wealth ///
Illiquid_assets_houseH Illiquid_assets_houseNH Illiquid_assets Risky_assets_houseH Risky_assets_houseNH Risky_assets net_Liquid_assets2 ///
net_Liquid_assets Liquid_assets net_Safe_assets2 net_Safe_assets Safe_assets checking_accounts savings_accounts money_market_accounts ///
CDS savings_bonds_safe mutual_safe life_insurance misc_assets annuity_safe trust_safe IRA_safe total_pension_safe credit_cards_debt /// 
other_consumer_loans education_loan brokerage stocks mutual_risky annuity_risky savings_bonds_risky trust_risky ///
IRA_risky total_pension_risky net_house_wealth house_wealth income 
foreach var of varlist `variables_to_deflate '{
qui replace `var' = `var'/`deflator'
label var `var' "``var' deflated by cpi"
}


save SCF_Data.dta, replace
