# ------------------------------------------- #
#       Fiscal reforms in Latin America       #
# ------------------------------------------- #

# code to generate panel of corporate taxes from EY
# author: Lucas Rosso
# date: 17-05-2021

'''
[New version]: simpler code and solves some bugs from the original
brief overview: this code intends to exploit structure of corporate tax guidees from EY.
In particular, the fact that all country sections start with a summary of tax rates named
"At a Glance". This allows to follow differnt tax rates across years and countries and build
a novel dataset on tax policy.
'''

#%% 
# required packages
import os
import pandas as pd
from os import chdir
main_dir = 'C:/Users/LR/Dropbox/Fiscal_Project/Python'
tax_guides = '/EY_tax_guides' 
data = '/Data'
chdir(main_dir + tax_guides)

import fitz

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
import time 
import numpy as np

url = 'https://www.ey.com/en_gl/tax-alerts'

# chrome webdriver options
options = Options()
options.page_load_strategy = 'normal'
options.add_argument("--start-maximized")
options.add_argument('--disable-extensions')

#%% 
'''
I will get the full list of countries by scraping the EY
website and then adding "troubling names" as an add list.
'''

driver = webdriver.Chrome(options=options) 
driver.get(url)

time.sleep(2)

# clicking in "I decline optional cookies"
try:
    WebDriverWait(driver, 10)\
        .until(EC.presence_of_element_located((By.XPATH, '//*[@id="cookiePolicy"]/div/div[2]/button[1]')))\
        .click()
except:
    pass
time.sleep(np.random.randint(0,5))

# opening the country filter
WebDriverWait(driver, 10)\
        .until(EC.presence_of_element_located((By.XPATH, '//*[@id="technical-content-search"]/div[1]/div[3]/div/div[2]/div[1]/div[4]/div[2]/label')))\
        .click()
time.sleep(np.random.randint(2,5))

# get the list of counties
jurisdiction = []
for i in range(1,152): 
    aux = driver.find_element_by_xpath('//*[@id="technical-content-search"]/div[1]/div[3]/div/div[2]/div[1]/div[4]/div[3]/div['+str(i)+']/label')
    aux.location_once_scrolled_into_view # locate element to get country name
    driver.execute_script("window.scrollTo(0, 400)")
    jurisdiction.append(aux.text)
    
# adding some countries or alternative country names to get cleaner data
add_countries = ['Armenia','Azerbaijan','British Virgin Islands','China Mainland','Congo, Democratic Republic of',
                 'Congo, Republic of','Côte d’Ivoire','Guernsey, Channel Islands','Guinea','Guyana',
                 'Moldova','Monaco','CHINA, PEOPLE’S REPUBLIC OF','CÔTE D’IVOIRE (IVORY COAST)',
                 'China, People\'s Republic of','China (mainland)','Congo, Democratic Republic of ',
                 'Hong Kong Special Administrative Region (SAR) of China',
                 'Hong Kong Special Administrative Region \n(SAR) of China','Bolivia','Eswatini','Georgia',
                 'Guam','Iran','Ireland, Republic of','Jersey, Channel Islands','Korea','Korea (South)','Laos',
                 'Lebanon','Lesotho','Macau','Macedonia, Former Yugoslav Republic of','Madagascar','Malawi',
                 'Maldives','Mauritania','Montenegro, Republic of','Nepal','Netherlands Antilles','Nicaragua',
                 'Northern Mariana Islands, \nCommonwealth of the','Northern Mariana Islands, \nCommonwealth ofthe',
                 'Northern Mariana Islands,\nCommonwealth ofthe','Northern Mariana Islands,\nCommonwealth of the',
                 'Palestinian Authority','St. Lucia','Saint Martin','são Tomé anD príncipE',
                 'Senegal','Serbia, Republic of','Slovak Republic','Sint Maarten','Sudan','Swaziland','Seychelles',
                 'Syria','Taiwan','Tanzania','US Virgin Islands','Venezuela','Kosovo','Macau Special Administrative Region (SAR) of China',
                 'Macau Special Administrative Region \n(SAR) of China','New Caledonia','Palestine','Serbia and Montenegro, Union of',
                 'Part 2: Republic of Montenegro','Saint-Martin','U.S. Virgin Islands','Zimbabwe']

jurisdiction_plus = jurisdiction + add_countries
    
driver.close()


#%% 
'''
This block extracts all paragraphs from EY corporate tax guides
'''

files = os.listdir(main_dir + tax_guides) # EY worldwide tax guide directory
corp_ind = 'ey-worldwide-corporate-tax-guide-' # corporate taxt guide indicator

text          = []
filename_aux  = []
for file in files:
    if file[0:len(corp_ind)] == corp_ind: # only corporate tax guides
        doc = fitz.open(file)
        num_pages = doc.pageCount
        for page in range(0,num_pages):
            text_ = doc.loadPage(page).getText("blocks")
            text_ = [x[4] for x in text_] # 4 = element with string in tuple
            text.extend(text_)        
            for te in text_:
                filename_aux.append(file) # allows to identify tax guide and year
                
                
#%% 
     
'''
In this block I will identify the country and the year from both the filename
and the text within each file.
'''
            
# we can identify country sections by titles "country+\n"
jurisdiction_aux = []
for j in range(0,len(jurisdiction_plus)):
    jurisdiction_aux.append(jurisdiction_plus[j].lower().replace(' ','')+'\n') # lower() to avoid problems with capital letters

# generate indicator for paragraph with country name title
ind = []
for t in range(0,len(text)):
    if text[t].lower().replace(' ','') in jurisdiction_aux:
        ind.append(t)
    else:
        pass

# Extract country name from chapter titles
countryname = []
for t in range(0,len(text)):
    if text[t].lower().replace(' ','') in jurisdiction_aux:
        countryname.append(text[t])
    else:
        try:
            countryname.append(countryname[-1])
        except:
            countryname.append('.')

# Get the year from the filename            
year_aux = []
for f in filename_aux:
    year_aux.append(f[-8:-4])
    
                
#%%
'''
Now we need to identify the section "At a glance" and extract all specific 
corporate tax rates.  
'''

withholding = ['Withholding Tax (%)', 'Withholding Tax Rate (%)', 'Withholding Tax* (%)', 'Withholding Tax (% )',
               'Withholding Tax (%) (a)','Withholding Tax (%)(a)','Withholding Tax (%) (b)','Withholding Tax (%)(b)','Yield Tax Rate (%)',
               'Withholding Tax (%) (c)','Withholding Tax (%)(c)']

# Some countries have both withholding tax and Capital Gains Withholding Tax Rate (%) which is troubling
for t in range(0,len(text)):
    text[t] = text[t].replace('Capital Gains Withholding Tax Rate (%)','Not_Interested')
    

# get indicator for "At a glance" and "Withholding Tax"
ind_1 = []
ind_2 = []
for t in range(0,len(text)):
    # both on the same element (the last condition is due to problems with Zimbawe 2016)
    if 'a. at a glance' in text[t].lower() and any(y in text[t] for y in withholding): # and 'Capital Gains Withholding Tax Rate (%)' not in text[t]:
        ind_1.append(t)
        ind_2.append(t)
    # only 'at a glance'
    elif 'a. at a glance' in text[t].lower():
        ind_1.append(t)
    # only 'Withholding Tax'
    elif any(y in text[t] for y in withholding): # and 'Capital Gains Withholding Tax Rate (%)' not in text[t]: 
        ind_2.append(t)
    elif countryname[t].replace(' ','') == 'Philippines\n' and 'Withholding Tax \n' in text[t]:
        ind_2.append(t)

# use both indicators to get the desired text, country and filename (the last one to get the year)
at_a_glance  = []
country_name = []
year         = []
for i in range(0,len(ind_1)-1):
    if len(ind_1) != len(ind_2):
        print('Different length')
        break
    else:
        at_a_glance.extend(text[ind_1[i]:ind_2[i]+1])
        country_name.extend(countryname[ind_1[i]:ind_2[i]+1])
        year.extend(year_aux[ind_1[i]:ind_2[i]+1])

#%%

'''
Given the lists at_a_glance, country_name and year we already have the data we need. 
Now I will try to give the data a panel structure.
'''

# first create list for extra words
not_needed = ['A. At a Glance','142 CHINA','Withholding Tax','Net Operating Losses']

special_cases = ['Formulas and Similar Items','Other Fees and Compensation for',
                  'Consulting Services','Rendered in Panama']

text_panel  = []
countryname = []
time        = []
for i in range(0,len(at_a_glance)):
    # delete countrynames and page indicators (e.g. 142 CHINA)
    
    # split the text
    aux = at_a_glance[i].split('\n') # separate blocks using \n
    
    # filter information
    # greater than zero drops empty obs, and not_needed list extra text
    aux = [x for x in aux if len(x) > 0 and not
              any(y in x for y in not_needed)]
    
    # deal with long tax entries (keep only first line)
    aux = [x for x in aux if not any(y == x for y in special_cases)]  
    
    # get the text    
    text_panel.extend(aux)
    
    # get the year
    time.extend([year[i]]*len(aux))
    
    # get the country (and get rid of \n in the country name)
    country_name[i] = country_name[i].replace('\n','').lower()
    countryname.extend([country_name[i]]*len(aux))
   
    
#%%

'''
We now have 2 columns (tax and rate) in one in text_panel. Thus, this block will 
split text into two variables.
'''

# first some cleaning
symbols = ['(%)', '(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '*'] 
for i in range(0,len(text_panel)):
    for s in symbols:
        text_panel[i] = text_panel[i].replace(s,'')

# tax categories we want 
taxes = ['CorporateIncomeTaxRate','CapitalGainsTaxRate','BranchTaxRate','Dividends']

tax     = []
rate    = []
country = []
year    = []
for t in range(0,len(text_panel)):
    if text_panel[t].replace(' ','') in taxes:
        if len(tax) == 0:
            # tax and rate
            tax.append(text_panel[t].replace(' ','')) # remove spaces will make the reshape easier later in Stata
            rate.append(text_panel[t+1].strip())
        
            # country and year
            year.append(time[t])
            country.append(countryname[t].strip())
        elif len(tax)>0 and text_panel[t].replace(' ','') != tax[len(tax)-1]:
            # tax and rate
            tax.append(text_panel[t].replace(' ','')) # remove spaces will make the reshape easier later in Stata
            rate.append(text_panel[t+1].strip())
        
            # country and year
            year.append(time[t])
            country.append(countryname[t].strip())
            
# final edits to some countrynames
for i in range(0,len(country)):
    country[i] = country[i].replace(', people’s republic of','')
    country[i] = country[i].replace('ireland, republic of','ireland')

#%%
'''
Finally, generate the dataframe and export to a csv file
'''

chdir(main_dir + data)
df = pd.DataFrame({
    'country':country,
    'year':year,
    'tax':tax,
    'rate':rate,
})

# export to csv file
df.to_csv('corp_tax_panel_full.csv',index=False,sep='|',encoding='utf-16')
