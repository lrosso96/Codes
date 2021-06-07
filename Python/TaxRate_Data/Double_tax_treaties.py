# ---------------------------------------------- #
#            Personal tax guide panel            #
# ---------------------------------------------- #

# code to generate panel of corporate taxes from EY
# author: Lucas Rosso
# date: 04-06-2021

'''
This code extracts information on personal taxes from EY worldwide guides.
After extracting and selecting relevant text, the code will collapse the data
and create panel data on this taxes to be merged with corporate income tax data.
For now, the code only gets the list of double tax treaties.
'''

# %% PRELIMINARIES

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


#%%  CREATING LIST OF COUNTRY NAMES
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
personal_ind = 'worldwide-personal-tax-guide-' # personal tax guide indicator

text          = []
filename_aux  = []
for file in files:
    if file[0:len(personal_ind)] == personal_ind: # only corporate tax guides
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

# %% GET DOUBLE TAX RELIEF AND TAX TREATIES

'''
This block looks for the list of countrues with 'Double tax relief 
and tax treaties' for all tax guides.
'''         

tax_treatie = ['double tax relief and tax treaties\n','tax treaties\n' ]
temp_visas  = ['temporary visas\n']

double_tax = []
country    = []
year       = []
for t in range(0,len(text)):
    if any(y in text[t].lower() for y in tax_treatie):
        for i in range(0,10):
            te = text[t+i].lower()
            # if you find a country list: stop
            if 'has not concluded any double tax treaties' in te:
                break
            elif te[:te.index("\n")+1].replace(' ','') in jurisdiction_aux:
                break
            # if you find temp. visas: you've gone to far
            elif any(y in text[t+i].lower() for y in temp_visas):
                break
        if 'has not concluded any double tax treaties' in te:
            double_tax.append('No double tax treaties')
            country.append(countryname[t].lower().replace('\n', ''))
            year.append(year_aux[t])          
        elif te[:te.index("\n")+1].replace(' ','') not in jurisdiction_aux:
            double_tax.append('.')
            country.append(countryname[t].lower().replace('\n', ''))
            year.append(year_aux[t])
        else:
            double_tax.append(te.replace('\n','.'))  
            country.append(countryname[t].lower().replace('\n', ''))
            year.append(year_aux[t])

# %% EXPORT DATA TO CSV

chdir(main_dir + data)
df = pd.DataFrame({
    'country':country,
    'year':year,
    'double_tax':double_tax,
})

# export to csv file
df.to_csv('double_tax_panel.csv',index=False,sep='|',encoding='utf-16')
                          
