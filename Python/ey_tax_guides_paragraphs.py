# ------------------------------------------- #
#       Fiscal reforms in Latin America       #
# ------------------------------------------- #

# code to process worldwide tax guides available from EY
# author: Lucas Rosso
# date: 20-04-2021

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
# I will get the full list of countries from EY webpage

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
add_countries = ['Armenia','British Virgin Islands','China Mainland','Congo, Democratic Republic of',
                 'Congo, Republic of','Côte d’Ivoire','Guernsey, Channel Islands','Guinea','Guyana',
                 'Moldova','Monaco','']

jurisdiction = jurisdiction + add_countries
    
driver.close()

#%% 

# extract paragraphs from EY tax guides

files = os.listdir(main_dir + tax_guides) # EY worldwide tax guides

text          = []
filename_aux  = []
for file in files:
    if file != 'chromedriver.exe': # only file that is not an art. iv:
        doc = fitz.open(file)
        num_pages = doc.pageCount
        for page in range(0,num_pages):
            text_ = doc.loadPage(page).getText("blocks")
            text_ = [x[4] for x in text_] # 4 = element with string in tuple
            text.extend(text_)        
            for te in text_:
                filename_aux.append(file) # allows to identify tax guide and year

# we can identify country sections by titles "country+\n"
jurisdiction_aux = []
for j in range(0,len(jurisdiction)):
    jurisdiction_aux.append(jurisdiction[j]+'\n')

# generate indicator for paragraph with country name title
ind = []
for t in range(0,len(text)):
    if text[t] in jurisdiction_aux:
        ind.append(t)
    else:
        pass

# filter by sample of countries of interest
countries = ['Argentina','Bolivia','Brazil','Chile','Colombia','Costa Rica',
             'Ecuador','El Salvador','Guatemala','Honduras','Jamaica',
             'Mexico','Nicaragua','Panama','Paraguay','Peru','Dominican Republic',
             'Trinidad and Tobago','Uruguay','Venezuela, Bolivarian Republic of']

# several aux variables to clean data
country_aux1 = []
country_aux2 = []
country_aux3 = []
country_aux4 = []
country_aux5 = []
for c in range(0,len(countries)):
    country_aux2.append(countries[c]+' \n')
    country_aux3.append('  '+countries[c].lower()+' \n')
    country_aux4.append(' \n'+countries[c].lower()+'  ')
    country_aux5.append(countries[c].lower()+' (continued)\n')
    countries[c] = countries[c]+'\n' 
    country_aux1.append('\n'+countries[c])

# to check if it is identifying the countries correctly
# for i in range(0,40):
#     print(text[ind[i]])
#     print(text[ind[i]] in countries)

final_text  = []
countryname = []
filename    = []
for i in range(0,len(ind)-1):
    if text[ind[i]] in countries:
        #text
        final_text.extend(text[ind[i]:ind[i+1]])
        
        #country
        countryname.extend([text[ind[i]]]*(ind[i+1]-ind[i]))
        
        #filename (allows to identify tax guide and year)
        filename.extend(filename_aux[ind[i]:ind[i+1]])
if text[ind[len(ind)-1]] in countries:
    final_text.extend(text[ind[len(ind)-1]:])
    countryname.extend([text[ind[i]]]*(ind[i+1]-ind[i]))
    filename.extend(filename_aux[ind[i]:ind[i+1]])
    
# adding filters to delete unnecessary info
not_needed = ['Email:','E-mail:','e-mail','Fax','Mobile:','ey.com','GMT','Tel ','Country code',
              'Executive contacts','Ernst & Young','Mail address','Street Address','EY contact',
              '\uf0a8','\uf0fe']
not_needed = country_aux1 + country_aux2 + not_needed
filtered_text = [x for x in final_text if not
              any(y in x for y in not_needed)]

# drop page indicator
country_auxs = country_aux3 + country_aux4 + country_aux5
filtered_text = [x for x in filtered_text if not
              any(y in x.lower() for y in country_auxs)]



# indicators for filename and countryname
ind_1 = []
ind_2 = []
for te in final_text:
    ind_1_ = (not any(y in te for y in not_needed))
    ind_1.append(ind_1_)
    
    ind_2_ = (not any(y in te.lower() for y in country_auxs))
    ind_2.append(ind_2_)

# edit length of filename to match filtered data
final_filename    = []   
final_countryname = []
for i in range(0,len(ind_1)):
    if (ind_1[i] == 1 and ind_2[i] == 1):
        final_filename.append(filename[i])
        final_countryname.append(countryname[i])
        
        

# get year of the report 
year = []
for file in final_filename:
    if file[0:5] == 'ey-20':
        year.append(file[3:7])
    else:
        year.append('20'+file[-6:-4])

# -------------------------------------------------------------

#%%

# first, delete \n to avoid problems in csv
for t in range(0,len(filtered_text)):
    filtered_text[t] = filtered_text[t].replace('\n','')
    filtered_text[t] = filtered_text[t].replace(';',',')
    final_countryname[t] = final_countryname[t].replace('\n','')
    

# Export data as csv file
chdir(main_dir + data)
df = pd.DataFrame({
    'filename':final_filename,
    'iso3c':final_countryname,
    'year':year,
    'text':filtered_text,
})

# export to csv file
df.to_csv('EY_tax_guides.csv',index=False,sep='|',encoding='utf-8')
