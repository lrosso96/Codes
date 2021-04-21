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
main_dir = 'C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python'
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
    
driver.close()

#%% 

# extract paragraphs from EY tax guides

files = os.listdir(main_dir + tax_guides) # EY worldwide tax guides

text = []
for file in files:
    if file != 'chromedriver.exe': # only file that is not an art. iv:
        doc = fitz.open(file)
        num_pages = doc.pageCount
        for page in range(0,num_pages):
            text_ = doc.loadPage(page).getText("blocks")
            text_ = [x[4] for x in text_] # 4 = element with string in tuple
            text_aux = [t.split('\n \n') for t in text_]
            flat_list = [item for sublist in text_aux for item in sublist]
            text.extend(flat_list)  

# we can identify country sections by titles "country+\n"
jurisdiction_aux = []
for j in range(0,len(jurisdiction)):
    jurisdiction_aux.append(jurisdiction[j]+'\n')

# generate indicator for paragraph with country name title
ind = []
for t in text:
    if t in jurisdiction_aux:
        ind.append(text.index(t))
    else:
        pass

# filter by sample of countries of interest
countries = ['Argentina','Bolivia','Brazil','Chile','Colombia','Costa Rica',
             'Ecuador','El Salvador','Guatemala','Honduras','Jamaica',
             'Mexico','Nicaragua','Panama','Paraguay','Peru','Dominican Republic',
             'Trinidad and Tobago','Uruguay','Venezuela, Bolivarian Republic of']

country_aux1 = []
country_aux2 = []
for c in range(0,len(countries)):
    country_aux2.append(countries[c]+' \n')
    countries[c] = countries[c]+'\n' 
    country_aux1.append('\n'+countries[c])

final_text = []
for i in range(0,len(ind)-1):
    if text[ind[i]] in countries:
        final_text.extend(text[ind[i]:ind[i+1]])
if text[ind[len(ind)-1]] in countries:
    final_text.extend(text[ind[len(ind)-1]:])
    
# adding filters to delete unnecessary info
not_needed = ['Email:','E-mail:','e-mail','Fax','Mobile:','ey.com','GMT','Tel ','Country code',
              'Executive contacts','Ernst & Young','Mail address','Street Address']
not_needed = country_aux1 + country_aux1 + not_needed
filtered_text = [x for x in final_text if not
              any(y in x for y in not_needed)]


# -------------------------------------------------------------
