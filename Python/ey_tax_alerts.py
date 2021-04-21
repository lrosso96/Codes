# ------------------------------------------- #
#       Fiscal reforms in Latin America       #
# ------------------------------------------- #

# code to download and store tax alert news from EY
# author: Lucas Rosso
# date: 07-04-2021

# %%
# loading libraries
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
import time 
import numpy as np
import os
from os import chdir
chdir("C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python")

url = 'https://www.ey.com/en_gl/tax-alerts'

# chrome webdriver options
options = Options()
options.page_load_strategy = 'normal'
options.add_argument("--start-maximized")
options.add_argument('--disable-extensions')

# %%

driver = webdriver.Chrome(options=options) 
driver.get(url)

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

# filter by sample of countries
countries = ['Argentina','Bolivia','Brazil','Chile','Colombia','Costa Rica',
             'Ecuador','El Salvador','Guatemala','Honduras','Jamaica',
             'Mexico','Nicaragua','Panama','Paraguay','Peru','Dominican Republic',
             'Trinidad and Tobago','Uruguay','Venezuela, Bolivarian Republic of']

# scroll up the page to leave space for the button to be clicked
driver.execute_script("window.scrollTo(0, 400)")

# list of countries to scroll down
country_list = driver.find_element_by_xpath('//*[@id="technical-content-search"]/div[1]/div[3]/div/div[2]/div[1]/div[4]/div[3]')

# loop over countries to click the ones in the list
for i in range(1,152): 
    country_name = driver.find_element_by_xpath('//*[@id="technical-content-search"]/div[1]/div[3]/div/div[2]/div[1]/div[4]/div[3]/div['+str(i)+']/label')
    country_name.location_once_scrolled_into_view # locate element to get country name
    driver.execute_script("window.scrollTo(0, 400)")
    country_name = country_name.text
    time.sleep(np.random.randint(0,1))
    if any(y in country_name for y in countries):
        print(country_name)
        time.sleep(np.random.randint(2,5))
        WebDriverWait(driver, 20)\
            .until(EC.presence_of_element_located((By.XPATH, '//*[@id="technical-content-search"]/div[1]/div[3]/div/div[2]/div[1]/div[4]/div[3]/div['+str(i)+']/label')))\
            .click()
        # driver.find_element_by_xpath('//*[@id="technical-content-search"]/div[1]/div[3]/div/div[2]/div[1]/div[4]/div[3]/div['+str(i)+']/label').click()
    else:
        pass


# %%
# with countries selected, now we need to extract tax alerts

driver.execute_script("window.scrollTo(0, 400)")

# creating lists
country     = []
tax_type    = []
date_alert  = []
paragraphs  = []

while True:
    #elements in page
    elements = driver.find_elements_by_class_name('item')
    for i in range(1,len(elements)+1):
        # extract date, country and type of tax
        # date_alert.append(driver.find_element_by_xpath('//*[@id="technical-content-search"]/div[1]/div[4]/div/div/div/div/div['+str(i)+']/div/span[1]').text)
        # tax_type.append(driver.find_element_by_xpath('//*[@id="technical-content-search"]/div[1]/div[4]/div/div/div/div/div['+str(i)+']/div/span[2]').text)
        # country.append(driver.find_element_by_xpath('//*[@id="technical-content-search"]/div[1]/div[4]/div/div/div/div/div['+str(i)+']/div/span[3]').text)
    
        driver.execute_script("window.scrollTo(0, 400)")
        time.sleep(np.random.randint(3,6))
    
        # click to enter the news    
        WebDriverWait(driver, 20)\
            .until(EC.presence_of_element_located((By.XPATH, '//*[@id="technical-content-search"]/div[1]/div[4]/div/div/div/div/div['+str(i)+']/h4/a')))\
            .click()   
    
        # getting paragraphs, date, country and type of tax
        p_elements = driver.find_elements_by_tag_name('p')
        for p in p_elements:
            time.sleep(np.random.randint(0,2))
            # all paragraphs
            paragraphs.append(p.text)
        
            #date, tax type and country info
            date_alert.append(driver.find_element_by_xpath('//*[@id="technical-content-left"]/div/div[1]/div/div/div[2]/div[1]/span').text)
            tax_type.append(driver.find_element_by_xpath('//*[@id="technical-content-left"]/div/div[1]/div/div/div[2]/div[3]/a').text)
            country.append(driver.find_element_by_xpath('//*[@id="technical-content-left"]/div/div[1]/div/div/div[2]/div[4]/a').text)

        time.sleep(np.random.randint(3,5))
        driver.back()
    
    time.sleep(np.random.randint(3,5))
    
    # scroll the page up to the right spot
    # driver.find_element_by_xpath('//*[@id="technical-content-search"]/div[2]/button[2]').location_once_scrolled_into_view # locate button for next page
    driver.execute_script("window.scrollTo(0, 2100)")
    
    # click to go to next page
    WebDriverWait(driver, 20)\
        .until(EC.presence_of_element_located((By.XPATH, '//*[@id="technical-content-search"]/div[2]/button[2]')))\
        .click()
    time.sleep(np.random.randint(5,10))
    
driver.close()

