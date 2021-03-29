# ------------------------------------------- #
#       Fiscal reforms in Latin America       #
# ------------------------------------------- #

# code to download and store art. IV from latin american countries
# author: Lucas Rosso
# date: 24-03-2021

# loading libraries
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
import time 
import numpy as np
import pandas as pd
import os
from os import chdir
chdir("C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python")

# function to take care of downloading file
def enable_download_headless(browser,download_dir):
    browser.command_executor._commands["send_command"] = ("POST", '/session/$sessionId/chromium/send_command')
    params = {'cmd':'Page.setDownloadBehavior', 'params': {'behavior': 'allow', 'downloadPath': download_dir}}
    browser.execute("send_command", params)
    
def Rename_file(new_name, Dl_path):  
    filename = max([f for f in os.listdir(Dl_path)])
    os.rename(os.path.join(Dl_path, filename), os.path.join(Dl_path, new_name+'.pdf'))


countries = ['Argentina', 'Bolivia', 'Brasil', 'Chile', 'Colombia',
             'Costa Rica', 'Ecuador', 'El Salvador','Guatemala', 
             'Honduras', 'Jamaica','México', 'Nicaragua','Panamá', 
             'Paraguay', 'Perú','República Dominicana',
             'Trinidad and Tobago','Uruguay', 'Venezuela']

iso3c = ['ARG', 'BOL', 'BRA', 'CHL', 'COL', 'CRI', 'ECU', 'SLV', 
         'GTM', 'HND', 'JAM', 'MEX', 'NIC', 'PAN', 'PRY', 'PER',
         'DOM', 'TTO', 'URY', 'VEN']

links = []
for iso in iso3c:
    links.append('https://www.imf.org/en/countries/' 
                    +iso +'?selectedfilters=Article%20IV%20Staff%20Reports#whatsnew')

# chrome webdriver options
options = Options()
options.page_load_strategy = 'normal'
#options.add_argument("--start-maximized")
options.add_argument('--disable-extensions')
options.add_argument("--disable-notifications")
options.add_argument('--no-sandbox')
options.add_argument('--verbose')
options.add_experimental_option("prefs", {
        "download.default_directory": "C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python/raw_article_iv",
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing_for_trusted_sources_enabled": False,
        "safebrowsing.enabled": False
})
options.add_argument('--disable-gpu')
options.add_argument('--disable-software-rasterizer')

driver = webdriver.Chrome(options=options) 

# directory to store art. iv files
download_dir = "C:\\Users\\LR\\Desktop\\ME\\Ayudantía Wagner\\Fiscal_Project\\Python\\raw_article_iv"

# function to handle setting up headless download
enable_download_headless(driver, download_dir)

for link in links:
    driver.get(link)
    #driver.get(links[5])
    time.sleep(np.random.randint(0,10))

    for i in range(1,15):
        try:
            # art. iv year (for name of the file)
            year = driver.find_element_by_xpath('//*[@id="docSearch_GUID"]/div/div[2]/div/div['+str(i)+']/p').text
            year = year[year.find(",")+2:year.find(",")+6]    
           
            WebDriverWait(driver, 10)\
                .until(EC.presence_of_element_located((By.XPATH, '//*[@id="docSearch_GUID"]/div/div[2]/div/div['+str(i)+']/h6/a')))\
                .click() 
            
            # countrycode
            countrycode = link[link.find("s/")+2:link.find("s/")+5]
            
            # art. iv year (for name of the file)
            # year = driver.find_element_by_xpath('/html/body/div[3]/main/article/div[1]/div/section[1]/h2').text
            # year = year[year.find(":")+1:year.find(":")+6]
             
            # download art. iv
            WebDriverWait(driver, 10)\
                .until(EC.presence_of_element_located((By.XPATH, '/html/body/div[3]/main/article/div[1]/div/section[1]/p[6]/a[1]')))\
                .click()
            
            # rename file with common structure    
            Rename_file(countrycode+'_'+str(year),download_dir)
            
            # back to main page
            driver.back()
            
            time.sleep(np.random.randint(0,3))
            
        except:
            try:
                # go to next page
                WebDriverWait(driver, 10)\
                    .until(EC.presence_of_element_located((By.XPATH, '//*[@id="docSearch_GUID"]/div/div[3]/div/a')))\
                    .click() 
                for i in range(1,15):
                    try:
                        # art. iv year (for name of the file)
                        year = driver.find_element_by_xpath('//*[@id="docSearch_GUID"]/div/div[2]/div/div['+str(i)+']/p').text
                        year = year[year.find(",")+2:year.find(",")+6]    
           
                        WebDriverWait(driver, 10)\
                            .until(EC.presence_of_element_located((By.XPATH, '//*[@id="docSearch_GUID"]/div/div[2]/div/div['+str(i)+']/h6/a')))\
                            .click() 

                        # art. iv year (for name of the file)
                        # year = driver.find_element_by_xpath('/html/body/div[3]/main/article/div[1]/div/section[1]/h2').text
                        # year = year[year.find(":")+1:year.find(":")+6]
             
                        # download art. iv
                        WebDriverWait(driver, 10)\
                            .until(EC.presence_of_element_located((By.XPATH, '/html/body/div[3]/main/article/div[1]/div/section[1]/p[6]/a[1]')))\
                            .click()
            
                        # rename file with common structure    
                        Rename_file(countrycode+str(year),download_dir)
            
                        # back to main page
                        driver.back()
                    except IndexError:
                        break
            except:
                print('on to the next')
                pass