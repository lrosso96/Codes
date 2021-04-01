# ------------------------------------------- #
#       Fiscal reforms in Latin America       #
# ------------------------------------------- #

# code to download and store art. IV from latin american countries
# author: Lucas Rosso
# date: 24-03-2021

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
chdir("C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python/raw_article_ivs")

# function to take care of downloading file
def enable_download_headless(browser,download_dir):
    browser.command_executor._commands["send_command"] = ("POST", '/session/$sessionId/chromium/send_command')
    params = {'cmd':'Page.setDownloadBehavior', 'params': {'behavior': 'allow', 'downloadPath': download_dir}}
    browser.execute("send_command", params)

# function to rename file once downloaded
def Rename_file(new_name, Dl_path):  
    filename = max([f for f in os.listdir(Dl_path)], key=os.path.getctime)
    os.rename(os.path.join(Dl_path, filename), os.path.join(Dl_path, new_name+'.pdf'))


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

# %%

driver = webdriver.Chrome(options=options) 

# directory to store art. iv files
download_dir = "C:\\Users\\LR\\Desktop\\ME\\Ayudantía Wagner\\Fiscal_Project\\Python\\raw_article_ivs"

# function to handle setting up headless download
enable_download_headless(driver, download_dir)

# emptying the directory (except the chromedriver file)
files = os.listdir(download_dir) 
for file in files:
    if file != 'chromedriver.exe':
        os.remove(download_dir +'\\' + file)
    else:
        pass


# some countries appear on the wrong countrypage
not_wanted = ['Republic of Croatia', 'Republic of the Marshall Islands',
              'Thailand']

for link in links:
    driver.get(link)
    time.sleep(np.random.randint(0,10))

    #elements in page
    elements = driver.find_elements_by_class_name('result-item')
    for i in range(1,len(elements)+1):             
            # click to get art iv (and pass for "not wanted" countries)
            country_check = driver.find_element_by_xpath('//*[@id="docSearch_GUID"]/div/div[2]/div/div['+str(i)+']/h6/a').text
            country_check = country_check[0:country_check.find(":")-1]
            if any(y in country_check for y in not_wanted):
                continue
            else:
                WebDriverWait(driver, 10)\
                    .until(EC.presence_of_element_located((By.XPATH, '//*[@id="docSearch_GUID"]/div/div[2]/div/div['+str(i)+']/h6/a')))\
                    .click() 
            
            # art. iv year (for name of the file)
            try:
                year = driver.find_element_by_xpath('/html/body/div[3]/main/article/div[1]/div/section[1]/p[4]').text
                year = year[year.find(",")+2:year.find(",")+6] 
                ind = 1
            except:
                year = driver.find_element_by_xpath('/html/body/div[3]/main/article/div[2]/div/section[1]/p[4]').text
                year = year[year.find(",")+2:year.find(",")+6] 
                ind = 2
            
            time.sleep(np.random.randint(1,3))
            
            # countrycode (for name of the file)
            countrycode = link[link.find("s/")+2:link.find("s/")+5]
             
            # download art. iv
            try:
                WebDriverWait(driver, 10)\
                    .until(EC.presence_of_element_located((By.XPATH, '/html/body/div[3]/main/article/div['+str(ind)+']/div/section[1]/p[6]/a[1]')))\
                    .click()
            except:
                WebDriverWait(driver, 10)\
                    .until(EC.presence_of_element_located((By.XPATH, '/html/body/div[3]/main/article/div['+str(ind)+']/div/section[1]/p[4]/a[1]')))\
                    .click()
                
            time.sleep(np.random.randint(4,7))
            
            # back to main page
            driver.back()
            
            # rename file with common structure 
            try:
                Rename_file(countrycode+'_'+str(year),download_dir)
            except: 
                pass
            
            time.sleep(np.random.randint(4,7))
    
    # try to loop again for elements in second page
    try:
        # go to next page
        WebDriverWait(driver, 15)\
            .until(EC.presence_of_element_located((By.XPATH, '//*[@id="docSearch_GUID"]/div/div[3]/div/a')))\
            .click() 
                        
        time.sleep(np.random.randint(4,7))
        
        #elements in page
        elements = driver.find_elements_by_class_name('result-item')
        
        for i in range(1,len(elements)+1):               
            # click to get art iv (and pass for "not wanted" countries)
            country_check = driver.find_element_by_xpath('//*[@id="docSearch_GUID"]/div/div[2]/div/div['+str(i)+']/h6/a').text
            country_check = country_check[0:country_check.find(":")-1]
            if any(y in country_check for y in not_wanted):
                continue
            else:
                WebDriverWait(driver, 10)\
                    .until(EC.presence_of_element_located((By.XPATH, '//*[@id="docSearch_GUID"]/div/div[2]/div/div['+str(i)+']/h6/a')))\
                    .click()
            
            # art. iv year (for name of the file) 
            year = driver.find_element_by_xpath('/html/body/div[3]/main/article/div[1]/div/section[1]/p[4]').text
            year = year[year.find(",")+2:year.find(",")+6]
            
            time.sleep(np.random.randint(1,3))
            
            # download art. iv                                   
            WebDriverWait(driver, 10)\
                .until(EC.presence_of_element_located((By.XPATH, '/html/body/div[3]/main/article/div[1]/div/section[1]/p[6]/a[1]')))\
                .click()
                            
            time.sleep(np.random.randint(5,8))
            
            # rename file with common structure    
            Rename_file(countrycode+'_'+str(year),download_dir)
            
            # back to main page
            driver.get(link)
            time.sleep(np.random.randint(5,8))
              
            # go to next page (browser goes back to page 1)
            WebDriverWait(driver, 10)\
                .until(EC.presence_of_element_located((By.XPATH, '//*[@id="docSearch_GUID"]/div/div[1]/div/span/a')))\
                .click()
    except:
        # click "no feedback"
        # try:
        #     WebDriverWait(driver, 10)\
        #         .until(EC.presence_of_element_located((By.XPATH, '//*[@id="fsrInvite"]/section[3]/button[2]')))\
        #         .click()
                
        time.sleep(np.random.randint(3,6))
    
driver.close()
