# ------------------------------------------- #
#       Fiscal reforms in Latin America       #
# ------------------------------------------- #

# code to download worldwide tax guides available from EY
# author: Lucas Rosso
# date: 08-04-2021

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
chdir("C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python/EY_tax_guides")

# function to take care of downloading file
def enable_download_headless(browser,download_dir):
    browser.command_executor._commands["send_command"] = ("POST", '/session/$sessionId/chromium/send_command')
    params = {'cmd':'Page.setDownloadBehavior', 'params': {'behavior': 'allow', 'downloadPath': download_dir}}
    browser.execute("send_command", params)

# function to rename file once downloaded
def Rename_file(new_name, Dl_path):  
    filename = max([f for f in os.listdir(Dl_path)], key=os.path.getctime)
    os.rename(os.path.join(Dl_path, filename), os.path.join(Dl_path, new_name+'.pdf'))

# chrome webdriver options
options = Options()
options.page_load_strategy = 'normal'
#options.add_argument("--start-maximized")
options.add_argument('--disable-extensions')
options.add_argument("--disable-notifications")
options.add_argument('--no-sandbox')
options.add_argument('--verbose')
options.add_experimental_option("prefs", {
        "download.default_directory": "C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python/EY_tax_guides",
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing_for_trusted_sources_enabled": False,
        "safebrowsing.enabled": False,
        "plugins.always_open_pdf_externally": True
})
options.add_argument('--disable-gpu')
options.add_argument('--disable-software-rasterizer')

# %%
# scraping the webpage and downloading files
url = 'https://www.ey.com/en_gl/tax-guides/tax-guide-library-archive'
driver = webdriver.Chrome(options=options) 

# directory to store art. iv files
download_dir = "C:\\Users\\LR\\Desktop\\ME\\Ayudantía Wagner\\Fiscal_Project\\Python\\EY_tax_guides"

# function to handle setting up headless download
enable_download_headless(driver, download_dir)

# emptying the directory (except the chromedriver file)
files = os.listdir(download_dir) 
for file in files:
    if file != 'chromedriver.exe':
        os.remove(download_dir +'\\' + file)
    else:
        pass

# opening the link    
driver.get(url)

# clicking in "I decline optional cookies"
try:
    WebDriverWait(driver, 10)\
        .until(EC.presence_of_element_located((By.XPATH, '//*[@id="cookiePolicy"]/div/div[2]/button[1]')))\
        .click()
except:
    pass
time.sleep(np.random.randint(0,5))

# loop for clicking specific tax guide types
for i in range(1,9):
    WebDriverWait(driver, 10)\
        .until(EC.presence_of_element_located((By.XPATH, '//*[@id="body-components-wrapper"]/div/div[1]/div/div[2]/div/ul/li['+str(i)+']/header/h3/a')))\
        .click()  
    
    time.sleep(np.random.randint(3,5))
    # number of files to download
    num_files = driver.find_elements_by_class_name('richText-content')
    for j in range(1,len(num_files)+1):
        #click to download pdf
        WebDriverWait(driver, 30)\
            .until(EC.presence_of_element_located((By.XPATH, '//*[@id="accordion-content-02094958150-'+str(i-1)+'"]/div/div/div/div/div/p['+str(j)+']/a')))\
            .click() 
        time.sleep(np.random.randint(10,20))
    
    time.sleep(10)
    
    # close detail for tax guide type already done
    WebDriverWait(driver, 10)\
        .until(EC.presence_of_element_located((By.XPATH, '//*[@id="body-components-wrapper"]/div/div[1]/div/div[2]/div/ul/li['+str(i)+']/header/h3/a')))\
        .click()
        
    time.sleep(np.random.randint(20,30))

driver.close()
   
