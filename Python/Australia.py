# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:57:39 2021

@author: LR
"""

# installing packages
#pip install requests_html
#pip install selenium

# required libraries
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
import time
from os import chdir

# setting wd for web driver
chdir ("C:/Users/LR/Desktop/ME/Ayudantía Wagner/Holidays and Growth/Python")

# year list (to update url)
year = list(range(2000,2020))

url_a = "https://www.timeanddate.com/date/workdays.html?d1=01&m1=01&y1="
url_b = "&d2=31&m2=12&y2="

URL = []
for y in year:
    url_aux = url_a + str(y) + url_b + str(y) + "&ti=on&"
    URL.append(url_aux)
print(URL)

options = webdriver.ChromeOptions()
options.add_argument("--start-maximized")
options.add_argument('--disable-extensions')
driver = webdriver.Chrome(chrome_options=options) # web driver must be on the current wd

url = URL[0]
driver.get(url)
link = driver.find_element_by_link_text("Change Country")
link.send_keys("\n")

try:
    element = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, "country"))
        )
    element.click()   
except:
    driver.quit
    
all_options = element.find_elements_by_tag_name("option")   

country = all_options[0].text
print(country)

ele = all_options[0].click()

WebDriverWait(driver, 10)\
    .until(EC.presence_of_element_located((By.ID, "tzq_save")))\
        .click()
        
time.sleep(5)
            
WebDriverWait(driver, 10)\
    .until(EC.presence_of_element_located((By.ID, "ti")))\
        .click() 
        
html = driver.page_source
soup = BeautifulSoup(html, "lxml")

driver.quit()
days = soup.find("div", class_ = "re-result five columns").h2.text
days = days[-9:-5]
print(days)

year_aux = url[-11:-7]
print(year_aux)

weekend = soup.find("div", class_ = "re-details five columns offset-2") 

weekend1 = weekend.find_all("h4")[0].text
aux = str.split(weekend1)
weekend1 = aux[2]   
print(weekend1)    

weekend2 = weekend.find_all("h4")[1].text
aux = str.split(weekend2)
weekend2 = aux[2]       
print(weekend2)    

print(aux)
