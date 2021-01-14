# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:57:39 2021

@author: LR
"""

# installing packages
#pip install requests_html
#pip install selenium

# required libraries
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from os import chdir
import time
import pandas as pd

# setting wd for web driver
chdir ("C:/Users/LR/Desktop/ME/Ayudantía Wagner/Holidays and Growth/Python")

# chrome webdriver options
options = Options()
options.page_load_strategy = 'normal'
options.add_argument("--start-maximized")
options.add_argument('--disable-extensions')
driver = webdriver.Chrome(options=options) # web driver must be on the current wd

# pre-allocation
working_days = []
country      = []
year         = []
weekend1     = []
weekend2     = []
holidays     = []
    
driver = webdriver.Chrome(options=options)

url = "https://www.timeanddate.com/date/workdays.html?d1=01&m1=01&y1=2019&d2=31&m2=12&y2=2019&ti=on&"

driver.get(url)
    
# link = driver.find_element_by_link_text("Change Country")
# link.send_keys("\n")

# element = WebDriverWait(driver, 10).until(
#     EC.presence_of_element_located((By.ID, "country"))
#     )
# element.click()   

# all_options = element.find_elements_by_tag_name("option")

# 8: Afganistan, 239: Zimbabwe    
for o in range(8,240):

    WebDriverWait(driver, 10)\
        .until(EC.presence_of_element_located((By.ID, "chco2")))\
        .send_keys("\n")

    WebDriverWait(driver, 10)\
        .until(EC.presence_of_element_located((By.ID, "country")))\
        .click() 
            
    # year 
    year_aux = url[-11:-7]
    year.append(year_aux)
    print(year_aux)
            
    # country 
    time.sleep(10) 
    # xpath = "//*[@id='country']/option[%o]" %(o)
    xpath = "//*[@id='country']/option[" + str(o) + "]"
    driver.find_element_by_xpath(xpath).click()
    country_aux = driver.find_element_by_xpath(xpath).text
    country.append(country_aux)
    print(country_aux)
        
            
    WebDriverWait(driver, 10)\
        .until(EC.presence_of_element_located((By.ID, "tzq_save")))\
        .click()
               
    time.sleep(5) 
    
    # working days    
    try:
        days = driver.find_element_by_xpath("//*[@id='weekday_resall']/div[1]/div/div[1]/h2").text
        days = days[-8:-5]
        working_days.append(days)
    except:
        days = driver.find_element_by_xpath("//*[@id='weekday_resall']/div[2]/div/div[1]/h2").text
        days = days[-8:-5]
        working_days.append(days)
        
    time.sleep(2)
    
    # first weekend day    
    try:
        weekend1_aux = driver.find_element_by_xpath("//*[@id='weekday_resall']/div[1]/div/div[2]/h4[1]").text
        aux1 = str.split(weekend1_aux)
        weekend1_ = aux1[2]
        weekend1.append(weekend1_)
    except:
        weekend1_aux = driver.find_element_by_xpath("//*[@id='weekday_resall']/div[2]/div/div[2]/h4[1]").text
        aux1 = str.split(weekend1_aux)
        weekend1_ = aux1[2]
        weekend1.append(weekend1_)
        
    time.sleep(2)
    
    # second weekend day    
    try:
        weekend2_aux = driver.find_element_by_xpath("//*[@id='weekday_resall']/div[1]/div/div[2]/h4[2]").text
        aux2 = str.split(weekend2_aux)
        weekend2_ = aux2[2]
        if weekend2_ == "holidays:":
            weekend2_ = aux1[2]
        weekend2.append(weekend2_) 
    except:
        try:
            weekend2_aux = driver.find_element_by_xpath("//*[@id='weekday_resall']/div[2]/div/div[2]/h4[2]").text
            aux2 = str.split(weekend2_aux)
            weekend2_ = aux2[2]
            weekend2.append(weekend2_) # only one weekend day
        except:
            weekend2_ = aux1[2] # only one weekend day
            weekend2.append(weekend2_)
    
    # holidays                
    # try:
    #     holidays_aux = driver.find_element_by_xpath("//*[@id='weekday_resall']/div[1]/div/div[2]/h4[3]").text
    #     aux3 = str.split(holidays_aux)
    #     holidays_ = aux3[1]
    #     holidays.append(holidays_)
    # except:
    #     try:
    #         holidays_aux = driver.find_element_by_xpath("//*[@id='weekday_resall']/div[1]/div/div[2]/h4[2]").text
    #         aux3 = str.split(holidays_aux)
    #         holidays_ = aux3[1]
    #         holidays.append(holidays_)        
    #     except:
    #         holidays_ = "."
    #         holidays.append(holidays_)
            
    time.sleep(5)
        
driver.close()
        
df = pd.DataFrame({
    'country':country,
    'year':year,
    'working_days':working_days,
    'weekend1':weekend1,
    'weekend2':weekend2,
})

print(df)

df.to_csv(path_or_buf='workweek.csv',na_rep='.',sep=',',index=False)

# print(all_options)

# print(country_aux)


# print(option)
# print(country)
# print(working_days)
# print(weekend2)

# print(country_list)
