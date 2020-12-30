# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 18:08:01 2020

@author: LR
"""

# installing packages
#pip install requests_html
#pip install selenium

from os import chdir
chdir ("C:/Users/LR/Desktop/ME/Programación/Python")

from bs4 import BeautifulSoup
import requests

url = "https://www.timeanddate.com/holidays/"
page = requests.get(url, 'html.parser')

soup = BeautifulSoup(page.text)

print(soup.find_all("ul", class_='category-list__list'))

ul=BeautifulSoup(str(soup.find_all("ul", class_='category-list__list')))

# we can display them with a loop
for i in ul.find_all("li"):
    print(i.find("a").get_text())
    
# we can store them in a loop
COUNTRIES=[]
for i in ul.find_all("li"):
    print(i.find("a").get_text())
    COUNTRIES.append(i.find("a").get_text())
print(COUNTRIES)

# generating the list of urls to scrape information
URLS=[]
for i in ul.find_all("li"):
    for year in range(2000, 2019):
        print(i.find("a").get_text())
        URLS.append(i.find("a").get("href")+ str(year) + '?hol=9')
print(URLS)
print(URLS[0])

# Example: Chile 2000
from selenium import webdriver
URL="https://www.timeanddate.com%s"%URLS[741]
browser = webdriver.Chrome()

import time

browser.get(URL)
time.sleep(3)
html = browser.page_source
soup = BeautifulSoup(html, "lxml")

browser.quit()

table=soup.find("table", class_="table table--left table--inner-borders-rows table--full-width table--sticky table--holidaycountry").tbody

date = [] # date of the holiday day, month
dow  = [] # day of the week
name = [] # name of the holiday
det  = [] # type (e.g. national holiday)
for i in range(1,100):
    try:
        date_ = table.find_all("tr", class_="showrow")[i].find("th", class_="nw")
        date.append(date_.text)
        
        dow_ = table.find_all("tr", class_="showrow")[i].find("td", class_="nw")
        dow.append(dow_.text)
        
        name_ = table.find_all("tr", class_="showrow")[i].find_all("td")[1]
        name.append(name_.text)
        
        det_ = table.find_all("tr", class_="showrow")[i].find_all("td")[2]
        det.append(det_.text)
    except IndexError:
        break
    
import pandas as pd
df = pd.DataFrame({'date':date,'dow':dow,'name':name,'type':det})
df.to_csv(path_or_buf='data.csv',na_rep='.',sep=',',index=False)
print(df)

#final_data = []
date = [] # date of the holiday day, month
dow  = [] # day of the week
name = [] # name of the holiday
det  = [] # type (e.g. national holiday)
year = [] 
country = [] 
# Now for all countries in every year
for i in URLS:
    URL="https://www.timeanddate.com%s"%i
    print(URL)
    browser = webdriver.Chrome()
    browser.get(URL)
    time.sleep(30)
    html = browser.page_source
    soup = BeautifulSoup(html, "lxml")
    browser.quit()
    try:
        table=soup.find("table", class_="table table--left table--inner-borders-rows table--full-width table--sticky table--holidaycountry").tbody

        for j in range(1,365):
            try:
                date_ = table.find_all("tr", class_="showrow")[j].find("th", class_="nw")
                date.append(date_.text)
        
                dow_ = table.find_all("tr", class_="showrow")[j].find("td", class_="nw")
                dow.append(dow_.text)
        
                name_ = table.find_all("tr", class_="showrow")[j].find_all("td")[1]
                name.append(name_.text)
        
                det_ = table.find_all("tr", class_="showrow")[j].find_all("td")[2]
                det.append(det_.text)
                
                year.append(URL[-10:-6]) # URL always ends with '/year?hol=9'
                
                start = URL.find("s/") # path always of the form '/holidays/country/'
                end   = URL.find("/2")
                country.append(URL[start+2:end])
            except IndexError:
                break
        
        #final_data.append(date, dow, name, det)    
        #df = pd.DataFrame({'date':date,'dow':dow,'name':name,'type':det}, columns=None)
        #final_data.append(df)
    except: 
        print(soup.find("section", class_="table-data__table").find("p").text)
        pass      
    #DATA = pd.DataFrame(final_data)

df = pd.DataFrame({'country':country,'year':year,'date':date,'dow':dow,'name':name,'type':det})    
df.to_csv(path_or_buf='data.csv',na_rep='.',sep=',',index=False)
