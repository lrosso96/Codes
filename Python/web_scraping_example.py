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

print(len(soup.find_all("table")))
table=soup.find("table", class_="table table--left table--inner-borders-rows table--full-width table--sticky table--holidaycountry").tbody
print(table)

table2 = table.select("[class=showrow]")
print(table2)

browser.quit()

Chile = []
for tr_tag in table2:
    print(tr_tag.text)
    Chile.append(tr_tag.text)
print(Chile)

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

print(df)