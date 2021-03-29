# ------------------------------------------- #
#       Fiscal reforms in Latin America       #
# ------------------------------------------- #

# code to do textmine art. ivs from latin american countries
# author: Lucas Rosso
# date: 24-03-2021

#%% 

# pip install fitz

import os
import pandas as pd
from os import chdir
main_dir = 'C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python'
raw_article_iv = '/raw_article_iv' 
data = '/Data'
chdir(main_dir + raw_article_iv)

import fitz
files = os.listdir(main_dir + raw_article_iv) # art. ivs

# pre-allocation of vars
text      = []
iso3c     = []
year      = [] 

# loop to extract all paragraphs (aka "blocks")
for file in files:
    doc = fitz.open(file)
    num_pages = doc.pageCount
    for page in range(0,num_pages):
        text_ = doc.loadPage(page).getText("blocks")
        text_ = [x[4] for x in text_] # 4 = element with string in tuple
        text_aux = [t.split('\n \n') for t in text_]
        flat_list = [item for sublist in text_aux for item in sublist]
        for te in flat_list:
            iso3c_ = file[0:3]
            iso3c.append(iso3c_)
            
            year_  = file[4:8]
            year.append(year_)
        text.extend(flat_list)  

# keywords indicating tax policy
keywords = ['fiscal consolidation', 'tax reform', 'corporate tax'
            'value added tax', 'income tax', 'property tax', 'VAT', 'tax bill']

# Filtering the blocks that contain "keywords"
key_para = [x for x in text if
              any(y in x for y in keywords)]


# generating variable with list of keywords that appear on selected paragraph
strip_keywords = []
for k in keywords:
    strip_keywords.append(k.replace(" ", ""))
    
found_words = []
for key in key_para:
    found_words_ = [word for word in strip_keywords 
                    if word in key.replace(" ", "")]
    found_words.append(found_words_)
    
for word in found_words:
    for l in range(0,len(word)):
        for k in keywords:
            if word[l] == k.replace(" ", ""):
                word[l] = k
            else:
                pass

# indicators for year and country code
inds = []
for te in text:
    inds_ = (any(y in te for y in keywords))
    inds.append(inds_)
    
key_year  = [] 
key_iso3c = []
for ind in range(0,len(inds)):
    if (inds[ind] == 1):
        key_year.append(year[ind])
        key_iso3c.append(iso3c[ind])
        
# export key paragraphs as csv
chdir(main_dir + data)
df = pd.DataFrame({
    'iso3c':key_iso3c,
    'year':key_year,
    'key_paragraph':key_para,
    'keyword':found_words,
})
    
df.to_csv('key_paragraphs.csv',index=False,encoding='utf-16')
