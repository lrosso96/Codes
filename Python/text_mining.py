# ------------------------------------------- #
#       Fiscal reforms in Latin America       #
# ------------------------------------------- #

# code to do textmine art. iv from latin american countries
# author: Lucas Rosso
# date: 24-03-2021

# installing packages
# !pip install nltk
#pip install pdfminer.six
#pip install PyPDF2
#pip install pdfkit
# pip install fitz
# pip install PyMuPDF



from os import chdir
main_dir = 'C:/Users/LR/Desktop/ME/Ayudantía Wagner/Fiscal_Project/Python' 
raw_article_iv = '/raw_article_iv'
chdir(main_dir)

import PyPDF2
import mining_func


chdir(main_dir + raw_article_iv)
reader = PyPDF2.PdfFileReader(
    'CHL_18.pdf')

print(reader.documentInfo)

num_of_pages = reader.numPages
print('Number of pages: ' + str(num_of_pages))

reader = PyPDF2.PdfFileReader('CHL_18.pdf')
num_of_pages = reader.numPages # for the loop
reader.getPage(11-1).extractText().find('fiscal')

import nltk
from nltk.tokenize import sent_tokenize

text = reader.getPage(12).extractText()
tokenized_text=sent_tokenize(text)

import io
CHI_18 = open('CHL_18.pdf','r')

pdfFileObj = open('CHL_18.pdf','rb')     #'rb' for read binary mode
reader = PyPDF2.PdfFileReader(pdfFileObj)


import fitz 
doc = fitz.open("CHL_18.pdf")
num_pages = doc.pageCount
page = doc.loadPage(11).getText() 
text = doc.loadPage(11).get_text("html")

page.searchFor("fiscal")
doc.loadPage(11).searchFor("fiscal")

paa = []
para = doc.loadPage(0).get_text("blocks")
para = [x[4] for x in para]
paa.extend(para)


##################################
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
text =  []
iso3c = []
year  = [] 

# loop to extract all paragraphs (aka "blocks")
for file in files:
    doc = fitz.open(file)
    num_pages = doc.pageCount
    for page in range(0,num_pages):
        text_ = doc.loadPage(page).getText("blocks")
        text_ = [x[4] for x in text_] # 4 = element with string in tuple
        for te in text_:
            iso3c_ = file[0:3]
            iso3c.append(iso3c_)
            
            year_  = str(20) + file[4:6]
            year.append(year_)
        text.extend(text_)  

# keywords indicating tax policy
keywords = ['fiscal consolidation', 'tax reform', 'corporate tax'
            'value added tax', 'income tax', 'property tax', 'VAT', 'tax bill']

# Filtering the blocks that contain "keywords"
key_para = [x for x in text if
              any(y in x for y in keywords)]

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
})
    
df.to_csv(path_or_buf='key_paragraphs.csv',sep='|',index=False)



        
