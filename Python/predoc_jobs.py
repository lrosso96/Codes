# ------------------------------------ #
#           Pre-Doc Openings           #
# ------------------------------------ #

# Created by: Lucas Rosso
# Created on: 19/02/2021
# any bugs please report to lrosso@fen.uchile.cl
# closely linked to "adopt a pet" from Alvaro Carril (https://acarril.github.io/posts/adopt-dog-python)

from os import chdir
chdir ("C:/Users/LR/Desktop/ME/Programación/Python") #must be changed

from bs4 import BeautifulSoup
import requests
import pandas as pd
import datetime
import sys

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

# Not-NBER positions
url = "https://www.nber.org/career-resources/research-assistant-positions-not-nber"
page = requests.get(url, 'html.parser')

soup = BeautifulSoup(page.text, features="lxml")

jobs = soup.find('div', class_='page-header__intro-inner').find_all('p')
jobs = jobs[1:len(jobs)-2] # delete introduction and closing paragraphs

id          = []
pos_name    = []
researcher  = []
institution = []
field       = []
link        = []
for job in jobs:    
    # position name
    pos_name_ = job.text.split('\n')[0]
    pos_name.append(pos_name_)
        
    # researcher
    researcher_ = job.text.split("\n")[1].split(":")[1]
    researcher.append(researcher_)
        
    # institution
    institution_ = job.text.split("\n")[2].split(":")[1]
    institution.append(institution_)
        
    # field
    field_ = job.text.split("\n")[3].split(":")[1]
    field.append(field_)
        
    # job posting link
    link_ = job.find('a').get('href')
    link.append(link_)
    
    # id
    id_ = int(str(len(pos_name_)) + str(len(researcher_)) + str(len(institution_)))
    id.append(id_)
    
# Create dataframe with newly collected data 
df = pd.DataFrame({
    'id': id,
    'pos_name': pos_name,
    'researcher': researcher,
    'institution': institution,
    'field': field,
    'link': link,
})

# 'old' dataset with jobs. must be runned onced to be stored in cd and then compared to new
# df.to_csv(path_or_buf='current_jobs.csv',na_rep='.',sep=';',index=False)

# read 'old' dataset for comparison
df_old = pd.read_csv('current_jobs.csv',sep=';',index_col=False)

# check differences between datasets (and exit if there is no diff)
# diffs = sum((df.id != df_old.id))
diffs = list(set(df.id) ^ set(df_old.id))

# report execution with no updates and exit
if diffs == []:
    now = datetime.datetime.now()
    date_time = now.strftime('%m/%d/%Y %H:%M:%S')
    sys.stdout.write(date_time + '\n')
    print('no new openings :(')
    sys.exit()    

# compute sets of added and removed jobs, create corresponding dataframes (email)
added = list(set(df.id) - set(df_old.id))
gone  = list(set(df_old.id) - set(df.id))
dfAdded  = df[df['id'].isin(added)]
dfGone = df_old[df_old['id'].isin(gone)]



# send new jobs

# subject
if len(added) > 0:
    subject = 'New Opening Alert!'
else:
    subject = 'Closed positions.'

# sender
sender = 'your_email@gmail.com'
password = 'your_password'
smtp_server = "smtp.gmail.com"
port = 587

# recepient
recipient  = 'your_email@gmail.com'

msg = MIMEMultipart()
msg['From'] = sender
msg['To'] = recipient
msg['Subject'] = subject  


message = f"""\

Check changes in RA openings below 

Opened:
{dfAdded[['pos_name', 'researcher', 'institution', 'field', 'link']].stack().to_string(index=False)}


Closed:
{dfGone[['pos_name', 'researcher', 'institution', 'field', 'link']].stack().to_string(index=False)}
"""
msg.attach(MIMEText(message,'plain')) 

smtpObj = smtplib.SMTP(smtp_server, port)
smtpObj.ehlo()
smtpObj.starttls()
smtpObj.login(sender, password)
smtpObj.sendmail(sender, recipient, msg.as_string())
smtpObj.quit()

fileName = 'current_jobs.csv'
df.to_csv(fileName,na_rep='.',sep=';',index=False)

# report execution with updates
now = datetime.datetime.now()
date_time = now.strftime("%m/%d/%Y %H:%M:%S")
sys.stdout.write(date_time + ' Run with updates!\n')

