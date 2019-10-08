#!/usr/bin/env python
# coding: utf-8

# In[122]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys


# In[136]:


Filename =sys.argv[1] #Input file from User

df = pd.read_excel(Filename)#  r'/home/amritpal/Desktop/uniprot-human-filtered-reviewed:yes.xlsx') #Read EXCEL Files
new_data= pd.DataFrame()  #NEW DATA FRAME WITH EMPTY DATA


# In[137]:


# Storing all information
new_data['Entry']=df['Entry']
new_data['Entry name']=df['Entry name']
new_data['Protein names']=df['Protein names']
new_data['Gene names']=df['Gene names']
new_data['Length']=df['Length']


# In[138]:


new_data['Organism']=df['Organism']
new_data['Mass']=df['Mass']
columns=new_data[new_data['Organism'].str.contains("Human")]
new_data= columns


# In[140]:


new_data['Disulfide bond']=df['Disulfide bond'] #STORE 'ENTRY' COLUMN IN NEW DATA
new_data['Position']= df['Glycosylation'].str.extract(r'(\d+)')[0] # GET THE GLYCOSYLATION VALUE


# In[141]:


new_data=new_data.dropna().reset_index(drop=True) # DROP ALL THOSE ROWS WHOSE DATA IS NOT REQUIRED


# In[142]:


Tem_data= new_data['Disulfide bond']

def getdata(data):
    return re.findall('\d+ \d+',data)

new_data['Disulfide bond']= Tem_data.apply(getdata) # get postions for all 

print(new_data)


# In[114]:


new_data.to_excel('output.xlsx')


# In[ ]:





# In[ ]:





# In[ ]:




