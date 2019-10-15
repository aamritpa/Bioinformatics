#!/usr/bin/env python
# coding: utf-8

# In[185]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys


# In[186]:


#Input file must be given to run the code.
#For Example Run(in terminal) 'python3 Bioinformatics input.xlsx'

Filename =sys.argv[1] # Read input file 

df = pd.read_excel('input.xlsx')# Read EXCEL Files and storing in database 'df'
data= pd.DataFrame()# Creating New Data Frame 


# In[187]:


# Storing Usefull information From the database 'df' 
data['Entry']=df['Entry'] #Storing Entry
data['Entry name']=df['Entry name']  #Storing Entry Name
data['Protein names']=df['Protein names']  #Storing Protein Names
data['Gene names']=df['Gene names'] # Storing Gene Names
data['Length']=df['Length'] # Storing Length


# In[188]:


data['Organism']=df['Organism'] # Storing Type of Organism
data['Mass']=df['Mass'] #Storing Mass
 


# In[189]:


#Now removing the extra data we dont need and containg only organism 'Human'
columns=data[data['Organism'].str.contains("Human")] #Selecting only organism which contains HUMANS
data= columns


# In[190]:


data['Disulfide bond']=df['Disulfide bond'] # Storing Disulfide bond 
data['Glycosylation']= df['Glycosylation']# Storing Glycosylation Position


# In[191]:


data=data.dropna().reset_index(drop=True) # Drop all rows which contains Not a number and reset the index 


# In[192]:


# Working on getting the relative positions of the Disulfide bond
disulfide_column= data['Disulfide bond'] # Storing data temporary as 'disulfide_column'
def getdata(newdata):
    return re.findall('\d+ \d+',newdata)

data['Disulfide bond']= disulfide_column.apply(getdata) # Function call which gives all the positons 


# In[193]:


#Working on getting the positions of N-linked Glycosylation
Glyco_data= data['Glycosylation'] # Storing data temporary as 'Glyco_data'

def getGlycoNLinked(data):   #Making a Function to get the positions
    return re.findall('\d+\s* N-linked',data) # The return data will contain 'position and N-linked pattern'
data['Glycosylation']= Glyco_data.apply(getGlycoNLinked) # the return data will be like '49 N-linked'


# In[ ]:


#Now Removing the extra word 'N-linked' and getting all positions 
temp_data= data['Glycosylation'].astype(str)
def getValue(data):
    return re.findall('\d+',data) # return all positions
data['Glycosylation']= temp_data.apply(getValue) 


# In[199]:


data.to_excel('output.xlsx') # Output an Excel File with refined data.


# In[ ]:




