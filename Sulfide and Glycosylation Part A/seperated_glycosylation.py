#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys
import time


#Input file must be given to run the code.
#For Example Run(in terminal) 'python3 Bioinformatics input.xlsx'

Filename =sys.argv[1] # Read input file 

df = pd.read_excel(Filename)# Read EXCEL Files and storing in database 'df'
data= pd.DataFrame()# Creating New Data Frame 
glyco_sulfide_frame= pd.DataFrame()   #Data Frame which contains Glycosilation and sulfide.


# In[368]:


# Storing Usefull information From the database 'df' 
data['Entry']=df['Entry'] #Storing Entry
data['Entry name']=df['Entry name']  #Storing Entry Name
data['Protein names']=df['Protein names']  #Storing Protein Names
data['Gene names']=df['Gene names'] # Storing Gene Names
data['Length']=df['Length'] # Storing Length
data['Organism']=df['Organism'] # Storing Type of Organism
data['Mass']=df['Mass'] #Storing Mass
data['Status']=df['Status']

#Now removing the extra data we dont need and containg only organism 'Human'
human_data=data[data['Organism'].str.contains("Homo sapiens")] #Selecting only organism which contains 'Homo sapiens'
data= human_data

data['Disulfide bond']=df['Disulfide bond'] # Storing Disulfide bond 
data['Glycosylation']= df['Glycosylation']# Storing Glycosylation Position
data=data.dropna().reset_index(drop=True) # Drop all rows which contains Not a number and reset the index 

# Working on getting the relative positions of the Disulfide bond
temp_disulfide= data['Disulfide bond'] # Storing data temporary as 'disulfide_column'
disulfide_column= data['Disulfide bond']

def get_sulfide_value(newdata):
    return re.findall('\d+\..\d+',newdata)

bond= temp_disulfide.apply(get_sulfide_value) # Function call which gives all the positons 
glyco_sulfide_frame['bond']=bond
data['Disulfide bond']=disulfide_column.apply(get_sulfide_value)# Storing extracted disulfide value in data Frame.


def get_int_value(data):        #function defined for global purpose use to get integer values for any column.
    return re.findall('\d+',data)
def get_float_value(data):
    return re.findall('\d+\.\d+',data)


# In[369]:


#CARBOHYD 110;  /note="N-linked (GlcNAc...) asparagine'
#CARBOHYD 36; /note="N-linked (Glc) (glycation...
# In[32]:
def get_Nlinked_GlcNAc(data):
    return re.findall('\d+;  /note="N-linked \(GlcNAc...\) asparagine"',data)

def get_Nlinked_Glc(data):
    return re.findall('\d+;  /note="N-linked \(Glc\)',data)

def get_OLinked(data):
    return re.findall('\d+\;  /note=\"O-linked',data)
    
data1= data['Glycosylation'] # Storing data temporary as 'Glyco_data'
data['Glycosylation_GlcNac']= data1.apply(get_Nlinked_GlcNAc) # the return data with 'N-linked (GlcNAc...)'

data2= data['Glycosylation']
data['Glycosylation_Nlinked_Glc']=data2.apply(get_Nlinked_Glc)

data3= data['Glycosylation']
data['Glycosylation_Olinked']=data3.apply(get_OLinked)

#Now Removing the extra word 'CARBOHYD' and getting all positions of Glycosylation.
temp_data= data['Glycosylation_GlcNac'].astype(str)
data['Glycosylation_GlcNac']= temp_data.apply(get_int_value) 

temp_data= data['Glycosylation_Nlinked_Glc'].astype(str)
data['Glycosylation_Nlinked_Glc']= temp_data.apply(get_int_value) 

temp_data= data['Glycosylation_Olinked'].astype(str)
data['Glycosylation_Olinked']= temp_data.apply(get_int_value) 


# In[370]:


data.to_excel('seperated_glycosylation.xlsx')


# In[ ]:





# In[ ]:





# In[ ]:




