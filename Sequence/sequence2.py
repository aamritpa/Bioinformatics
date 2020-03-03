#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys
import time


# In[77]:


#Filename =sys.argv[1] # Read input file 
df = pd.read_excel('sequenceInput.xlsx')


data= pd.DataFrame()
data=df

def get_integer_value(newdata):
    return re.findall('\d+',newdata)


human_data=data[data['Organism'].str.contains("Homo sapiens")] #Selecting only organism which contains 'Homo sapiens'
data= human_data.reset_index(drop=True)


# In[2]:




#Calculating length of sequences
Sequence_Length=pd.Series(data['Sequence']).str.len()
position = pd.Series(len(data),dtype=np.str)
for i in range(len(data)):
    position[i]=''

for i in range(len(data)):
    for j in range((Sequence_Length[i])):
        if (data['Sequence'][i][j])=='C':
            position[i]=position[i]+','+str(j)


data['Sequence C position']= position.apply(get_integer_value)
count =pd.Series(data['Sequence C position'])
length=count.str.len()
data['Sequence C count']=length


pos = pd.Series(len(data),dtype=np.str)
for i in range(len(data)):
    pos[i]=''


for i in range(len(data)):
    for j in range(len(data['Sequence'][i])-2):
        if data['Sequence'][i][j]=='N' and data['Sequence'][i][j+1]!='P' and (data['Sequence'][i][j+2]=='S' or data['Sequence'][i][j+2]=='T'):
            pos[i]=pos[i]+str(j)+','
data['N_X_S/T']=pos.apply(get_integer_value)
data['Sequence N_X_S/T count']=pd.Series(data['N_X_S/T']).str.len()           

data.to_excel('sequenceOutput2.xlsx')


# In[3]:





# In[4]:





# In[ ]:




