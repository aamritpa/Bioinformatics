#!/usr/bin/env python
# coding: utf-8

# In[12]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys
import time


# In[13]:


Filename =sys.argv[1] # Read input file 
df = pd.read_excel(Filename)


# In[14]:


df=df.dropna().reset_index(drop=True)
data= pd.DataFrame()
data=df


# In[ ]:





# In[15]:


sulfide=df['Disulfide bond']
def get_sulfide_value(newdata):
    return re.findall('\d+\..\d+',newdata)

values= sulfide.apply(get_sulfide_value)
data['Disulfide bond']=values


# In[16]:


#Calculating length of sequences
length =pd.Series(data['Sequence'])
Sequence_Length=length.str.len()
position = pd.Series(len(data),dtype=np.str)
for i in range(len(data)):
    position[i]=''


# In[17]:


for i in range(len(data)):
    for j in range((Sequence_Length[i])):
        if (data['Sequence'][i][j])=='C':
            position[i]=position[i]+','+str(j)


# In[19]:


def get_integer_value(newdata):
    return re.findall('\d+',newdata)

data['Sequence positions']= position.apply(get_integer_value)
count =pd.Series(data['Sequence positions'])
length=count.str.len()
data['Sequence C count']=length


# In[20]:


data.to_excel('sequenceOutput.xlsx')


# In[21]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




