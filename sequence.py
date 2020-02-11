#!/usr/bin/env python
# coding: utf-8

# In[76]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys
import time


# In[77]:


Filename =sys.argv[1] # Read input file 
df = pd.read_excel(Filename)


# In[78]:


df=df.dropna().reset_index(drop=True)
data= pd.DataFrame()
data=df

def get_integer_value(newdata):
    return re.findall('\d+',newdata)

# In[79]:


sulfide=df['Disulfide bond']
def get_sulfide_value(newdata):
    return re.findall('\d+\..\d+',newdata)

values= sulfide.apply(get_sulfide_value)
data['Disulfide bond']=values


# In[80]:


#Intrabond distance formation

lenB=pd.Series(data['Disulfide bond']).str.len()
intrabond=pd.Series(len(data),dtype=np.str)
for i in range(len(data)):
    intrabond[i]=''
j=0;
while j<len(data):
    for k in range(lenB[j]):
        split = re.findall('\d+',data['Disulfide bond'][j][k]) 
        firstPair= int(split[0])
        secondPair=int(split[1])
        intrabond[j]=intrabond[j]+str(secondPair-firstPair)+','    
    j=j+1
data['intrabond']=intrabond.apply(get_integer_value)


# In[81]:


#Calculating length of sequences
Sequence_Length=pd.Series(data['Sequence']).str.len()
position = pd.Series(len(data),dtype=np.str)
for i in range(len(data)):
    position[i]=''


# In[82]:


for i in range(len(data)):
    for j in range((Sequence_Length[i])):
        if (data['Sequence'][i][j])=='C':
            position[i]=position[i]+','+str(j)


# In[83]:




data['Sequence C position']= position.apply(get_integer_value)
count =pd.Series(data['Sequence C position'])
length=count.str.len()
data['Sequence C count']=length


# In[84]:


pos = pd.Series(len(data),dtype=np.str)
for i in range(len(data)):
    pos[i]=''


# In[ ]:





# In[85]:


for i in range(len(data)):
    for j in range(len(data['Sequence'][i])-2):
        if data['Sequence'][i][j]=='N' and data['Sequence'][i][j+1]!='P' and (data['Sequence'][i][j+2]=='S' or data['Sequence'][i][j+2]=='T'):
            pos[i]=pos[i]+str(j)+','
data['N_X_S/T']=pos.apply(get_integer_value)
data['Sequence N_X_S/T count']=pd.Series(data['N_X_S/T']).str.len()           


# In[86]:


data.to_excel('sequenceOutput.xlsx')


# In[87]:


data


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




