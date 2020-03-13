#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

# In[7]:


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


#Filename =sys.argv[1] # Read input file 
df = pd.read_excel('sequenceInput.xlsx')


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
human_data=data[data['Organism'].str.contains("Homo sapiens")] #Selecting only organism which contains 'Homo sapiens'
data= human_data.reset_index(drop=True)

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

data['intrabond']=data['intrabond'].astype(str)
#Calculate average distance of intrabond distance
average= data['intrabond'].apply(get_integer_value)
lenIntra=average.str.len()

averageintrabond= np.zeros(len(data))
for i in range(len(data)):
    totalsum=0
    for j in range(lenIntra[i]):
        totalsum= totalsum+int(average[i][j])
    averageintrabond[i]=totalsum/lenIntra[i]
data['average_Intrabond']=averageintrabond


# In[ ]:





# In[ ]:





# In[2]:


def get_float_value(data):
    return re.findall('\d+\.\d+',data)


# In[3]:


#Pair Score
array = pd.Series(data.size,dtype=np.str)
for x in range(data.size):
    array[x]=''

for i in range(0,len(data)):
    temp_array = np.zeros(lenB[i])
    for j in range(0,lenB[i]):
        for k in range(j+1,lenB[i]):
            if k!=lenB[i]:
                split1 = re.findall('\d+',data['Disulfide bond'][i][j]) 
                p1e1=int(split1[0])
                p1e2=int(split1[1])
                split2 = re.findall('\d+',data['Disulfide bond'][i][k]) 
                p2e1=int(split2[0])
                p2e2=int(split2[1])
            
                if(p2e1>p1e2):
                    temp_array[j]=temp_array[j]+0
                    temp_array[k]=temp_array[k]+0
                elif (p2e1<p1e2) and (p2e2>p1e2):
                    temp_array[j]=temp_array[j]+0.5
                    temp_array[k]=temp_array[k]+0.5    
                elif (p2e1<p1e2) and (p2e2<p1e2):
                    temp_array[j]=temp_array[j]+1.0
                    temp_array[k]=temp_array[k]+0.0
        array[i]=array[i]+str(temp_array[j])+','


# In[4]:


data['Disulfide_score']= array.apply(get_float_value)

#Claculating disulfide pair score average 
disulfie_pair_score_sum = np.zeros(len(data),dtype=np.float)
for i in range(len(data)):
    for j in range(len(data['Disulfide_score'][i])):
        disulfie_pair_score_sum[i]= disulfie_pair_score_sum[i]+ float(data['Disulfide_score'][i][j])
data['disulfie_pair_score_sum']=disulfie_pair_score_sum


# In[5]:


## Calculating Interbond Distance
interbond_distance= pd.Series(len(data),dtype= np.str)

for i in range(len(data)):
    interbond_distance[i]=''

for i in range(len(data)):
    
    if lenB[i]==0:
        interbond_distance[i]='No pair'
        
    if lenB[i]==1:
        interbond_distance[i]='Only One pair'
    
    for r in range(lenB[i]-1):
        split1= re.findall('\d+',data['Disulfide bond'][i][r])
        split2= re.findall('\d+',data['Disulfide bond'][i][r+1])
        interbond_distance[i]= interbond_distance[i]+ str(int(split2[0])-int(split1[1]))+ '|'
        
data['interbond_distance']=interbond_distance


# In[6]:




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


for i in range(len(data)):
    for j in range(len(data['Sequence'][i])-2):
        if data['Sequence'][i][j]=='N' and data['Sequence'][i][j+1]!='P' and (data['Sequence'][i][j+2]=='S' or data['Sequence'][i][j+2]=='T'):
            pos[i]=pos[i]+str(j)+','
data['N_X_S/T']=pos.apply(get_integer_value)
data['Sequence N_X_S/T count']=pd.Series(data['N_X_S/T']).str.len()           


# In[86]:


data.to_excel('sequenceOutput.xlsx')


# In[7]:





# In[ ]:





# In[ ]:




