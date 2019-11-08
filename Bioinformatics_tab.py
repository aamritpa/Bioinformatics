#!/usr/bin/env python
# coding: utf-8

# In[128]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys


# In[129]:


df = pd.read_csv("uniprot-human-filtered-reviewed%3Ayes.tab", sep='\t')


# In[132]:


# Storing Usefull information From the database 'df' 
data=pd.DataFrame()
data['Entry']=df['Entry'] #Storing Entry
data['Entry name']=df['Entry name']  #Storing Entry Name
data['Protein names']=df['Protein names']  #Storing Protein Names
data['Gene names']=df['Gene names'] # Storing Gene Names
data['Length']=df['Length'] # Storing Length


# In[133]:


data['Organism']=df['Organism'] # Storing Type of Organism
data['Mass']=df['Mass'] #Storing Mass
 


# In[134]:


#Now removing the extra data we dont need and containg only organism 'Human'
columns=data[data['Organism'].str.contains("Human")] #Selecting only organism which contains HUMANS
data= columns


# In[135]:


data['Disulfide bond']= df['Disulfide bond'] # Storing Disulfide bond 
data['Glycosylation']= df['Glycosylation']# Storing Glycosylation Position
data


# In[136]:


data=data.dropna().reset_index(drop=True) # Drop all rows which contains Not a number and reset the index 
data


# In[137]:


# Working on getting the relative positions of the Disulfide bond
disulfide_column= data['Disulfide bond'] # Storing data temporary as 'disulfide_column'
def getdata(newdata):
    return re.findall('\d+ \d+',newdata)

data['Disulfide bond']= disulfide_column.apply(getdata) # Function call which gives all the positons 


# In[138]:


#Working on getting the positions of N-linked Glycosylation
Glyco_data= data['Glycosylation'] # Storing data temporary as 'Glyco_data'

def getGlycoNLinked(data):   #Making a Function to get the positions
    return re.findall('\d+\s* N-linked',data) # The return data will contain 'position and N-linked pattern'
data['Glycosylation']= Glyco_data.apply(getGlycoNLinked) # the return data will be like '49 N-linked'


# In[139]:


#Now Removing the extra word 'N-linked' and getting all positions 
temp_data= data['Glycosylation'].astype(str)
def getValue(data):
    return re.findall('\d+',data) # return all positions
data['Glycosylation']= temp_data.apply(getValue) 


# In[ ]:





# In[141]:


a= pd.DataFrame()
b= pd.DataFrame()

a =pd.Series(data['Glycosylation'])
b =pd.Series(data['Disulfide bond'])

lenB=b.str.len()
lenA=a.str.len()
rel_pos = pd.Series(a.size,dtype=np.str)

for i in range(a.size):
    rel_pos[i]=''


# In[142]:


i=0;
while i<b.size:
    for j in range(lenA[i]):
        Position=int(a[i][j])
        rel_pos[i]=rel_pos[i]+'{'
        for k in range(lenB[i]):
            split = re.findall('\d+',b[i][k]) 
            firstPair= int(split[0])
            secondPair=int(split[1])
            if Position< firstPair:
                rel_pos[i] = rel_pos[i]+'o'+str(Position-firstPair)
            else:
                rel_pos[i]=  rel_pos[i]+'i' + str(Position-firstPair)
            if Position<=secondPair:
                rel_pos[i] = rel_pos[i]+'i'+ str(secondPair-Position)
            else:
                rel_pos[i] = rel_pos[i]+'o'+str(Position-secondPair)
        rel_pos[i]=rel_pos[i]+'}'
    i=i+1


# In[143]:


data['rel_pos']= rel_pos


# In[50]:


#Calculating Interbond Distance
interbond_distance= pd.Series(a.size,dtype= np.str)

for i in range(a.size):
    interbond_distance[i]=''

for i in range(a.size):
    
    if lenB[i]==0:
        interbond_distance[i]='No pair'
        
    if lenB[i]==1:
        interbond_distance[i]='Only One pair'
    
    for r in range(lenB[i]-1):
        split1= re.findall('\d+',data['Disulfide bond'][i][r])
        split2= re.findall('\d+',data['Disulfide bond'][i][r+1])
        interbond_distance[i]= interbond_distance[i]+ str(int(split2[0])-int(split1[1]))+ '|'
        
data['interbond_distance']=distance  


# In[51]:


#Intrabond distance formation
intrabond=pd.Series(b.size,dtype=np.str)
for i in range(a.size):
    intrabond[i]=''
j=0;
while j<b.size:
    for k in range(lenB[j]):
        split = re.findall('\d+',b[j][k]) 
        firstPair= int(split[0])
        secondPair=int(split[1])
        intrabond[j]=intrabond[j]+str(secondPair-firstPair)+'|'    
    j=j+1
data['intrabond']=intrabond


# In[52]:


#number of pairs of disuphide bonds
data['No. Disulphide bonds']=lenB


# In[53]:


#1st disulphide pair to N-terminus distance
first_pos = pd.Series(a.size,dtype=np.int)
for r in range(a.size): 
    split = re.findall('\d+',data['Disulfide bond'][r][0]) 
    firstPair= int(split[0])
    first_pos[r]= firstPair
data['first_position_occurance']=first_pos


# In[54]:


#last disulphide pair to C-terminus distance
last_pos= pd.Series(a.size,dtype=np.int)
for r in range(a.size):
    split= re.findall('\d+',data['Disulfide bond'][r][(lenB[r]-1)])
    lastPair= int(split[1])
    last_pos[r] = lastPair 
data['last_position_occurance']=last_pos


# In[55]:


data.to_excel('output.xlsx') # Output an Excel File with refined data.


# In[36]:





