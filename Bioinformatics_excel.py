#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

# In[38]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys


#Input file must be given to run the code.
#For Example Run(in terminal) 'python3 Bioinformatics input.xlsx'

Filename =sys.argv[1] # Read input file 

df = pd.read_excel(Filename)# Read EXCEL Files and storing in database 'df'
data= pd.DataFrame()# Creating New Data Frame 


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
columns=data[data['Organism'].str.contains("Human")] #Selecting only organism which contains HUMANS
data= columns

data['Disulfide bond']=df['Disulfide bond'] # Storing Disulfide bond 
data['Glycosylation']= df['Glycosylation']# Storing Glycosylation Position
data=data.dropna().reset_index(drop=True) # Drop all rows which contains Not a number and reset the index 

# Working on getting the relative positions of the Disulfide bond
disulfide_column= data['Disulfide bond'] # Storing data temporary as 'disulfide_column'
def get_sulfide_value(newdata):
    return re.findall('\d+ \d+',newdata)

data['Disulfide bond']= disulfide_column.apply(get_sulfide_value) # Function call which gives all the positons 


# In[2]:



#Working on getting the positions of N-linked Glycosylation
Glyco_data= data['Glycosylation'] # Storing data temporary as 'Glyco_data'

def getGlycoNLinked(data):   #Making a Function to get the positions
    return re.findall('\d+\s* N-linked',data) # The return data will contain 'position and N-linked pattern'
data['Glycosylation']= Glyco_data.apply(getGlycoNLinked) # the return data will be like '49 N-linked'


#Now Removing the extra word 'N-linked' and getting all positions 
temp_data= data['Glycosylation'].astype(str)
def get_gly_value(data):
    return re.findall('\d+',data) # return all positions
data['Glycosylation']= temp_data.apply(get_gly_value) 


# Creating a temporary dataframe for total length in columns and length in rows(like number of pairs in each index).
a= pd.DataFrame()
b= pd.DataFrame()
bond=pd.DataFrame()
a =pd.Series(data['Glycosylation'])
b =pd.Series(data['Disulfide bond'])
bond=pd.DataFrame()
bond['data']=data['Disulfide bond']
lenB=b.str.len()
lenA=a.str.len()

rel_pos = pd.Series(a.size,dtype=np.str)
for i in range(a.size):
    rel_pos[i]=''
    
InsidePairs=pd.Series(b.size,dtype=np.str)
InsidePairsLength=np.zeros(b.size,dtype=np.int)
for i in range(a.size):
    InsidePairs[i]=''

# Now using loop and properties of Python library 're' to get values and finding positions between disulphide bond and glycosation.

i=0
while i<b.size:
    for j in range(lenA[i]):
        Position=int(a[i][j])
        rel_pos[i]=rel_pos[i]+str(a[i][j])+'{|'
        for k in range(lenB[i]):
            split = re.findall('\d+',b[i][k]) 
            firstPair= int(split[0])
            secondPair=int(split[1])
            if Position< firstPair:
                rel_pos[i] = rel_pos[i]+'o'+str(Position-firstPair)
                rel_pos[i] = rel_pos[i]+'o'+str(Position-secondPair)
            elif Position>= firstPair and Position<=secondPair:
                rel_pos[i]=  rel_pos[i]+'i' + str(Position-firstPair)
                rel_pos[i] = rel_pos[i]+'i'+ str(Position-secondPair)
                InsidePairs[i]= InsidePairs[i] + str('('+ str(firstPair)+','+str(secondPair)+')')
                InsidePairsLength[i]=InsidePairsLength[i]+ (secondPair-firstPair)
            else:
                rel_pos[i] = rel_pos[i]+'o'+str(Position-firstPair)
                rel_pos[i] = rel_pos[i]+'o'+str(Position-secondPair)
            rel_pos[i]=rel_pos[i]+'|'    
        rel_pos[i]=rel_pos[i]+'}'
    i=i+1
data['rel_pos']= rel_pos
data['In_Gly_Dis_Pair']=InsidePairs


# In[3]:


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
        
data['interbond_distance']=interbond_distance

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


#number of pairs of disuphide bonds
data['No. Disulphide bonds']=lenB


#1st disulphide pair to N-terminus distance
first_pos = pd.Series(a.size,dtype=np.int)
for r in range(a.size): 
    split = re.findall('\d+',data['Disulfide bond'][r][0]) 
    firstPair= int(split[0])
    first_pos[r]= firstPair
data['first_position_occurance']=first_pos


#last disulphide pair to C-terminus distance
last_pos= pd.Series(a.size,dtype=np.int)
for r in range(a.size):
    split= re.findall('\d+',data['Disulfide bond'][r][(lenB[r]-1)])
    lastPair= int(split[1])
    last_pos[r] = lastPair 
data['last_position_occurance']=last_pos

#Calculate average distance of intrabond distance
def intrabond_Average(data):
    return re.findall('\d+',data)
average= data['intrabond'].apply(intrabond_Average)
lenIntra=average.str.len()

averageintrabond= np.zeros(a.size)
for i in range(a.size):
    totalsum=0
    for j in range(lenIntra[i]):
        totalsum= totalsum+int(average[i][j])
    averageintrabond[i]=totalsum/lenIntra[i]
data['average_Intrabond']=averageintrabond


# In[4]:


#calculating average distance of those pairs of disulphide bonds that have Glycosylation Inside.
temp_length=pd.Series(b.size,dtype=np.str)
for i in range(a.size):
    temp_length[i]=''
inside_length = np.zeros(a.size,dtype=np.int)
for r in range(a.size):
    inside_length[r]=len(re.findall('\d+',data['In_Gly_Dis_Pair'][r]))/2
    temp_length[r]=re.findall('\d+',data['In_Gly_Dis_Pair'][r])
    if inside_length[r]!=0:
        InsidePairsLength[r]=InsidePairsLength[r]/inside_length[r]
data['avg_Inside_Pairs_Length']=InsidePairsLength


# In[5]:


#If the glyco is inside then putting the sulphide position 'nil' at that point.
i=0
while i<b.size:
    for j in range(lenA[i]):
        Position=int(a[i][j])
        for k in range(lenB[i]):
            if b[i][k]!='nil':
                split = re.findall('\d+',b[i][k]) 
                firstPair= int(split[0])
                secondPair=int(split[1])
                if Position>= firstPair and Position<=secondPair:
                    bond['data'][i][k]='nil'
    i=i+1


# In[6]:


#calculating length between sulfide bonds.
temp=pd.Series(b.size,dtype=np.str)
for i in range(a.size):
    temp[i]=''

for i in range(len(bond)):
    for j in range(lenB[i]):
        if bond['data'][i][j]!='nil':
            split= re.findall('\d+',bond['data'][i][j])
            temp[i]=str(temp[i] + str(int(split[1])-int(split[0]))+'|')
data['glyco_outside_bond']=bond['data']


# In[9]:


#get the integer values.
def get_gly_outside_values(temp):
    return re.findall('\d+',temp)
intrabond_outside= temp.apply(get_gly_outside_values)
data['intrabond_glyco_outside_bond']=intrabond_outside


# In[10]:


#calculating avegrage lenth of sulphide bonds which do not have glyco inside.
new =pd.Series(intrabond_outside)
newlen=new.str.len()
average= np.zeros(a.size)
for i in range(a.size):
    totalsum=0
    for j in range(newlen[i]):
        totalsum= totalsum+int(intrabond_outside[i][j])
    if newlen[i]!=0:
        average[i]=totalsum/newlen[i]
data['average_Intrabond_outside_glyco_bond']=average




# In[12]:


data.to_excel('output.xlsx')


# In[ ]:




