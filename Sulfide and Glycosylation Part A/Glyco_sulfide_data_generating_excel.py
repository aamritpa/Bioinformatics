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

#1.1--------------------------------------------------------------------------------
#filename= sys.argv[1]
df = pd.read_excel('Input.xlsx')
data= pd.DataFrame()# Creating New Data Frame 
glyco_sulfide_frame= pd.DataFrame()   #Data Frame which contains Glycosilation and sulfide.

#1.2-----------------------------------------------------------------------------------
# Storing Usefull information From the database 'df' 
data['Entry']=df['Entry'] #Storing Entry
data['Entry name']=df['Entry name']  #Storing Entry Name
data['Protein names']=df['Protein names']  #Storing Protein Names
data['Gene names']=df['Gene names'] # Storing Gene Names
data['Length']=df['Length'] # Storing Length
data['Organism']=df['Organism'] # Storing Type of Organism
data['Mass']=df['Mass'] #Storing Mass
data['Status']=df['Status']
data['Disulfide bond']=df['Disulfide bond'] # Storing Disulfide bond 
data['Glycosylation']= df['Glycosylation']# Storing Glycosylation Position

#Now removing the extra data we dont need and containg only organism 'Human'
data=data[data['Organism'].str.contains("Homo sapiens")] #Selecting only organism which contains 'Homo sapiens'

data=data.dropna().reset_index(drop=True) # Drop all rows which contains Not a number and reset the index 

# Working on getting the relative positions of the Disulfide bond
temp_disulfide= data['Disulfide bond'] # Storing data temporary as 'disulfide_column' to resolve pointer problem.
disulfide_column= data['Disulfide bond']


# In[2]:



#1.3------------------------------------------------------------------------------------------
#Functions to use globally to get integer and float values.
def get_sulfide_value(newdata):
    return re.findall('\d+\..\d+',newdata)
def get_int_value(data):        #function defined for global purpose use to get integer values for any column.
    return re.findall('\d+',data)
def get_float_value(data):
    return re.findall('\d+\.\d+',data)
def get_Nlinked(data):
    return re.findall('\d+\;\s*\/note=\"N-linked',data)
#1.4-------------------------------------------------------------------------------------------
bond= temp_disulfide.apply(get_sulfide_value) # Function call which gives all the positons 
glyco_sulfide_frame['bond']=bond
data['Disulfide bond']=disulfide_column.apply(get_sulfide_value)# Storing extracted disulfide value in data Frame.

Glyco_data= data['Glycosylation'] # Storing data temporary as 'Glyco_data'
data['Glycosylation']= Glyco_data.apply(get_Nlinked) # the return data will be like 'CARBOHYD 110'
glyco_sulfide_frame['glyco']=data['Glycosylation']
data=data.dropna().reset_index(drop=True)
glyco_sulfide_frame=glyco_sulfide_frame.dropna().reset_index(drop=True)

#Now Removing the extra word 'CARBOHYD' and getting all positions of Glycosylation.
temp_data= data['Glycosylation'].astype(str)
data['Glycosylation']= temp_data.apply(get_int_value) 


# In[10]:



#1.5----------------------------------------------------------------------------------------------

# Creating a temporary dataframe for total length in columns and length in rows(like number of pairs in each index).

glyco_data =pd.Series(data['Glycosylation'])
sulfide_data =pd.Series(data['Disulfide bond'])
length_glyco=glyco_data.str.len() #length of Glycosylation
length_sulfide=sulfide_data.str.len() # length of Disulfide

rel_pos = pd.Series(glyco_data.size,dtype=np.str)
for i in range(glyco_data.size):
    rel_pos[i]=''
    
InsidePairs=pd.Series(sulfide_data.size,dtype=np.str)
InsidePairsLength=np.zeros(sulfide_data.size,dtype=np.int)
for i in range(glyco_data.size):
    InsidePairs[i]=''

# Now using loop and properties of Python library 're' to get values and finding positions between disulphide bond and glycosation.
i=0
while i<sulfide_data.size:
    for j in range(length_glyco[i]):
        Position=int(glyco_data[i][j])
        rel_pos[i]=rel_pos[i]+str(glyco_data[i][j])+'{|'
        for k in range(length_sulfide[i]):
            split = re.findall('\d+',sulfide_data[i][k]) 
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
#1.6-------------------------------------------------------------------------------------------
#Caluclating Score of each pairs of disulfide bonds
# O score if other pair outside.
# 0.5 score if other pair half inside and half outside.
# 1.0 score if other pair full inside.
array = pd.Series(glyco_data.size,dtype=np.str)
for x in range(glyco_data.size):
    array[x]=''

for i in range(0,len(data)):
    temp_array = np.zeros(length_sulfide[i])
    for j in range(0,length_sulfide[i]):
        for k in range(j+1,length_sulfide[i]):
            if k!=length_sulfide[i]:
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
data['Disulfide_score']= array.apply(get_float_value)

#Claculating disulfide pair score average 
disulfie_pair_score_sum = np.zeros(len(data),dtype=np.float)
for i in range(len(data)):
    for j in range(len(data['Disulfide_score'][i])):
        disulfie_pair_score_sum[i]= disulfie_pair_score_sum[i]+ float(data['Disulfide_score'][i][j])
data['disulfie_pair_score_sum']=disulfie_pair_score_sum

#1.7---------------------------------------------------------------------------------------
## Calculating Interbond Distance
interbond_distance= pd.Series(glyco_data.size,dtype= np.str)

for i in range(glyco_data.size):
    interbond_distance[i]=''

for i in range(glyco_data.size):
    
    if length_sulfide[i]==0:
        interbond_distance[i]='No pair'
        
    if length_sulfide[i]==1:
        interbond_distance[i]='Only One pair'
    
    for r in range(length_sulfide[i]-1):
        split1= re.findall('\d+',data['Disulfide bond'][i][r])
        split2= re.findall('\d+',data['Disulfide bond'][i][r+1])
        interbond_distance[i]= interbond_distance[i]+ str(int(split2[0])-int(split1[1]))+ '|'
        
data['interbond_distance']=interbond_distance
#1.8-------------------------------------------------------------------------------------------
#Intrabond distance formation
intrabond=pd.Series(sulfide_data.size,dtype=np.str)
for i in range(glyco_data.size):
    intrabond[i]=''

j=0;
while j<sulfide_data.size:
    for k in range(length_sulfide[j]):
        split = re.findall('\d+',sulfide_data[j][k]) 
        firstPair= int(split[0])
        secondPair=int(split[1])
        intrabond[j]=intrabond[j]+str(secondPair-firstPair)+'|'    
    j=j+1
data['intrabond']=intrabond
#1.9-----------------------------------------------------------------------------------------

#number of pairs of disuphide bonds
data['total sulphide bonds']=length_sulfide

#2.0-------------------------------------------------------------------------------------------
#number of glcosylation positions
data['total glyco positions']=length_glyco

#2.1-------------------------------------------------------------------------------------
#1st disulphide pair to N-terminus distance
first_pos = pd.Series(glyco_data.size,dtype=np.int)
for x in range(glyco_data.size):
    if(data['Disulfide bond'][x]!=[]):
        split = re.findall('\d+',data['Disulfide bond'][x][0]) 
        firstPair= int(split[0])
        first_pos[x]= firstPair
data['first_position_occurance']=first_pos

#last disulphide pair to C-terminus distance
last_pos= pd.Series(glyco_data.size,dtype=np.int)
for x in range(glyco_data.size):
    if(data['Disulfide bond'][x]!=[]):
        split= re.findall('\d+',data['Disulfide bond'][x][(length_sulfide[x]-1)])
        lastPair= int(split[1])
        last_pos[x] = lastPair 
data['last_position_occurance']=last_pos

#2.2--------------------------------------------------------------------------------
#Calculate average distance of intrabond distance
average= data['intrabond'].apply(get_int_value)
lenIntra=average.str.len()

averageintrabond= np.zeros(glyco_data.size)
for i in range(glyco_data.size):
    totalsum=0
    for j in range(lenIntra[i]):
        totalsum= totalsum+int(average[i][j])
    averageintrabond[i]=totalsum/lenIntra[i]
data['average_Intrabond']=averageintrabond

#2.3----------------------------------------------------------------------------------------
#calculating average distance of those pairs of disulphide bonds that have Glycosylation Inside.
temp_length=pd.Series(sulfide_data.size,dtype=np.str)
#for i in range(a.size):
#    temp_length[i]=''
inside_length = np.zeros(glyco_data.size,dtype=np.int)
for r in range(glyco_data.size):
    inside_length[r]=len(re.findall('\d+',data['In_Gly_Dis_Pair'][r]))/2
    temp_length[r]=re.findall('\d+',data['In_Gly_Dis_Pair'][r])
    if inside_length[r]!=0:
        InsidePairsLength[r]=InsidePairsLength[r]/inside_length[r]
data['avg_Inside_Pairs_Length']=InsidePairsLength

#2.4-------------------------------------------------------------------------------------
#Find those pairs which does not have glycosylation inside it.

#If the glyco is inside then putting the sulphide position 'nil' at that point.
i=0
while i<sulfide_data.size:
    for j in range(length_glyco[i]):
        Position=int(glyco_data[i][j])
        for k in range(length_sulfide[i]):
            if glyco_sulfide_frame['bond'][i][k]!='nil':
                split = re.findall('\d+',sulfide_data[i][k]) 
                firstPair= int(split[0])
                secondPair=int(split[1])
                if Position>= firstPair and Position<=secondPair:
                    glyco_sulfide_frame['bond'][i][k]='nil'
    i=i+1

#Calculating length between sulfide bonds.
temp=pd.Series(sulfide_data.size,dtype=np.str)
for i in range(glyco_data.size):
    temp[i]=''

for i in range(len(glyco_sulfide_frame['bond'])):
    for j in range(length_sulfide[i]):
        if glyco_sulfide_frame['bond'][i][j]!='nil':
            split= re.findall('\d+',glyco_sulfide_frame['bond'][i][j])
            temp[i]=str(temp[i] + str(int(split[1])-int(split[0]))+'|')
data['glyco_outside_bond']=glyco_sulfide_frame['bond']

#get the integer values.
intrabond_outside= temp.apply(get_int_value)
data['intrabond_glyco_outside_bond']=intrabond_outside
#2.5--------------------------------------------------------------------------------------


#calculating avegrage lenth of sulphide bonds which do not have glyco inside.

new =pd.Series(intrabond_outside)
newlen=new.str.len()
average= np.zeros(glyco_data.size)
for i in range(glyco_data.size):
    totalsum=0
    for j in range(newlen[i]):
        totalsum= totalsum+int(intrabond_outside[i][j])
    if newlen[i]!=0:
        average[i]=totalsum/newlen[i]
data['average_Intrabond_outside_glyco_bond']=average

#2.6-------------------------------------------------------------------------------------------
data.to_excel('Glyco_sulfide_data_generated.xlsx')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




