import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys
import time



#Input file must be given to run the code.
#For Example Run(in terminal) 'python3 Bioinformatics input.xlsx'

filename= sys.argv[1]
df = pd.read_csv(filename, sep='\t')


data= pd.DataFrame()# Creating New Data Frame 

glyco_sulfide_frame= pd.DataFrame()   #Data Frame which contains Glycosilation and sulfide.

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


# In[ ]:





# In[31]:


def get_int_value(data):        #function defined for global purpose use to get integer values for any column.
    return re.findall('\d+',data)
def get_float_value(data):
    return re.findall('\d+\.\d+',data)


# In[ ]:





# In[32]:


N_linked = re.compile(r'N-linked')
def get_Nlinked(data):
    match = N_linked.search(data)
    if match:
        return re.findall('CARBOHYD \d+',data)
    else:
        return None
    
Glyco_data= data['Glycosylation'] # Storing data temporary as 'Glyco_data'
data['Glycosylation']= Glyco_data.apply(get_Nlinked) # the return data will be like 'CARBOHYD 110'
glyco_sulfide_frame['glyco']=data['Glycosylation']
data=data.dropna().reset_index(drop=True)
glyco_sulfide_frame=glyco_sulfide_frame.dropna().reset_index(drop=True)

#Now Removing the extra word 'CARBOHYD' and getting all positions of Glycosylation.
temp_data= data['Glycosylation'].astype(str)
data['Glycosylation']= temp_data.apply(get_int_value) 


# In[ ]:





# In[33]:


# Creating a temporary dataframe for total length in columns and length in rows(like number of pairs in each index).
a= pd.DataFrame()
b= pd.DataFrame()
a =pd.Series(data['Glycosylation'])
b =pd.Series(data['Disulfide bond'])
lenA=a.str.len() #length of Glycosylation
lenB=b.str.len() # length of Disulfide


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


# In[34]:


lenB[110]


# In[35]:


#Caluclating Score of each pairs of disulfide bonds
# O score if other pair outside.
# 0.5 score if other pair half inside and half outside.
# 1.0 score if other pair full inside.
array = pd.Series(a.size,dtype=np.str)
for x in range(a.size):
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
data['Disulfide_score']= array.apply(get_float_value)


# In[36]:


#Claculating disulfide pair score average 
disulfie_pair_score_sum = np.zeros(len(data),dtype=np.float)
for i in range(len(data)):
    for j in range(len(data['Disulfide_score'][i])):
        disulfie_pair_score_sum[i]= disulfie_pair_score_sum[i]+ float(data['Disulfide_score'][i][j])
data['disulfie_pair_score_sum']=disulfie_pair_score_sum


# In[ ]:





# In[37]:


## Calculating Interbond Distance
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
data['total sulphide bonds']=lenB

#number of glcosylation positions
data['total glyco positions']=lenA


#1st disulphide pair to N-terminus distance
first_pos = pd.Series(a.size,dtype=np.int)
for x in range(a.size):
    if(data['Disulfide bond'][x]!=[]):
        split = re.findall('\d+',data['Disulfide bond'][x][0]) 
        firstPair= int(split[0])
        first_pos[x]= firstPair
data['first_position_occurance']=first_pos


#last disulphide pair to C-terminus distance
last_pos= pd.Series(a.size,dtype=np.int)
for x in range(a.size):
    if(data['Disulfide bond'][x]!=[]):
        split= re.findall('\d+',data['Disulfide bond'][x][(lenB[x]-1)])
        lastPair= int(split[1])
        last_pos[x] = lastPair 
data['last_position_occurance']=last_pos

#Calculate average distance of intrabond distance
average= data['intrabond'].apply(get_int_value)
lenIntra=average.str.len()

averageintrabond= np.zeros(a.size)
for i in range(a.size):
    totalsum=0
    for j in range(lenIntra[i]):
        totalsum= totalsum+int(average[i][j])
    averageintrabond[i]=totalsum/lenIntra[i]
data['average_Intrabond']=averageintrabond


# In[38]:


#calculating average distance of those pairs of disulphide bonds that have Glycosylation Inside.
temp_length=pd.Series(b.size,dtype=np.str)
#for i in range(a.size):
#    temp_length[i]=''
inside_length = np.zeros(a.size,dtype=np.int)
for r in range(a.size):
    inside_length[r]=len(re.findall('\d+',data['In_Gly_Dis_Pair'][r]))/2
    temp_length[r]=re.findall('\d+',data['In_Gly_Dis_Pair'][r])
    if inside_length[r]!=0:
        InsidePairsLength[r]=InsidePairsLength[r]/inside_length[r]
data['avg_Inside_Pairs_Length']=InsidePairsLength


# In[ ]:





# In[ ]:





# In[ ]:





# In[39]:


#If the glyco is inside then putting the sulphide position 'nil' at that point.
i=0
while i<b.size:
    for j in range(lenA[i]):
        Position=int(a[i][j])
        for k in range(lenB[i]):
            if glyco_sulfide_frame['bond'][i][k]!='nil':
                split = re.findall('\d+',b[i][k]) 
                firstPair= int(split[0])
                secondPair=int(split[1])
                if Position>= firstPair and Position<=secondPair:
                    glyco_sulfide_frame['bond'][i][k]='nil'
    i=i+1


# In[40]:


#Calculating length between sulfide bonds.
temp=pd.Series(b.size,dtype=np.str)
for i in range(a.size):
    temp[i]=''

for i in range(len(glyco_sulfide_frame['bond'])):
    for j in range(lenB[i]):
        if glyco_sulfide_frame['bond'][i][j]!='nil':
            split= re.findall('\d+',glyco_sulfide_frame['bond'][i][j])
            temp[i]=str(temp[i] + str(int(split[1])-int(split[0]))+'|')
data['glyco_outside_bond']=glyco_sulfide_frame['bond']


# In[41]:


#get the integer values.
intrabond_outside= temp.apply(get_int_value)
data['intrabond_glyco_outside_bond']=intrabond_outside


# In[42]:


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


# In[ ]:





# In[43]:


data.to_excel('finaloutput.xlsx')
