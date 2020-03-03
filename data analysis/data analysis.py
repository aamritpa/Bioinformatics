#!/usr/bin/env python
# coding: utf-8

# In[4]:


from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 


# In[5]:


df = pd.read_excel('output.xlsx')


# In[6]:


case1=df['average_Intrabond']
temp=df[df['avg_Inside_Pairs_Length']!=0].reset_index(drop=True)
case2=temp['avg_Inside_Pairs_Length']
c=df['average_Intrabond_outside_glyco_bond']
case3=c.dropna().reset_index(drop=True)


# In[28]:


plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.hist(case1,bins=60)
plt.xlabel('Number of Entries')
plt.ylabel('Average length')
plt.subplot(1, 2, 2)
plt.plot(case1)
plt.xlabel('Number of Entries')
plt.ylabel('Average length')

plt.savefig('plt1.png')


# In[29]:


plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.hist(case2,bins=60)
plt.xlabel('Number of Entries')
plt.ylabel('Average length')
plt.subplot(1, 2, 2)
plt.plot(case2)
plt.xlabel('Number of Entries')
plt.ylabel('Average length')
plt.savefig('plt2.png')


# In[30]:


plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.hist(case3,bins=60)
plt.xlabel('Number of Entries')
plt.ylabel('Average length')
plt.subplot(1, 2, 2)
plt.plot(case3)
plt.xlabel('Number of Entries')
plt.ylabel('Average length')
plt.savefig('plt3.png')


# In[96]:


#Statistical Tests
ttest1 = stats.ttest_ind(case1,case2)
ttest2 = stats.ttest_ind(case1,case3)
ttest3 = stats.ttest_ind(case2,case3)
print('ttest1 between Case1 and Case2', ttest1)
print('ttest1 between Case1 and Case3', ttest2)
print('ttest1 between Case2 and Case3', ttest3)
print('\n')

#CHECKING MEAN
print('The mean for average_Intrabond is',case1.mean(axis=0))
print('The mean for avg_Inside_Pairs_Length is',case2.mean(axis=0))
print('The mean for intrabond_glyco_outside_bond is',case3.mean(axis=0))
print('\n')

#EQUAL VARIANCE TEST
print('Equal Variance between Case1 and Case2',stats.levene(case1,case2).pvalue)
print('Equal Variance between Case1 and Case3',stats.levene(case2,case3).pvalue)
print('Equal Variance between Case2 and Case3',stats.levene(case3,case1).pvalue)


# In[31]:


plt.scatter(df['average_Intrabond'],df['avg_Inside_Pairs_Length'],c='blue')
plt.savefig('plt3.png')
plt.xlabel('Case1')
plt.ylabel('Case2')
plt.savefig('case1case2.png')
#Correlation Coeeficiant
print(stats.linregress(df['average_Intrabond'],df['avg_Inside_Pairs_Length']).rvalue)


# In[32]:


plt.scatter(df['average_Intrabond'],df['average_Intrabond_outside_glyco_bond'],c='blue')
plt.savefig('plt3.png')
plt.xlabel('Case1')
plt.ylabel('Case3')
plt.savefig('case1case3.png')
#Correlation Coeeficiant
print(stats.linregress(df['average_Intrabond'],df['average_Intrabond_outside_glyco_bond']).rvalue)


# In[33]:


plt.scatter(df['avg_Inside_Pairs_Length'],df['average_Intrabond_outside_glyco_bond'],c='blue')
plt.savefig('plt3.png')
plt.xlabel('Case2')
plt.ylabel('Case3')
plt.savefig('case3case2.png')
#Correlation Coeeficiant
print(stats.linregress(df['average_Intrabond_outside_glyco_bond'],df['avg_Inside_Pairs_Length']).rvalue)


# In[ ]:




