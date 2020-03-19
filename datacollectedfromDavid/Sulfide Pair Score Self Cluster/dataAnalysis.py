#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import re 
import matplotlib.pyplot as plt 
import sys
import time


import seaborn as sns
sns.set()


# In[3]:


df1 = pd.read_excel('data.xlsx')
len(df1)


# In[4]:


plt.figure(figsize=(10, 10))
plt.plot(df1['disulfie_pair_score_sum'],df1['Length'],'.')
plt.xlabel('Score')
plt.ylabel('Length')
plt.title('Score data Vs Length')
plt.savefig('allScore.jpg')


# In[5]:


frame1= df1[(df1['disulfie_pair_score_sum']==0)]
plt.figure(figsize=(15, 10))

plt.subplot(2,3,1)
plt.plot(frame1['disulfie_pair_score_sum'],frame1['Length'],'.')
plt.xlabel('Score')
plt.ylabel('Length')
plt.title('Score=0')
print('Total length with Score=0 are: ',len(frame1))

plt.subplot(2,3,2)
plt.plot(frame1['total glyco positions'],frame1['Length'],'.')
plt.title('Score=0')
plt.xlabel('Total Number of Glycosylation positions')
plt.ylabel('Length')

plt.subplot(2,3,3)
plt.plot(frame1['total sulphide bonds'],frame1['Length'],'.')
plt.title('Score=0')
plt.xlabel('Total Number of Sulfide positions')
plt.ylabel('Length')

plt.tight_layout()



plt.subplot(2,3,4)
plt.plot(frame1['disulfie_pair_score_sum'],frame1['total sulphide bonds'],'.')
plt.xlabel('Score')
plt.ylabel('total sulphide bonds')
plt.title('Score=0')

plt.subplot(2,3,5)
plt.plot(frame1['disulfie_pair_score_sum'],frame1['total glyco positions'],'.')
plt.title('Score=0')
plt.ylabel('Total Number of Glycosylation positions')
plt.xlabel('Score')

plt.subplot(2,3,6)
plt.plot(frame1['total sulphide bonds'],frame1['total glyco positions'],'.')
plt.title('Score=0')
plt.xlabel('total sulphide bonds')
plt.ylabel('total glyco positions')

plt.tight_layout()
plt.savefig('score0.jpg')


# In[6]:


frame2= df1[(df1['disulfie_pair_score_sum']>0) & (df1['disulfie_pair_score_sum']<=10)]
plt.figure(figsize=(15, 10))

plt.subplot(2,3,1)
plt.plot(frame2['disulfie_pair_score_sum'],frame2['Length'],'.')
plt.xlabel('Score')
plt.ylabel('Length')
plt.title('Score>0 and Score<=10')
print('Total Length with Score >0 and Score <= 10 are: ',len(frame2))

plt.subplot(2,3,2)
plt.plot(frame2['total glyco positions'],frame2['Length'],'.')
plt.title('Score>0 and Score<=10')
plt.xlabel('Total Number of Glycosylation positions')
plt.ylabel('Length')

plt.subplot(2,3,3)
plt.plot(frame2['total sulphide bonds'],frame2['Length'],'.')
plt.title('Score>0 and Score<=10')
plt.xlabel('Total Number of Sulfide positions')
plt.ylabel('Length')

plt.tight_layout()


plt.subplot(2,3,4)
plt.plot(frame2['disulfie_pair_score_sum'],frame2['total sulphide bonds'],'.')
plt.xlabel('Score')
plt.ylabel('total sulphide bonds')
plt.title('Score>0 and Score<=10')

plt.subplot(2,3,5)
plt.plot(frame2['disulfie_pair_score_sum'],frame2['total glyco positions'],'.')
plt.title('Score>0 and Score<=10')
plt.ylabel('Total Number of Glycosylation positions')
plt.xlabel('Score')

plt.subplot(2,3,6)
plt.plot(frame2['total sulphide bonds'],frame2['total glyco positions'],'.')
plt.title('Score>0 and Score<=10')
plt.xlabel('total sulphide bonds')
plt.ylabel('total glyco positions')

plt.tight_layout()
plt.savefig('score_0_10.jpg')


# In[7]:


frame3= df1[(df1['disulfie_pair_score_sum']>10) & (df1['disulfie_pair_score_sum']<=20)]
plt.figure(figsize=(15, 10))

plt.subplot(2,3,1)
plt.plot(frame3['disulfie_pair_score_sum'],frame3['Length'],'.')
plt.xlabel('Score')
plt.ylabel('Length')
plt.title('Score >10 and Score <= 20')
print('Total Length with Score >10 and Score <= 20 are: ',len(frame3))

plt.subplot(2,3,2)
plt.plot(frame3['total glyco positions'],frame3['Length'],'.')
plt.title('Score >10 and Score <= 20')
plt.xlabel('Total Number of Glycosylation positions')
plt.ylabel('Length')

plt.subplot(2,3,3)
plt.plot(frame3['total sulphide bonds'],frame3['Length'],'.')
plt.title('Score >10 and Score <= 20')
plt.xlabel('Total Number of Sulfide positions')
plt.ylabel('Length')

plt.tight_layout()

plt.subplot(2,3,4)
plt.plot(frame3['disulfie_pair_score_sum'],frame3['total sulphide bonds'],'.')
plt.xlabel('Score >10 and Score <= 20')
plt.ylabel('total sulphide bonds')
plt.title('Score=0')

plt.subplot(2,3,5)
plt.plot(frame3['disulfie_pair_score_sum'],frame3['total glyco positions'],'.')
plt.title('Score >10 and Score <= 20=0')
plt.ylabel('Total Number of Glycosylation positions')
plt.xlabel('Score')

plt.subplot(2,3,6)
plt.plot(frame3['total sulphide bonds'],frame3['total glyco positions'],'.')
plt.title('Score >10 and Score <= 20')
plt.xlabel('total sulphide bonds')
plt.ylabel('total glyco positions')

plt.tight_layout()
plt.savefig('score_10_20.jpg')


# In[8]:


frame4= df1[(df1['disulfie_pair_score_sum']>20)]
plt.figure(figsize=(15, 10))

plt.subplot(2,3,1)
plt.plot(frame4['disulfie_pair_score_sum'],frame4['Length'],'.')
plt.xlabel('Score')
plt.ylabel('Length')
plt.title('Score > 20')
print('Total Length with Score > 20 are: ',len(frame4))

plt.subplot(2,3,2)
plt.plot(frame4['total glyco positions'],frame4['Length'],'.')
plt.title('Score > 20')
plt.xlabel('Total Number of Glycosylation positions')
plt.ylabel('Length')

plt.subplot(2,3,3)
plt.plot(frame4['total sulphide bonds'],frame4['Length'],'.')
plt.title('Score > 20')
plt.xlabel('Total Number of Sulfide positions')
plt.ylabel('Length')

plt.tight_layout()


plt.subplot(2,3,4)
plt.plot(frame4['disulfie_pair_score_sum'],frame4['total sulphide bonds'],'.')
plt.xlabel('Score')
plt.ylabel('total sulphide bonds')
plt.title('Score=0')

plt.subplot(2,3,5)
plt.plot(frame4['disulfie_pair_score_sum'],frame4['total glyco positions'],'.')
plt.title('Score=0')
plt.ylabel('Total Number of Glycosylation positions')
plt.xlabel('Score')

plt.subplot(2,3,6)
plt.plot(frame4['total sulphide bonds'],frame4['total glyco positions'],'.')
plt.title('Score=0')
plt.xlabel('total sulphide bonds')
plt.ylabel('total glyco positions')

plt.tight_layout()
plt.savefig('score_20.jpg')


# In[12]:


frame1.to_excel('Cluster1(score=0).xlsx')
frame2.to_excel('Cluster2(score=0_10).xlsx')
frame3.to_excel('Cluster3(score_10_20).xlsx')
frame4.to_excel('Cluster4(score_20).xlsx')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




