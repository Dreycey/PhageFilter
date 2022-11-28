#!/usr/bin/env python
# coding: utf-8

# # Bench Plotting
# Given the output from the benchmarking, this notebook plots relevant figures.
# 

# In[90]:


import pandas
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set(style="whitegrid", palette="pastel", rc={"figure.dpi":300, 'savefig.dpi':300})


# ## Helper Methods

# In[84]:


def csv2dataframe(csv_file):
    """
    Converts the csv file to a pandas dataframe.
    """
    return pandas.read_csv(csv_file)


# ## Genomes Benchmarking(Plotting)

# In[342]:


genomes_csv = "results/res_genomes.csv"


# In[343]:


# open dataframe
df = csv2dataframe(genomes_csv)
df["Memory (GB)"] = df[" memory (bytes)"] / (10)**6
df["Time (minutes)"] = df[" time (ns)"] / (10**9 * 60)


# ### plotting

# In[344]:


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
# add data
ax1.scatter(df["genome count"], df["Memory (GB)"], s=100, color='#1357a6', label="memory")
ax2.scatter(df["genome count"], df["Time (minutes)"], color='red', edgecolor="black", label="time")
# configure axis labels
plt.title('Performance building gSBT (non-compiled; using JSON) \n')
ax1.set_xlabel('Number of genomes added to gSBT')
ax1.set_ylabel('Memory to build gSBT (MB)', color='#1357a6')
ax2.set_ylabel('Time to build gSBT (minutes)', color='red')
plt.show()


# ## Parameterization Benchmarking (Plotting)

# In[348]:


pameterization_csv = "results/res_parameterization.csv"


# In[357]:


# open dataframe
df = csv2dataframe(pameterization_csv)
df["F1"] = 2* df[" precision"]*df[" recall"] / (df[" precision"]+df[" recall"])
df[" read count error"] = pandas.to_numeric(df[" avg read count error"].replace(" nan", np.nan))
#df =  df.groupby(["kmer size", " theta"], as_index=False).mean()

df2 =  df.groupby(["kmer size", " number of genomes"], as_index=False).mean()

df3 =  df.groupby(["kmer size", " error rate"], as_index=False).mean()

df4 =  df.groupby([" theta", " error rate"], as_index=False).mean()

df =  df.groupby(["kmer size", " theta"], as_index=False).mean()


# In[358]:


sns.heatmap(data=df2.pivot(" number of genomes", "kmer size", "F1"), annot=True, linewidth=.5, cmap=sns.color_palette("ch:start=.2,rot=-.3"))


# In[359]:


ax = sns.violinplot(data=df2, x=" number of genomes", y="F1", hue=" number of genomes")


# In[360]:


ax = sns.scatterplot(data=df3, x="kmer size", y="F1", hue=" error rate", s=100)


# In[361]:


ax = sns.scatterplot(data=df4, x=" theta", y="F1", hue=" error rate", s=50)


# # ax = sns.heatmap(data=df.pivot("kmer size", " theta", " recall"), annot=True, linewidth=.5, cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
# , vmin=0, vmax=1.0)
# ax.set_title("Parameterization \n")

# In[364]:


ax = sns.heatmap(data=df.pivot("kmer size", " theta", "F1"), annot=True, linewidth=.5, cmap=sns.color_palette("ch:start=.2,rot=-.3"))
ax.set_title("F1 \n")


# In[365]:


ax = sns.heatmap(data=df.pivot("kmer size", " theta", " read count error"), annot=True, linewidth=.5, cmap=sns.color_palette("ch:start=.2,rot=-.3"))
ax.set_title("Read Count Error \n")


# ## Relative Performance Benchmarking (Plotting)

# In[315]:


relative_performance_csv = "results/res_relative_performance.csv"


# In[316]:


# open dataframe
df = csv2dataframe(relative_performance_csv)
df["F1"] = 2* df[" precision"]*df[" recall"] / (df[" precision"]+df[" recall"])
df["Mean read count error"] = pandas.to_numeric(df[" read count error"].replace(" nan", np.nan))
df["Memory (GB)"] = df[" memory"] / (10)**6
df["Time (minutes)"] = df[" time"] / (10**9 * 60)


# In[317]:


sns.violinplot(data=df, x="tool name", y=" recall", cut=0, inner="point", palette=sns.color_palette("ch:start=.2,rot=-.3"))


# In[318]:


sns.violinplot(data=df, x="tool name", y=" precision", cut=0, inner="point", palette=sns.color_palette("ch:start=.2,rot=-.3"))


# In[319]:


sns.violinplot(data=df, x="tool name", y="F1", cut=0, inner="point", palette=sns.color_palette("ch:start=.2,rot=-.3"))


# In[282]:


sns.violinplot(data=df, x="tool name", y="Mean read count error", cut=0, palette=sns.color_palette("ch:start=.2,rot=-.3"))


# In[283]:


sns.violinplot(data=df, x="tool name", y="Time (minutes)", cut=0, palette=sns.color_palette("ch:start=.2,rot=-.3"))


# In[284]:


sns.violinplot(data=df, x="tool name", y="Memory (GB)", cut=0, palette=sns.color_palette("ch:start=.2,rot=-.3"))


# In[ ]:





# In[ ]:




