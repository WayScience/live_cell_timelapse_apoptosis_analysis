#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd

# In[2]:


input_dir = pathlib.Path("../results/mAP_cell_number_subsampled/").resolve(strict=True)
output_file = pathlib.Path("../results/mAP_cell_number_subsampled.parquet").resolve()
# get a list of all files in the directory
df = pd.concat([pd.read_parquet(file) for file in list(input_dir.glob("*"))])
df.to_parquet(output_file)
df.head()
