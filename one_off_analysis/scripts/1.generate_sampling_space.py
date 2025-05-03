#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import random

import numpy as np
import pandas as pd

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# In[2]:


number_of_cells_path = pathlib.Path("../combinations/number_of_cells.txt").resolve()

number_of_cells_path.parent.mkdir(parents=True, exist_ok=True)


# In[3]:


data_file_path = pathlib.Path(
    "../../data/CP_feature_select/profiles/features_selected_profile.parquet"
).resolve(strict=True)
df = pd.read_parquet(data_file_path)

df.head()


# In[ ]:


# In[4]:


random.seed(0)
num_of_iterations = 100
seed_space = np.random.randint(0, 1000000, num_of_iterations).tolist()  # min, max
min_number_of_cells = df["Metadata_number_of_singlecells"].min()
cells_to_sample_per_well = np.linspace(
    1, min_number_of_cells, min_number_of_cells
).tolist()  # start, stop, num
cells_to_sample_per_well = [int(x) for x in cells_to_sample_per_well]


with open(number_of_cells_path, "w") as f:
    for cells_to_sample in cells_to_sample_per_well:
        f.write(f"{cells_to_sample}\n")
