#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd

# In[2]:


stats_path = pathlib.Path("../../data/cell_tracks_data").resolve(strict=True)
output_combined_stats_file_path = pathlib.Path(
    "../data/combined_stats.parquet"
).resolve()
output_combined_stats_file_path.parent.mkdir(parents=True, exist_ok=True)
# get the list of all files in the directory
files = sorted(stats_path.glob("*_stats.parquet"))
stats_df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
stats_df.to_parquet(output_combined_stats_file_path, index=False)
print(stats_df.shape)
stats_df.head()
