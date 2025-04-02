#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
from copairs import map
from copairs.matching import assign_reference_index

# In[2]:


endpoint_path = pathlib.Path(
    "../../data/CP_feature_select/endpoint_whole_image/feature_selected_whole_image.parquet"
).resolve(strict=True)

map_results_path = pathlib.Path("../data/0.ground_truth/map.parquet").resolve()
map_results_path.parent.mkdir(parents=True, exist_ok=True)

endpoint_df = pd.read_parquet(endpoint_path)
endpoint_df.head()


# In[3]:


reference_col = "Metadata_reference_index"
df_activity = assign_reference_index(
    endpoint_df,
    "Metadata_dose == 0.0",
    reference_col=reference_col,
    default_value=-1,
)
df_activity.head()


# In[4]:


pos_sameby = ["Metadata_dose", reference_col]
pos_diffby = []
neg_sameby = []
neg_diffby = ["Metadata_dose", reference_col]
metadata = df_activity.filter(regex="Metadata")
profiles = df_activity.filter(regex="^(?!Metadata)").values

activity_ap = map.average_precision(
    metadata, profiles, pos_sameby, pos_diffby, neg_sameby, neg_diffby
)

activity_ap = activity_ap.query("Metadata_dose != 0.0")
activity_ap.head()


# In[5]:


activity_map = map.mean_average_precision(
    activity_ap, pos_sameby, null_size=1000000, threshold=0.05, seed=0
)
activity_map.reset_index(inplace=True)
activity_map.to_parquet(map_results_path)
activity_map.head()
