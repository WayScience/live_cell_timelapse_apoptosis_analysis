#!/usr/bin/env python
# coding: utf-8

# This notebook generates the umap embeddings of the images in the dataset. The embeddings are saved in a parquet file.

# In[1]:


import argparse
import pathlib

import numpy as np
import pandas as pd
import umap

# In[ ]:


# set the arg parser
parser = argparse.ArgumentParser(description="UMAP on a matrix")

parser.add_argument("--data_mode", type=str, default="CP", help="data mode to use")

# get the args
args = parser.parse_args()

# set data mode to either "CP" or "scDINO" or "combined"
data_mode = args.data_mode


# In[2]:


# set the paths to the data
CP_fs_sc_profiles_path = pathlib.Path(
    "../../data/feature_selected_data/run_20230920ChromaLiveTL_24hr4ch_MaxIP_norm_fs.parquet"
).resolve(strict=True)
scDINO_sc_profiles_path = pathlib.Path(
    "../../data/outputdir/apoptosis_timelapse/CLS_features/CLS_features_annotated.parquet"
).resolve(strict=True)
combined_profiles_path = pathlib.Path(
    "../../data/20231017ChromaLive_6hr_4ch_MaxIP_normalized_combined_data_feature_selected.parquet"
).resolve(strict=True)

output_path = pathlib.Path(f"../../data/umap/").resolve()
output_path.mkdir(parents=True, exist_ok=True)


# In[3]:


if data_mode == "CP":
    # read the data
    profiles_df = pd.read_parquet(CP_fs_sc_profiles_path)
elif data_mode == "scDINO":
    # read the data
    profiles_df = pd.read_parquet(scDINO_sc_profiles_path)
elif data_mode == "combined":
    # read the data
    profiles_df = pd.read_parquet(combined_profiles_path)
else:
    raise ValueError("data_mode must be either 'CP' or 'scDINO' or 'combined'")
print(profiles_df.shape)
# show all columns
pd.set_option("display.max_columns", None)
profiles_df.head()


# In[4]:


# filter the data and drop nan values
print(profiles_df.shape)
# drop nan values in non metadata columns
profiles_df = profiles_df.dropna(
    subset=profiles_df.columns[~profiles_df.columns.str.contains("Meta")]
)
print(profiles_df.shape)


# In[5]:


# get the metadata columns
metadata_cols = profiles_df.columns.str.contains("Metadata_")
metadata_df = profiles_df.loc[:, metadata_cols]
features_df = profiles_df.loc[:, ~metadata_cols]

# set the umap parameters
umap = umap.UMAP(
    n_neighbors=15,
    n_components=2,
    metric="euclidean",
    random_state=42,
    min_dist=0.1,
    n_epochs=500,
    learning_rate=1,
)

# fit the umap model
umap.fit(features_df)

# transform the data
umap_transformed = umap.transform(features_df)

# create a dataframe with the transformed data
umap_df = pd.DataFrame(
    umap_transformed, columns=["UMAP0", "UMAP1"], index=features_df.index
)

# combine the metadata and umap dataframes
umap_df = pd.concat([metadata_df, umap_df], axis=1)
print(umap_df.shape)
umap_df.head()


# In[ ]:


# save the umap dataframe
umap_df.to_parquet(f"../../data/umap/{data_mode}_umap_transformed.parquet")
