#!/usr/bin/env python
# coding: utf-8

# This notebook performs a fuzzy match on single-cells from the profile data to the endpoint data.
# This is necessary because  the endpoint data was not included in the tracking module.
# Further this notebook provides information about how long the cell track is.

# In[1]:


import pathlib
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean

try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False
if in_notebook:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm


# In[2]:


sc_profile_file_path = pathlib.Path(
    "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs.parquet"
).resolve(strict=True)
endpoint_sc_profile_file_path = pathlib.Path(
    "../../data/CP_feature_select/endpoints/features_selected_profile.parquet"
).resolve(strict=True)
sc_profile_df = pd.read_parquet(sc_profile_file_path)
endpoint_sc_profile_df = pd.read_parquet(endpoint_sc_profile_file_path)
endpoint_sc_profile_df["Metadata_Well_FOV"] = (
    endpoint_sc_profile_df["Metadata_Well"].astype(str)
    + "_"
    + endpoint_sc_profile_df["Metadata_FOV"].astype(str)
)
print(endpoint_sc_profile_df.shape)
sc_profile_df["Metadata_sc_unique_track_id"] = (
    sc_profile_df["Metadata_Well"].astype(str)
    + "_"
    + sc_profile_df["Metadata_FOV"].astype(str)
    + "_"
    + sc_profile_df["Metadata_track_id"].astype(str)
)
sc_profile_df["Metadata_Well_FOV"] = (
    sc_profile_df["Metadata_Well"].astype(str)
    + "_"
    + sc_profile_df["Metadata_FOV"].astype(str)
)
print(sc_profile_df.shape)
sc_profile_df.head()


# In[3]:


# drop all nan values in the location columns
endpoint_sc_profile_df = endpoint_sc_profile_df.dropna(
    subset=["Metadata_Nuclei_Location_Center_X", "Metadata_Nuclei_Location_Center_Y"]
)


# In[4]:


last_time_point_df = sc_profile_df.loc[
    sc_profile_df["Metadata_Well_FOV"].isin(
        endpoint_sc_profile_df["Metadata_Well_FOV"].unique()
    )
]


# In[ ]:


# chunk the dataframe by well_fov so that there are fewer indexed records to search at once during fuzzy matching
dict_of_sc_well_fovs = {}
for well_fov in last_time_point_df["Metadata_Well_FOV"].unique():
    dict_of_sc_well_fovs[well_fov] = last_time_point_df[
        last_time_point_df["Metadata_Well_FOV"] == well_fov
    ].copy()
    # get only the last timepoint of the track

    dict_of_sc_well_fovs[well_fov].reset_index(drop=True, inplace=True)
dict_of_sc_well_fovs_endpoint = {}
for well_fov in endpoint_sc_profile_df["Metadata_Well_FOV"].unique():
    dict_of_sc_well_fovs_endpoint[well_fov] = endpoint_sc_profile_df[
        endpoint_sc_profile_df["Metadata_Well_FOV"] == well_fov
    ].copy()
    dict_of_sc_well_fovs_endpoint[well_fov].reset_index(drop=True, inplace=True)


# In[7]:


start_time = time.time()


# In[8]:


for well_fov in tqdm(list(dict_of_sc_well_fovs.keys()), desc="Processing Well-FOVs"):
    for i, row in tqdm(
        dict_of_sc_well_fovs[well_fov].iterrows(),
        total=len(dict_of_sc_well_fovs[well_fov]),
        desc="Outer Loop",
        leave=False,
    ):
        for j, row2 in dict_of_sc_well_fovs_endpoint[well_fov].iterrows():
            # check that the well_fov is the same
            if row["Metadata_Well_FOV"] == row2["Metadata_Well_FOV"]:
                distance = abs(
                    euclidean(
                        [
                            row["Metadata_Nuclei_Location_Center_X"],
                            row["Metadata_Nuclei_Location_Center_Y"],
                        ],
                        [
                            row2["Metadata_Nuclei_Location_Center_X"],
                            row2["Metadata_Nuclei_Location_Center_Y"],
                        ],
                    )
                )
                if distance < 10:
                    dict_of_sc_well_fovs_endpoint[well_fov].at[
                        j, "Metadata_sc_unique_track_id"
                    ] = row["Metadata_sc_unique_track_id"]


# In[9]:


print("Fuzzy matching completed!")
print(f"Took {time.time() - start_time} seconds")
print(f"Took {round((time.time() - start_time) / 60, 2)} minutes")
print(f"Took {round((time.time() - start_time) / 3600, 2)} hours")


# In[10]:


sc_well_fovs_endpoint_df = pd.concat(
    dict_of_sc_well_fovs_endpoint.values(), ignore_index=True
)
# drop the rows where Metadata_sc_unique_track_id is NaN
sc_well_fovs_endpoint_df = sc_well_fovs_endpoint_df.dropna(
    subset=["Metadata_sc_unique_track_id"]
)
print(sc_well_fovs_endpoint_df.shape)
sc_well_fovs_endpoint_df.reset_index(drop=True, inplace=True)
sc_well_fovs_endpoint_df["Metadata_Time"] = 13.0
sc_well_fovs_endpoint_df.head()


# In[ ]:


# map the value counts to a new column for each Metadata_sc_unique_track_id
sc_profile_df["Metadata_sc_unique_track_id_count"] = sc_profile_df[
    "Metadata_sc_unique_track_id"
].map(sc_profile_df["Metadata_sc_unique_track_id"].value_counts())
sc_profile_df.head()


# In[ ]:


# write the cleaned dataframe to a parquet file
output_sc_file_path = pathlib.Path("../results/cleaned_sc_profile.parquet").resolve(
    strict=False
)
output_sc_endpoint_file_path = pathlib.Path(
    "../results/cleaned_endpoint_sc_profile.parquet"
).resolve(strict=False)
output_sc_file_path.parent.mkdir(parents=True, exist_ok=True)

# we save the two profiles separately because they have different feature spaces

sc_profile_df.to_parquet(output_sc_file_path, index=False)
sc_well_fovs_endpoint_df.to_parquet(output_sc_endpoint_file_path, index=False)
