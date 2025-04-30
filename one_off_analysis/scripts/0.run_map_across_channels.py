#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import random
import sys
import warnings

import numpy as np
import pandas as pd
import pycytominer.aggregate

# Suppress all RuntimeWarnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import itertools

sys.path.append("../utils")
from mAP_utils import run_mAP_across_time

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# In[2]:


data_file_path = pathlib.Path(
    "../../data/CP_feature_select/profiles/features_selected_profile.parquet"
).resolve(strict=True)
# load the dose information


pathlib.Path("../results").mkdir(exist_ok=True)
# output filepath
output_file_path = pathlib.Path("../results/mAP_across_channels.parquet").resolve()
df = pd.read_parquet(data_file_path)

df.head()


# In[3]:


metadata_cols = [cols for cols in df.columns if "Metadata" in cols]
features_cols = [cols for cols in df.columns if "Metadata" not in cols]
features_cols = features_cols
aggregate_df = pycytominer.aggregate(
    population_df=df,
    strata=["Metadata_Well", "Metadata_Time"],
    features=features_cols,
    operation="median",
)
metadata_df = df[metadata_cols]
metadata_df = metadata_df.drop_duplicates(subset=["Metadata_Well", "Metadata_Time"])
metadata_df = metadata_df.reset_index(drop=True)
aggregate_df = pd.merge(
    metadata_df, aggregate_df, on=["Metadata_Well", "Metadata_Time"]
)
print(aggregate_df.shape)
aggregate_df.head()


# In[4]:


# get the non metadata columns
metadata_cols = [cols for cols in aggregate_df.columns if "Metadata" in cols]
features_cols = [cols for cols in aggregate_df.columns if "Metadata" not in cols]
# drop features that contain AreaShape in the feature name
non_feature_cols = [cols for cols in features_cols if "AreaShape" in cols]
non_feature_cols = (
    non_feature_cols
    + [col for col in features_cols if "Location" in col]
    + [col for col in features_cols if "Neighbors" in col]
)
correlation_cols = [col for col in features_cols if "Correlation" in col]
features_cols = [cols for cols in features_cols if cols not in non_feature_cols]
features_cols = [cols for cols in features_cols if cols not in correlation_cols]
df = pd.DataFrame(features_cols, columns=["features"])

df[
    [
        "Compartment",
        "Feature_type",
        "Measurement",
        "Channel",
        "misc1",
        "misc2",
        "misc3",
        "misc4",
        "misc5",
    ]
] = df["features"].str.split("_", expand=True)
# Update Channel column where it contains "CL"

df.loc[df["Channel"].str.contains("CL", na=False), "Channel"] = (
    df["Channel"] + "_" + df["misc1"]
)
# Update Channel column where it contains "CL_488"
df.loc[df["Channel"].str.contains("CL_488", na=False), "Channel"] = (
    df["Channel"] + "_" + df["misc2"]
)
df.drop(columns=["misc1", "misc2", "misc3", "misc4", "misc5"], inplace=True)
df.head()


# In[5]:


# create a dictionary of loadable features for each channel
loadable_features = {}
for channel in df.Channel.unique():
    channel_df = df.query("Channel == @channel")
    loadable_features[channel] = channel_df.features.tolist()
loadable_features["None"] = non_feature_cols
loadable_features["All"] = [
    cols for cols in aggregate_df.columns if "Metadata" not in cols
]


# In[6]:


unique_channels = df["Channel"].unique().tolist()
# unique_channels = unique_channels + ['None']
# channel combinations
channel_combinations = []
for i in range(1, len(unique_channels) + 1):
    channel_combinations.extend(list(itertools.combinations(unique_channels, i)))
channel_combinations = [list(comb) for comb in channel_combinations]
channel_combinations = channel_combinations + [["None"]]
channel_combinations


# In[7]:


channel_combo_output_dict = {}
# loop through the channel combinations and shuffle
for shuffle in [False, True]:
    for channel_combination in channel_combinations:
        features_to_load = []
        features_to_load.append(metadata_cols)
        # rename the all channels to all
        if len(channel_combination) == 5:
            channel_combination = ["All"]
        elif channel_combination == ["None"]:
            features_to_load.append(loadable_features["None"])
        else:
            for channel in channel_combination:
                features_to_load.append(loadable_features[channel])
            features_to_load.append(loadable_features["None"])
        # flatten the list
        features_to_load = list(itertools.chain(*features_to_load))
        temporary_df = aggregate_df[features_to_load]
        # shuffle the data
        if shuffle == True:
            print("shuffled")
            shuffle_status = "shuffled"
            random.seed(0)
            # permutate the data
            for col in temporary_df.columns:
                if "Metadata_Time" in col or "Metadata_dose" in col:
                    continue
                else:
                    temporary_df[col] = np.random.permutation(temporary_df[col])
        else:
            print("not shuffled")
            shuffle_status = "non_shuffled"

        # run mAP with the 0 dose as the reference
        dict_of_map_dfs = run_mAP_across_time(
            temporary_df,
            seed=0,
            time_column="Metadata_Time",
            reference_column_name="Metadata_dose",
            reference_group=aggregate_df["Metadata_dose"].min(),
        )
        # concat and rename the columns
        df = pd.concat(dict_of_map_dfs.values(), keys=dict_of_map_dfs.keys())
        df.reset_index(inplace=True)
        df.rename(
            columns={
                "level_0": "Metadata_Time",
            },
            inplace=True,
        )
        df["Channel"] = "_".join(channel_combination)
        df["shuffle"] = shuffle_status
        channel_combo_output_dict[
            f"{'_'.join(channel_combination)}_{shuffle_status}"
        ] = df


# In[8]:


final_df = pd.concat(
    channel_combo_output_dict.values(), keys=channel_combo_output_dict.keys()
)
final_df.reset_index(inplace=True, drop=True)
# save the output to a file
final_df.to_parquet(output_file_path)
final_df.head()
