#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import random

import numpy as np
import pandas as pd
from copairs import map
from copairs.matching import assign_reference_index

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False
import warnings

import pycytominer.aggregate

# Suppress all RuntimeWarnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import itertools
import sys

sys.path.append("../utils")
from mAP_utils import run_mAP_across_time

# In[2]:


data_file_path = pathlib.Path(
    "../../data/CP_aggregated/profiles/aggregated_profile.parquet"
).resolve(strict=True)
# load the dose information
well_dose_df = pd.read_parquet(
    pathlib.Path(
        "../../data/CP_feature_select/profiles/features_selected_profile.parquet"
    ).resolve(strict=True),
    columns=["Metadata_Well", "Metadata_dose"],
).drop_duplicates()


pathlib.Path("../results").mkdir(exist_ok=True)
# output filepath
output_file_path = pathlib.Path("../results/mAP_across_channels.parquet").resolve()
aggregate_df = pd.read_parquet(data_file_path)
print(aggregate_df.shape)
# drop rows with nan in features (non Metadata)
aggregate_df = aggregate_df.dropna(
    subset=aggregate_df.filter(regex="^(?!Metadata)").columns
)
# merge the dose information with the profile data
aggregate_df = aggregate_df.merge(
    well_dose_df,
    left_on="Metadata_Well",
    right_on="Metadata_Well",
    how="left",
)
dose = aggregate_df.pop("Metadata_dose")
aggregate_df.insert(
    1,
    "Metadata_dose",
    dose,
)
print(aggregate_df.shape)
aggregate_df.head()


# In[3]:


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


# In[4]:


# create a dictionary of loadable features for each channel
loadable_features = {}
for channel in df.Channel.unique():
    channel_df = df.query("Channel == @channel")
    loadable_features[channel] = channel_df.features.tolist()
loadable_features["None"] = non_feature_cols
loadable_features["All"] = [
    cols for cols in aggregate_df.columns if "Metadata" not in cols
]


# In[5]:


unique_channels = df["Channel"].unique().tolist()
# unique_channels = unique_channels + ['None']
# channel combinations
channel_combinations = []
for i in range(1, len(unique_channels) + 1):
    channel_combinations.extend(list(itertools.combinations(unique_channels, i)))
channel_combinations = [list(comb) for comb in channel_combinations]
channel_combinations = channel_combinations + [["None"]]
channel_combinations


# In[6]:


channel_combo_output_dict = {}
for shuffle in [False, True]:
    for channel_combination in channel_combinations:
        features_to_load = []
        features_to_load.append(metadata_cols)
        if len(channel_combination) == 5:
            channel_combination = ["All"]
        elif channel_combination == ["None"]:
            features_to_load.append(loadable_features["None"])
        else:
            for channel in channel_combination:
                features_to_load.append(loadable_features[channel])
            features_to_load.append(loadable_features["None"])
        features_to_load = list(itertools.chain(*features_to_load))
        temporary_df = aggregate_df[features_to_load]
        if shuffle:
            shuffle = "shuffled"
            random.seed(0)
            # permutate the data
            for col in temporary_df.columns:
                if "Metadata_Time" in col or "Metadata_dose" in col:
                    pass
                else:
                    temporary_df[col] = np.random.permutation(temporary_df[col])
        else:
            shuffle = "non_shuffled"

        dict_of_map_dfs = run_mAP_across_time(
            temporary_df,
            seed=0,
            time_column="Metadata_Time",
            reference_column_name="Metadata_dose",
            reference_group=aggregate_df["Metadata_dose"].min(),
        )
        df = pd.concat(dict_of_map_dfs.values(), keys=dict_of_map_dfs.keys())
        df.reset_index(inplace=True)
        df["Channel"] = "_".join(channel_combination)
        df["shuffle"] = shuffle
        channel_combo_output_dict["_".join(channel_combination)] = df


# In[7]:


final_df = pd.concat(
    channel_combo_output_dict.values(), keys=channel_combo_output_dict.keys()
)
final_df.reset_index(inplace=True, drop=True)
# save the output to a file
final_df.to_parquet("../results/mAP_across_channels.parquet")
final_df.head()
