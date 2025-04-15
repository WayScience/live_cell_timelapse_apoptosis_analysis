#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
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
import tqdm

# Suppress all RuntimeWarnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import sys

sys.path.append("../utils")
from mAP_utils import run_mAP_across_time

# In[2]:


if not in_notebook:
    # setup the argument parser
    parser = argparse.ArgumentParser(
        description="Generate a map for differing cell counts"
    )

    parser.add_argument(
        "--percentage", type=float, help="Percentage of wells to use for the map file"
    )
    parser.add_argument("--seed", type=int, help="Seed for the random number generator")
    parser.add_argument(
        "--shuffle", action="store_true", help="Shuffle the order of the wells"
    )
    # parse the arguments
    args = parser.parse_args()
    percentage = args.percentage
    set_seed = args.seed
    shuffle = args.shuffle
else:
    percentage = 0.1
    set_seed = 0
    shuffle = False

output_file = pathlib.Path(
    f"../results/mAP_cell_percentages/{percentage}_{set_seed}_{shuffle}.parquet"
)
output_file.parent.mkdir(exist_ok=True, parents=True)


# In[3]:


data_file_path = pathlib.Path(
    "../../data/CP_feature_select/profiles/features_selected_profile.parquet"
).resolve(strict=True)
df = pd.read_parquet(data_file_path)

df.head()


# In[4]:


random.seed(set_seed)
subset_df = df.groupby(["Metadata_Time", "Metadata_Well"]).apply(
    lambda x: x.sample(frac=percentage, random_state=set_seed),
    include_groups=True,
)


# In[5]:


subset_df.reset_index(drop=True, inplace=True)
if shuffle:
    # permutate the data
    for col in subset_df.columns:
        subset_df[col] = np.random.permutation(subset_df[col])
metadata_cols = [cols for cols in subset_df.columns if "Metadata" in cols]
features_cols = [cols for cols in subset_df.columns if "Metadata" not in cols]
features_cols = features_cols + ["Metadata_number_of_singlecells"]
aggregate_df = pycytominer.aggregate(
    population_df=subset_df,
    strata=["Metadata_Well", "Metadata_Time"],
    features=features_cols,
    operation="median",
)


# In[6]:


metadata_df = subset_df[metadata_cols]
metadata_df = metadata_df.drop_duplicates(subset=["Metadata_Well", "Metadata_Time"])
metadata_df = metadata_df.reset_index(drop=True)
aggregate_df = pd.merge(
    metadata_df, aggregate_df, on=["Metadata_Well", "Metadata_Time"]
)
aggregate_df.head()


# In[ ]:


dict_of_map_dfs = run_mAP_across_time(
    aggregate_df,
    seed=set_seed,
    time_column="Metadata_Time",
    reference_column_name="Metadata_dose",
    reference_group=aggregate_df["Metadata_dose"].min(),
)


output_df = pd.concat(dict_of_map_dfs.values(), keys=dict_of_map_dfs.keys())
output_df.reset_index(inplace=True)
output_df.rename(columns={"level_0": "Metadata_Time"}, inplace=True)
# add the percentage of cells to the keys
output_df["percentage_of_cells"] = percentage
output_df["seed"] = set_seed
output_df["shuffle"] = shuffle
output_df.reset_index(drop=True, inplace=True)
output_df.to_parquet(output_file)
output_df.head()
