#!/usr/bin/env python
# coding: utf-8

# This is quite the complex data splitting procedure.
# The data is split into holdout data, training, validation, and testing.
# The training and validation data only contains single-cells that have grond truth labels at the terminal time point.
# While testing and holdout data contains celles that do and no not have ground truth labels at the terminal time point.

# In[1]:


import pathlib

import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean
from sklearn.model_selection import train_test_split

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


# write the cleaned dataframe to a parquet file
sc_file_path = pathlib.Path("../results/cleaned_sc_profile.parquet").resolve(
    strict=True
)
sc_endpoint_file_path = pathlib.Path(
    "../results/cleaned_endpoint_sc_profile.parquet"
).resolve(strict=True)

train_test_wells_file_path = pathlib.Path(
    "../../5.bulk_timelapse_model/data_splits/train_test_wells.parquet"
).resolve(strict=True)

sc_profile = pd.read_parquet(sc_file_path)
sc_endpoint_profile = pd.read_parquet(sc_endpoint_file_path)
train_test_wells = pd.read_parquet(train_test_wells_file_path)
print(f"sc_profile shape: {sc_profile.shape}")
print(f"sc_endpoint_profile shape: {sc_endpoint_profile.shape}")


# In[3]:


sc_endpoint_profile.head()


# In[4]:


sc_profile.head()


# In[5]:


# add a ground truth column to the sc_profile dataframe based on if the track id is in the endpoint profile
sc_profile["Metadata_ground_truth_present"] = (
    sc_profile["Metadata_sc_unique_track_id"]
    .isin(sc_endpoint_profile["Metadata_sc_unique_track_id"])
    .astype(bool)
)
sc_profile


# At this point there are two subsets to the dataset and will be split into the following datasplits:
# - Single-cells that have a single cell ground truth
#      - Holdout wells: 1/3 of wells
#     - Train: 80% of the all cells with a single cell ground truth from the non-holdout wells
#     - Validation:  10% of the all cells with a single cell ground truth from the non-holdout wells
#     - Test: 10% of the all cells with a single cell ground truth from the non-holdout wells
# - Single-cells that do not have a single cell ground truth
#     - Holdout wells: 1/3 of wells
#     - Test: 100% of the all cells with a single cell ground truth from the non-holdout wells
#

# ### hold out wells regardless of ground truth

# In[6]:


index_data_split_and_ground_truth_dict = {
    "index": [],
    "data_split": [],
    "ground_truth": [],
}


# In[7]:


# map the data_split to the sc_profile dataframe via the well
sc_profile["Metadata_data_split"] = sc_profile["Metadata_Well"].map(
    train_test_wells.set_index("Metadata_Well")["data_split"]
)
sc_profile.loc[sc_profile["Metadata_data_split"] == "test", "Metadata_data_split"] = (
    "well_holdout"
)
holdout_df = sc_profile.loc[sc_profile["Metadata_data_split"] == "well_holdout"]
index_data_split_and_ground_truth_dict["index"].append(holdout_df.index.tolist())
index_data_split_and_ground_truth_dict["data_split"].append(
    holdout_df["Metadata_data_split"].tolist()
)
index_data_split_and_ground_truth_dict["ground_truth"].append(
    holdout_df["Metadata_ground_truth_present"].tolist()
)
# get the non holdout wells
non_holdout_wells = sc_profile.loc[sc_profile["Metadata_data_split"] != "well_holdout"]
print(f"sc_profile shape after mapping data_split: {non_holdout_wells.shape}")
print(f"holdout_df shape: {holdout_df.shape}")


# ### Cells that have a single cell ground truth

# In[8]:


cell_wout_ground_truth_df = non_holdout_wells.loc[
    non_holdout_wells["Metadata_ground_truth_present"] == False
].copy()
cell_w_ground_truth_df = non_holdout_wells.loc[
    non_holdout_wells["Metadata_ground_truth_present"] == True
].copy()

print(f"cell_w_ground_truth_df shape: {cell_w_ground_truth_df.shape}")
print(f"cell_wout_ground_truth_df shape: {cell_wout_ground_truth_df.shape}")


# ##

# In[9]:


# split the data into 80, 10, 10 stratified by the well
train_sc_w_ground_truth_df, test_sc_w_ground_truth_df = train_test_split(
    cell_w_ground_truth_df,
    test_size=0.2,
    stratify=cell_w_ground_truth_df["Metadata_Well"],
    random_state=0,
)
test_sc_w_ground_truth_df, val_sc_w_ground_truth_df = train_test_split(
    test_sc_w_ground_truth_df,
    test_size=0.5,
    stratify=test_sc_w_ground_truth_df["Metadata_Well"],
    random_state=0,
)

train_sc_w_ground_truth_df["Metadata_data_split"] = "train"
train_sc_w_ground_truth_df["Metadata_ground_truth_present"] = True
val_sc_w_ground_truth_df["Metadata_data_split"] = "val"
val_sc_w_ground_truth_df["Metadata_ground_truth_present"] = True
test_sc_w_ground_truth_df["Metadata_data_split"] = "test"
test_sc_w_ground_truth_df["Metadata_ground_truth_present"] = True

print(f"train_sc_w_ground_truth_df shape: {train_sc_w_ground_truth_df.shape[0]}")
print(f"val_sc_w_ground_truth_df shape: {val_sc_w_ground_truth_df.shape[0]}")
print(f"test_sc_w_ground_truth_df shape: {test_sc_w_ground_truth_df.shape[0]}")
assert (
    train_sc_w_ground_truth_df.shape[0]
    + val_sc_w_ground_truth_df.shape[0]
    + test_sc_w_ground_truth_df.shape[0]
    == cell_w_ground_truth_df.shape[0]
)
assert (
    np.round(train_sc_w_ground_truth_df.shape[0] / cell_w_ground_truth_df.shape[0], 2)
    == 0.8
)
assert (
    np.round(val_sc_w_ground_truth_df.shape[0] / cell_w_ground_truth_df.shape[0], 2)
    == 0.1
)
assert (
    np.round(test_sc_w_ground_truth_df.shape[0] / cell_w_ground_truth_df.shape[0], 2)
    == 0.1
)

# add to records
index_data_split_and_ground_truth_dict["index"].append(
    train_sc_w_ground_truth_df.index.tolist()
)
index_data_split_and_ground_truth_dict["data_split"].append(
    train_sc_w_ground_truth_df["Metadata_data_split"].tolist()
)
index_data_split_and_ground_truth_dict["ground_truth"].append(
    train_sc_w_ground_truth_df["Metadata_ground_truth_present"].tolist()
)
index_data_split_and_ground_truth_dict["index"].append(
    val_sc_w_ground_truth_df.index.tolist()
)
index_data_split_and_ground_truth_dict["data_split"].append(
    val_sc_w_ground_truth_df["Metadata_data_split"].tolist()
)
index_data_split_and_ground_truth_dict["ground_truth"].append(
    val_sc_w_ground_truth_df["Metadata_ground_truth_present"].tolist()
)
index_data_split_and_ground_truth_dict["index"].append(
    test_sc_w_ground_truth_df.index.tolist()
)
index_data_split_and_ground_truth_dict["data_split"].append(
    test_sc_w_ground_truth_df["Metadata_data_split"].tolist()
)
index_data_split_and_ground_truth_dict["ground_truth"].append(
    test_sc_w_ground_truth_df["Metadata_ground_truth_present"].tolist()
)


# #### Non tracked cells

# In[10]:


cell_wout_ground_truth_df["Metadata_data_split"] = "test"
# add to records
index_data_split_and_ground_truth_dict["index"].append(
    cell_wout_ground_truth_df.index.tolist()
)
index_data_split_and_ground_truth_dict["data_split"].append(
    cell_wout_ground_truth_df["Metadata_data_split"].tolist()
)
index_data_split_and_ground_truth_dict["ground_truth"].append(
    cell_wout_ground_truth_df["Metadata_ground_truth_present"].tolist()
)
print(f"test_sc_wo_ground_truth_df shape: {cell_wout_ground_truth_df.shape[0]}")


# ### Fetch the indices from each ground truth and data split and add the status back to sc_profile

# In[11]:


# flatten each list in the dictionary
import itertools

for key in index_data_split_and_ground_truth_dict.keys():
    index_data_split_and_ground_truth_dict[key] = list(
        itertools.chain.from_iterable(index_data_split_and_ground_truth_dict[key])
    )
data_split_data_df = pd.DataFrame.from_dict(
    index_data_split_and_ground_truth_dict,
    orient="columns",
)
assert data_split_data_df.shape[0] == sc_profile.shape[0]


# In[12]:


# sort the dataframe by index
data_split_data_df.sort_values(
    by=["index"],
    inplace=True,
)
# make the index the index column in data_split_data_df
data_split_data_df.set_index("index", inplace=True)
data_split_data_df.reset_index(drop=False, inplace=True)
data_split_data_df.head()


# In[13]:


# addthe data_split and ground truth columns to the sc_profile dataframe by index
sc_profile_with_data_splits_df = pd.concat(
    [sc_profile, data_split_data_df],
    axis=1,
)
sc_profile_with_data_splits_df.drop(
    columns=["Metadata_data_split", "Metadata_ground_truth_present"],
    inplace=True,
)
sc_profile_with_data_splits_df.rename(
    columns={
        "data_split": "Metadata_data_split",
        "ground_truth": "Metadata_ground_truth_present",
    },
    inplace=True,
)


# In[14]:


# final breakdown of the data
train_gt = sc_profile_with_data_splits_df[
    sc_profile_with_data_splits_df["Metadata_data_split"] == "train"
].copy()
train_gt = train_gt[train_gt["Metadata_ground_truth_present"] == True].copy()
val_gt = sc_profile_with_data_splits_df[
    sc_profile_with_data_splits_df["Metadata_data_split"] == "val"
].copy()
val_gt = val_gt[val_gt["Metadata_ground_truth_present"] == True].copy()
test_gt = sc_profile_with_data_splits_df[
    sc_profile_with_data_splits_df["Metadata_data_split"] == "test"
].copy()
test_gt = test_gt[test_gt["Metadata_ground_truth_present"] == True].copy()
test_wo_gt = sc_profile_with_data_splits_df[
    sc_profile_with_data_splits_df["Metadata_data_split"] == "test"
].copy()
test_wo_gt = test_wo_gt[test_wo_gt["Metadata_ground_truth_present"] == False].copy()
holdout_w_gt = sc_profile_with_data_splits_df[
    sc_profile_with_data_splits_df["Metadata_data_split"] == "well_holdout"
].copy()
holdout_w_gt = holdout_w_gt[
    holdout_w_gt["Metadata_ground_truth_present"] == True
].copy()
holdout_wo_gt = sc_profile_with_data_splits_df[
    sc_profile_with_data_splits_df["Metadata_data_split"] == "well_holdout"
].copy()
holdout_wo_gt = holdout_wo_gt[
    holdout_wo_gt["Metadata_ground_truth_present"] == False
].copy()
# assertion time :)
assert sc_profile_with_data_splits_df.shape[0] == sc_profile.shape[0]
assert (
    sc_profile_with_data_splits_df.shape[0]
    == train_gt.shape[0]
    + val_gt.shape[0]
    + test_gt.shape[0]
    + test_wo_gt.shape[0]
    + holdout_w_gt.shape[0]
    + holdout_wo_gt.shape[0]
)


# In[15]:


print(f"train_gt shape: {train_gt.shape[0]}")
print(f"val_gt shape: {val_gt.shape[0]}")
print(f"test_gt shape: {test_gt.shape[0]}")
print(f"test_wo_gt shape: {test_wo_gt.shape[0]}")
print(f"holdout_w_gt shape: {holdout_w_gt.shape[0]}")
print(f"holdout_wo_gt shape: {holdout_wo_gt.shape[0]}")


# In[16]:


# write the data splits to a parquet file
data_split_file_path = pathlib.Path("../results/data_splits.parquet").resolve()
data_split_data_df.to_parquet(
    data_split_file_path,
    index=False,
)
