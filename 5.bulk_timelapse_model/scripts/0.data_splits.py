#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import pathlib

import numpy as np
import pandas as pd

# In[2]:


bulk_data_file_path = pathlib.Path(
    "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs_aggregated.parquet"
).resolve(strict=True)
whole_image_final_data_file_path = pathlib.Path(
    "../../data/CP_aggregated/endpoints/aggregated_whole_image.parquet"
).resolve(strict=True)
ground_truth_file_path = pathlib.Path(
    "../../1.ground_truth/data/0.ground_truth/ground_truth.csv"
).resolve(strict=True)
data_splits_dir = pathlib.Path("../data_splits/").resolve()
data_splits_dir.mkdir(parents=True, exist_ok=True)

# Load the data
bulk_df = pd.read_parquet(bulk_data_file_path)
ground_truth_df = pd.read_csv(ground_truth_file_path)
whole_image_final_df = pd.read_parquet(whole_image_final_data_file_path)
bulk_df["Metadata_dose"] = bulk_df["Metadata_dose"].astype("float64")
bulk_df["Metadata_Time"] = bulk_df["Metadata_Time"].astype("float64")
# get the final_timepoint only for the bulk data
bulk_df = bulk_df[bulk_df["Metadata_Time"] == bulk_df["Metadata_Time"].max()]
bulk_df.drop(columns=["Metadata_Time"], inplace=True)
bulk_df.head()


# In[3]:


# prepend "Terminal" to all columns in the whole image final dataframe
for col in whole_image_final_df.columns:
    if col == "Metadata_dose":
        continue
    if col == "Metadata_Well":
        continue
    whole_image_final_df.rename(columns={col: "Terminal_" + col}, inplace=True)


# In[4]:


print("Bulk data shape: ", bulk_df.shape)
print("Whole image final data shape: ", whole_image_final_df.shape)


# In[5]:


bulk_df = pd.merge(
    bulk_df,
    ground_truth_df[["Metadata_dose", "apoptosis"]],
    how="left",
    left_on="Metadata_dose",
    right_on="Metadata_dose",
)
gt = bulk_df.pop("apoptosis")
bulk_df.insert(3, "Metadata_apoptosis_ground_truth", gt)

bulk_df = pd.merge(
    bulk_df,
    whole_image_final_df,
    how="left",
    left_on=["Metadata_dose", "Metadata_Well"],
    right_on=["Metadata_dose", "Metadata_Well"],
)
bulk_df.head()


# In[6]:


dose_wells = bulk_df.copy()
dose_wells = dose_wells[["Metadata_dose", "Metadata_Well"]]
dose_wells = dose_wells.drop_duplicates()
dose_wells = dose_wells.reset_index(drop=True)


# In[ ]:


# there are 10 doses, with three wells each
# one well is needed for each dose for training
# select one well per dose
test_wells = []
for dose in dose_wells["Metadata_dose"].unique():
    wells = dose_wells[dose_wells["Metadata_dose"] == dose]["Metadata_Well"].tolist()
    selected_well = np.random.choice(wells, 1)[0]
    print(f"Selected well {selected_well} for dose {dose}")
    test_wells.append(str(selected_well))

train_wells = dose_wells[~dose_wells["Metadata_Well"].isin(test_wells)][
    "Metadata_Well"
].tolist()


# In[8]:


train_df = bulk_df[bulk_df["Metadata_Well"].isin(train_wells)]
test_df = bulk_df[bulk_df["Metadata_Well"].isin(test_wells)]
train_df = train_df.reset_index(drop=True)
test_df = test_df.reset_index(drop=True)
# write the train and test dataframes to parquet files
train_df_file_path = data_splits_dir / "train.parquet"
train_df.to_parquet(train_df_file_path, index=False)
test_df_file_path = data_splits_dir / "test.parquet"
test_df.to_parquet(test_df_file_path, index=False)


# In[9]:


print("Train data shape: ", train_df.shape)
train_df.head()


# In[10]:


print("Test data shape: ", test_df.shape)
test_df.head()


# In[ ]:


# make a df with the wells used for training and testing with their respective doses
test_well_df = pd.DataFrame(test_wells, columns=["Metadata_Well"])
train_well_df = pd.DataFrame(train_wells, columns=["Metadata_Well"])
test_well_df["data_split"] = "test"
train_well_df["data_split"] = "train"
train_test_well_df = pd.concat([train_well_df, test_well_df], axis=0)
# save the train test well df to a parquet file
train_test_well_file_path = data_splits_dir / "train_test_wells.parquet"
train_test_well_df.to_parquet(train_test_well_file_path, index=False)
