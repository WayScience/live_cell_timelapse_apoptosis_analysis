#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import joblib
import numpy as np
import pandas as pd
import pycytominer

# In[2]:


model_file_dir = pathlib.Path(
    "../models/multi_regression_model_ntrees_1000.joblib"
).resolve()
shuffled_model_file_dir = pathlib.Path(
    "../models/shuffled_multi_regression_model_ntrees_1000.joblib"
).resolve()
train_test_wells_path = pathlib.Path(
    "../data_splits/train_test_wells.parquet"
).resolve()

predictions_save_path = pathlib.Path(
    "../results/predicted_terminal_profiles_from_all_time_points.parquet"
).resolve()

profile_data_path = pathlib.Path(
    "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs_aggregated.parquet"
).resolve()
terminal_column_names = pathlib.Path("../results/terminal_columns.txt").resolve(
    strict=True
)
terminal_column_names = [
    line.strip() for line in terminal_column_names.read_text().splitlines()
]

data_split_df = pd.read_parquet(train_test_wells_path)
df = pd.read_parquet(profile_data_path)
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


# In[3]:


# map the train/test wells to the aggregate data
aggregate_df["Metadata_data_split"] = aggregate_df["Metadata_Well"].map(
    data_split_df.set_index("Metadata_Well")["data_split"]
)
data_split = aggregate_df.pop("Metadata_data_split")
aggregate_df.insert(0, "Metadata_data_split", data_split)
aggregate_df["Metadata_Time"] = aggregate_df["Metadata_Time"].astype(float)
aggregate_df["Metadata_data_split"].unique()


# In[4]:


aggregate_df.head(15)


# In[5]:


# if the data_split is train and the time is not 12 then set to non_trained_pair
aggregate_df["Metadata_data_split"] = aggregate_df.apply(
    lambda x: (
        "non_trained_pair"
        if (x["Metadata_data_split"] == "train" and x["Metadata_Time"] != 12.0)
        else x["Metadata_data_split"]
    ),
    axis=1,
)


# In[6]:


# load the model
model = joblib.load(model_file_dir)

metadata_columns = [x for x in aggregate_df.columns if "Metadata_" in x]
# remove metadata columns
features = aggregate_df.drop(columns=metadata_columns)
metadata_df = aggregate_df[metadata_columns]
# predict the terminal feature space
predictions = model.predict(features)
predictions_df = pd.DataFrame(predictions, columns=terminal_column_names)
# insert the metadata columns
for col in metadata_columns:
    predictions_df.insert(0, col, metadata_df[col])
predictions_df["shuffled"] = False


# In[7]:


# load the model
shuffled_model = joblib.load(shuffled_model_file_dir)

metadata_columns = [x for x in aggregate_df.columns if "Metadata_" in x]
shuffled_profile_df = aggregate_df.copy()
for col in shuffled_profile_df.columns:
    shuffled_profile_df[col] = np.random.permutation(shuffled_profile_df[col])
# remove metadata columns
features = shuffled_profile_df.drop(columns=metadata_columns)
metadata_df = aggregate_df[metadata_columns]


# predict the terminal feature space
predictions = shuffled_model.predict(features)
shuffled_predictions_df = pd.DataFrame(predictions, columns=terminal_column_names)
# insert the metadata columns
for col in metadata_columns:
    shuffled_predictions_df.insert(0, col, metadata_df[col])
shuffled_predictions_df["shuffled"] = True


# In[8]:


final_predictions_df = pd.concat([predictions_df, shuffled_predictions_df], axis=0)
# save the predictions
final_predictions_df.to_parquet(predictions_save_path, index=False)
final_predictions_df
