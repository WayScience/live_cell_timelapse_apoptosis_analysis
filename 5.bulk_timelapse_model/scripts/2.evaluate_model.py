#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import joblib
import numpy as np
import pandas as pd

# In[2]:


# load the training data
profile_file_dir = pathlib.Path(
    "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs_aggregated.parquet"
).resolve(strict=True)
model_file_dir = pathlib.Path("../models/multi_regression_model.joblib").resolve(
    strict=True
)
shuffled_model_file_dir = pathlib.Path(
    "../models/shuffled_multi_regression_model.joblib"
).resolve(strict=True)
terminal_column_names = pathlib.Path("../results/terminal_columns.txt").resolve(
    strict=True
)
predictions_save_path = pathlib.Path(
    "../results/predicted_terminal_profiles.parquet"
).resolve()
terminal_column_names = [
    line.strip() for line in terminal_column_names.read_text().splitlines()
]
results_dir = pathlib.Path("../results/").resolve()
results_dir.mkdir(parents=True, exist_ok=True)
profile_df = pd.read_parquet(profile_file_dir)
print(profile_df.shape)
profile_df.head()


# ## Get the non-shuffled predictions

# In[3]:


# load the model
model = joblib.load(model_file_dir)
shuffled_model = joblib.load(shuffled_model_file_dir)

metadata_columns = [x for x in profile_df.columns if "Metadata_" in x]
# remove metadata columns
features = profile_df.drop(columns=metadata_columns)
metadata_df = profile_df[metadata_columns]
# predict the terminal feature space
predictions = model.predict(features)
predictions_df = pd.DataFrame(predictions, columns=terminal_column_names)
# insert the metadata columns
for col in metadata_columns:
    predictions_df.insert(0, col, metadata_df[col])
predictions_df["shuffled"] = False


# ## Get the shuffled predictions

# In[4]:


# load the model
shuffled_model = joblib.load(shuffled_model_file_dir)

metadata_columns = [x for x in profile_df.columns if "Metadata_" in x]
shuffled_profile_df = profile_df.copy()
for col in shuffled_profile_df.columns:
    shuffled_profile_df[col] = np.random.permutation(shuffled_profile_df[col])
# remove metadata columns
features = shuffled_profile_df.drop(columns=metadata_columns)
metadata_df = profile_df[metadata_columns]


# predict the terminal feature space
predictions = shuffled_model.predict(features)
shuffled_predictions_df = pd.DataFrame(predictions, columns=terminal_column_names)
# insert the metadata columns
for col in metadata_columns:
    shuffled_predictions_df.insert(0, col, metadata_df[col])
shuffled_predictions_df["shuffled"] = True


# In[5]:


final_predictions_df = pd.concat([predictions_df, shuffled_predictions_df], axis=0)
# save the predictions
final_predictions_df.to_parquet(predictions_save_path, index=False)
final_predictions_df
