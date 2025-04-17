#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
from typing import Tuple

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import (
    explained_variance_score,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor

# ## Import the data

# In[2]:


# load the training data
train_data_file_path = pathlib.Path("../data_splits/train.parquet").resolve(strict=True)
test_data_file_path = pathlib.Path("../data_splits/test.parquet").resolve(strict=True)
model_dir = pathlib.Path("../models/").resolve()
model_dir.mkdir(parents=True, exist_ok=True)
results_dir = pathlib.Path("../results/").resolve()
results_dir.mkdir(parents=True, exist_ok=True)
train_df = pd.read_parquet(train_data_file_path)
test_df = pd.read_parquet(test_data_file_path)
train_df.head()


# In[3]:


metadata_columns = [x for x in train_df.columns if "Metadata" in x]
terminal_columns = [x for x in train_df.columns if "Terminal" in x]


def x_y_data_separator(
    df: pd.DataFrame,
    y_columns: list,
    metadata_columns: list,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df_shuffled = df.copy()
    for col in df_shuffled.columns:
        # permute the columns
        df_shuffled[col] = np.random.permutation(df_shuffled[col])
    metadata = df[metadata_columns]
    df.drop(columns=metadata_columns, inplace=True)
    X = df.drop(columns=y_columns)
    y = df[y_columns]

    metadata_shuffled = df_shuffled[metadata_columns]
    df_shuffled.drop(columns=metadata_columns, inplace=True)
    X_shuffled = df_shuffled.drop(columns=y_columns)
    y_shuffled = df_shuffled[y_columns]
    return X, y, metadata, X_shuffled, y_shuffled, metadata_shuffled


(
    train_X,
    train_y,
    train_metadata,
    train_shuffled_X,
    train_shuffled_y,
    train_metadata_shuffled,
) = x_y_data_separator(
    df=train_df, y_columns=terminal_columns, metadata_columns=metadata_columns
)

(
    test_X,
    test_y,
    test_metadata,
    test_shuffled_X,
    test_shuffled_y,
    test_metadata_shuffled,
) = x_y_data_separator(
    df=test_df, y_columns=terminal_columns, metadata_columns=metadata_columns
)

# check the shape of the data
print(f"train_X shape: {train_X.shape}, train_y shape: {train_y.shape}")
print(
    f"train_shuffled_X shape: {train_shuffled_X.shape}, train_shuffled_y shape: {train_shuffled_y.shape}"
)

print(f"test_X shape: {test_X.shape}, test_y shape: {test_y.shape}")
print(
    f"test_shuffled_X shape: {test_shuffled_X.shape}, test_shuffled_y shape: {test_shuffled_y.shape}"
)


# ## Model training

# In[4]:


# train the multi-output regression model
model = MultiOutputRegressor(
    RandomForestRegressor(
        n_estimators=1000,
        random_state=0,
    )
)

model.fit(train_X, train_y)
# save the model

model_file_path = pathlib.Path("../models/multi_regression_model.joblib").resolve()
joblib.dump(model, model_file_path)


shuffled_model = MultiOutputRegressor(
    RandomForestRegressor(
        n_estimators=1000,
        random_state=0,
    )
)
shuffled_model.fit(train_shuffled_X, train_shuffled_y)
# save the model
shuffled_model_file_path = pathlib.Path(
    "../models/shuffled_multi_regression_model.joblib"
).resolve()
joblib.dump(shuffled_model, shuffled_model_file_path)


# ## Model Evaluation

# In[5]:


dict_of_train_tests = {
    "train": {
        "X": train_X,
        "y": train_y,
        "metadata": train_metadata,
        "model_path": model_file_path,
    },
    "train_shuffled": {
        "X": train_shuffled_X,
        "y": train_shuffled_y,
        "metadata": train_metadata_shuffled,
        "model_path": shuffled_model_file_path,
    },
    "test": {
        "X": test_X,
        "y": test_y,
        "metadata": test_metadata,
        "model_path": model_file_path,
    },
    "test_shuffled": {
        "X": test_shuffled_X,
        "y": test_shuffled_y,
        "metadata": test_metadata_shuffled,
        "model_path": shuffled_model_file_path,
    },
}


# In[6]:


output_dict_of_dfs = {"prediction_stats_df": [], "predictions_df": []}
for split in dict_of_train_tests.keys():
    if "shuffle" in split:
        shuffle = True
    else:
        shuffle = False
    X = dict_of_train_tests[split]["X"]
    y = dict_of_train_tests[split]["y"]
    metadata = dict_of_train_tests[split]["metadata"]
    model = joblib.load(dict_of_train_tests[split]["model_path"])

    # make predictions on the training data
    y_pred = model.predict(X)
    # calculate the mean absolute error
    mae = mean_absolute_error(y, y_pred)
    # calculate the mse
    mse = mean_squared_error(y, y_pred)
    # calculate the r2 score
    r2 = r2_score(y, y_pred)
    # calculate the explained variance score
    evs = explained_variance_score(y, y_pred)

    prediction_stats_df = pd.DataFrame(
        {"MAE": [mae], "MSE": [mse], "R2": [r2], "EVS": [evs]}
    )

    prediction_stats_df["data_split"] = split
    prediction_stats_df["shuffled"] = shuffle
    output_dict_of_dfs["prediction_stats_df"].append(prediction_stats_df)

    predictions_df = pd.DataFrame(y_pred, columns=terminal_columns)
    predictions_df.insert(0, "Metadata_data_split", split)
    predictions_df.insert(1, "Metadata_shuffled", shuffle)
    # add the metadata columns to the predictions_df
    for col in metadata.columns:
        predictions_df.insert(2, col, metadata[col])
    output_dict_of_dfs["predictions_df"].append(predictions_df)


prediction_stats_df = pd.concat(output_dict_of_dfs["prediction_stats_df"], axis=0)
predictions_df = pd.concat(output_dict_of_dfs["predictions_df"], axis=0)
print(prediction_stats_df.shape)
print(predictions_df.shape)


# In[8]:


# save the final training results
prediction_stats_df_file_path = results_dir / "prediction_stats_df.parquet"
prediction_stats_df.to_parquet(prediction_stats_df_file_path, index=False)
predictions_df_file_path = results_dir / "predictions_df_final_timepoint.parquet"
predictions_df.to_parquet(predictions_df_file_path, index=False)
# write the terminal column names to a file
terminal_columns_file_path = results_dir / "terminal_columns.txt"
with open(terminal_columns_file_path, "w") as f:
    for col in terminal_columns:
        f.write(f"{col}\n")
