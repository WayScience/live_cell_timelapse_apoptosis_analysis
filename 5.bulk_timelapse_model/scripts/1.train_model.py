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
from sklearn.model_selection import KFold, cross_val_score, train_test_split
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


def shuffle_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Shuffle the data in the DataFrame.
    """
    df_shuffled = df.copy()
    for col in df_shuffled.columns:
        # permute the columns
        df_shuffled[col] = np.random.permutation(df_shuffled[col])
    return df_shuffled


def x_y_data_separator(
    df: pd.DataFrame,
    y_columns: list,
    metadata_columns: list,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Separate the data into X, y, and metadata.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to separate. ASSUMPTION:
            The metadata columns contain the string "Metadata" and the y columns contain the string "Terminal".
            The column names are passed in as lists.
    y_columns : list
        The y columns to separate.
    metadata_columns : list
        The metadata columns to separate.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Three DataFrames: X, y, and metadata.
    """
    metadata = df[metadata_columns]
    df.drop(columns=metadata_columns, inplace=True)
    X = df.drop(columns=y_columns)
    y = df[y_columns]
    return X, y, metadata


shuffled_train_df = train_df.copy()
shuffled_train_df = shuffle_data(shuffled_train_df)
shuffled_test_df = test_df.copy()
shuffled_test_df = shuffle_data(shuffled_test_df)

# split the data into train and test sets
# train
(train_X, train_y, train_metadata) = x_y_data_separator(
    df=train_df, y_columns=terminal_columns, metadata_columns=metadata_columns
)
(train_shuffled_X, train_shuffled_y, train_metadata_shuffled) = x_y_data_separator(
    df=shuffled_train_df, y_columns=terminal_columns, metadata_columns=metadata_columns
)

# test
(test_X, test_y, test_metadata) = x_y_data_separator(
    df=test_df, y_columns=terminal_columns, metadata_columns=metadata_columns
)
(test_shuffled_X, test_shuffled_y, test_metadata_shuffled) = x_y_data_separator(
    df=shuffled_test_df, y_columns=terminal_columns, metadata_columns=metadata_columns
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


# In[4]:


dict_of_train_tests = {
    "train": {
        "X": train_X,
        "y": train_y,
        "metadata": train_metadata,
        "model_path": [],
        "n_estimators": [],
    },
    "train_shuffled": {
        "X": train_shuffled_X,
        "y": train_shuffled_y,
        "metadata": train_metadata_shuffled,
        "model_path": [],
        "n_estimators": [],
    },
    "test": {
        "X": test_X,
        "y": test_y,
        "metadata": test_metadata,
        "model_path": [],
        "n_estimators": [],
    },
    "test_shuffled": {
        "X": test_shuffled_X,
        "y": test_shuffled_y,
        "metadata": test_metadata_shuffled,
        "model_path": [],
        "n_estimators": [],
    },
}


# ## Model training

# In[5]:


# set the number of trees in the forest to search for the best number of trees
n_trees_list = [10, 100, 1000]


# In[6]:


# Define cross-validation strategy
cv = KFold(n_splits=5, shuffle=True, random_state=0)  # 5-fold cross-validation

for n_trees in n_trees_list:
    # Initialize the multi-output regression model
    model = MultiOutputRegressor(
        RandomForestRegressor(
            n_estimators=n_trees,
            random_state=0,
        )
    )

    # Perform cross-validation
    cv_scores = cross_val_score(
        model, train_X, train_y, cv=cv, scoring="r2"
    )  # Using R^2 as the scoring metric
    print(
        f"Mean cross-validation scores for n_trees={n_trees}:  {np.mean(cv_scores):.4f}"
    )

    # Train the model on the full training set
    model.fit(train_X, train_y)

    # Save the trained model
    model_file_path = pathlib.Path(
        f"../models/multi_regression_model_ntrees_{n_trees}.joblib"
    ).resolve()
    joblib.dump(model, model_file_path)
    dict_of_train_tests["train"]["model_path"].append(model_file_path)
    dict_of_train_tests["test"]["model_path"].append(model_file_path)
    dict_of_train_tests["train"]["n_estimators"].append(n_trees)
    dict_of_train_tests["test"]["n_estimators"].append(n_trees)

    # Repeat for shuffled data
    shuffled_model = MultiOutputRegressor(
        RandomForestRegressor(
            n_estimators=n_trees,
            random_state=0,
            n_jobs=-1,
        )
    )

    # Perform cross-validation on shuffled data
    shuffled_cv_scores = cross_val_score(
        shuffled_model, train_shuffled_X, train_shuffled_y, cv=cv, scoring="r2"
    )
    print(
        f"Mean cross-validation scores for shuffled data with n_trees={n_trees}: {np.mean(shuffled_cv_scores):.4f}"
    )

    # Train the shuffled model on the full shuffled training set
    shuffled_model.fit(train_shuffled_X, train_shuffled_y)

    # Save the shuffled model
    shuffled_model_file_path = pathlib.Path(
        f"../models/shuffled_multi_regression_model_ntrees_{n_trees}.joblib"
    ).resolve()
    joblib.dump(shuffled_model, shuffled_model_file_path)
    dict_of_train_tests["train_shuffled"]["model_path"].append(shuffled_model_file_path)
    dict_of_train_tests["test_shuffled"]["model_path"].append(shuffled_model_file_path)
    dict_of_train_tests["train_shuffled"]["n_estimators"].append(n_trees)
    dict_of_train_tests["test_shuffled"]["n_estimators"].append(n_trees)


# ## Model Evaluation

# In[7]:


output_dict_of_dfs = {"prediction_stats_df": [], "predictions_df": []}
for split in dict_of_train_tests.keys():
    if "shuffle" in split:
        shuffle = True
    else:
        shuffle = False
    X = dict_of_train_tests[split]["X"]
    y = dict_of_train_tests[split]["y"]
    metadata = dict_of_train_tests[split]["metadata"]
    for n_tree_model in enumerate(dict_of_train_tests[split]["model_path"]):
        model = joblib.load(dict_of_train_tests[split]["model_path"][n_tree_model[0]])

        n_estimators = dict_of_train_tests[split]["n_estimators"][n_tree_model[0]]

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
        prediction_stats_df["n_estimators"] = n_estimators
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
