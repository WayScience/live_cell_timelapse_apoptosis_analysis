#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import warnings
from typing import List, Tuple

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tqdm
from sklearn.decomposition import PCA
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNetCV
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
terminal_columns = [
    x for x in train_df.columns if "Terminal" in x and "Metadata" not in x
]


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

feature_columns = train_y.columns.tolist()


# In[4]:


single_feature = "Terminal_Cytoplasm_Intensity_UpperQuartileIntensity_AnnexinV"


# In[5]:


# Define cross-validation strategy
cv = KFold(n_splits=5, shuffle=True, random_state=0)  # 5-fold cross-validation
# elastic net parameters
elastic_net_params = {
    "alpha": [0.1, 1.0, 10.0, 100.0, 1000.0],  # Regularization strength
    "l1_ratio": [0.1, 0.25, 0.5, 0.75, 1.0],  # l1_ratio = 1.0 is Lasso
    "max_iter": 10000,  # Increase max_iter for convergence
}
elastic_net_all_annexinv_features_model = MultiOutputRegressor(
    ElasticNetCV(
        alphas=elastic_net_params["alpha"],
        l1_ratio=elastic_net_params["l1_ratio"],
        cv=cv,
        random_state=0,
        max_iter=elastic_net_params["max_iter"],
    )
)
elastic_net_all_annexinv_features_model_shuffled = (
    elastic_net_all_annexinv_features_model
)
elastic_net_single_terminal_features_model = ElasticNetCV(
    alphas=elastic_net_params["alpha"],
    l1_ratio=elastic_net_params["l1_ratio"],
    cv=cv,
    random_state=0,
    max_iter=elastic_net_params["max_iter"],
)
elastic_net_single_terminal_features_model_shuffled = (
    elastic_net_single_terminal_features_model
)


# In[6]:


dict_of_train_tests = {
    "single_feature": {
        "train": {
            "X": train_X,
            "y": train_y[single_feature],
            "metadata": train_metadata,
            "model": elastic_net_single_terminal_features_model,
            "model_name": "elastic_net_single_terminal_features_model",
        },
        "train_shuffled": {
            "X": train_shuffled_X,
            "y": train_shuffled_y[single_feature],
            "metadata": train_metadata_shuffled,
            "model": elastic_net_single_terminal_features_model_shuffled,
            "model_name": "elastic_net_single_terminal_features_model_shuffled",
        },
        "test": {
            "X": test_X,
            "y": test_y[single_feature],
            "metadata": test_metadata,
            "model": elastic_net_single_terminal_features_model,
            "model_name": "elastic_net_single_terminal_features_model",
        },
        "test_shuffled": {
            "X": test_shuffled_X,
            "y": test_shuffled_y[single_feature],
            "metadata": test_metadata_shuffled,
            "model": elastic_net_single_terminal_features_model_shuffled,
            "model_name": "elastic_net_single_terminal_features_model_shuffled",
        },
    },
    "annexinV_features": {
        "train": {
            "X": train_X,
            "y": train_y,
            "metadata": train_metadata,
            "model": elastic_net_all_annexinv_features_model,
            "model_name": "elastic_net_all_annexinv_features_model",
        },
        "train_shuffled": {
            "X": train_shuffled_X,
            "y": train_shuffled_y,
            "metadata": train_metadata_shuffled,
            "model": elastic_net_all_annexinv_features_model_shuffled,
            "model_name": "elastic_net_all_annexinv_features_model_shuffled",
        },
        "test": {
            "X": test_X,
            "y": test_y,
            "metadata": test_metadata,
            "model": elastic_net_all_annexinv_features_model,
            "model_name": "elastic_net_all_annexinv_features_model",
        },
        "test_shuffled": {
            "X": test_shuffled_X,
            "y": test_shuffled_y,
            "metadata": test_metadata_shuffled,
            "model": elastic_net_all_annexinv_features_model_shuffled,
            "model_name": "elastic_net_all_annexinv_features_model_shuffled",
        },
    },
}


# ## Model training

# In[7]:


# Define cross-validation strategy
cv = KFold(n_splits=5, shuffle=True, random_state=0)  # 5-fold cross-validation
# elastic net parameters
elastic_net_params = {
    "alpha": [0.1, 1.0, 10.0, 100.0, 1000.0],  # Regularization strength
    "l1_ratio": [0.1, 0.25, 0.5, 0.75, 1.0],  # l1_ratio = 1.0 is Lasso
    "max_iter": 10000,  # Increase max_iter for convergence
}
elastic_net_all_terminal_features_model = MultiOutputRegressor(
    ElasticNetCV(
        alphas=elastic_net_params["alpha"],
        l1_ratio=elastic_net_params["l1_ratio"],
        cv=cv,
        random_state=0,
        max_iter=elastic_net_params["max_iter"],
    )
)

elastic_net_single_terminal_features_model = ElasticNetCV(
    alphas=elastic_net_params["alpha"],
    l1_ratio=elastic_net_params["l1_ratio"],
    cv=cv,
    random_state=0,
    max_iter=elastic_net_params["max_iter"],
)

# train the model
for model_type in dict_of_train_tests.keys():
    for train_test_key, train_test_data in tqdm.tqdm(
        dict_of_train_tests[model_type].items()
    ):
        if "test" in train_test_key:
            print(f"Skipping {train_test_key} as it is a test set.")
            continue
        print(f"Training model for {train_test_key}...{model_type}")
        X = train_test_data["X"]
        y = train_test_data["y"]
        metadata = train_test_data["metadata"]
        print(
            f"X shape: {X.shape}, y shape: {y.shape}, metadata shape: {metadata.shape}"
        )
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ConvergenceWarning)
            train_test_data["model"].fit(X, y)

        # save the model
        model_path = (
            model_dir / f"{train_test_key}_{train_test_data['model_name']}.joblib"
        )
        joblib.dump(train_test_data["model"], model_path)
        dict_of_train_tests[model_type][train_test_key]["model_path"] = model_path


# In[8]:


# test the model
for model_type in dict_of_train_tests.keys():
    for train_test_key, train_test_data in tqdm.tqdm(
        dict_of_train_tests[model_type].items()
    ):
        if "train" in train_test_key:
            print(f"Skipping {train_test_key} as it is a training set.")
            continue
        print(model_type, train_test_key)
        X = train_test_data["X"]
        y = train_test_data["y"]
        metadata = train_test_data["metadata"]
        if "shuffled" in train_test_key:
            model_path = dict_of_train_tests[model_type]["train_shuffled"]["model_path"]
        else:
            model_path = dict_of_train_tests[model_type]["train"]["model_path"]

        # load the model
        model = joblib.load(model_path)

        # make predictions
        y_pred = model.predict(X)
        if model_type == "single_feature":
            model.alpha_
            model.l1_ratio_
        else:

            alphas = model.estimators_[0].alpha_
            l1_ratios = model.estimators_[0].l1_ratio_
            print(f"Model parameters for {train_test_key}:")
            print(f"Alphas: {alphas}, L1 Ratios: {l1_ratios}")

        # calculate metrics
        metrics = {
            "explained_variance": explained_variance_score(y, y_pred),
            "mean_absolute_error": mean_absolute_error(y, y_pred),
            "mean_squared_error": mean_squared_error(y, y_pred),
            "r2_score": r2_score(y, y_pred),
        }


# In[9]:


# write the feature columns to a file
with open("../results/terminal_columns.txt", "w") as f:
    for col in feature_columns:
        f.write(f"{col}\n")
