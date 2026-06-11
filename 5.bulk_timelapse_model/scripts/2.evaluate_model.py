#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

# In[2]:


def model_stats_grab(
    predicted_df: pd.DataFrame,
    actual_df: pd.DataFrame,
) -> tuple:
    """
    Calculate model statistics such as Mean Squared Error (MSE), Mean Absolute Error (MAE), and R-squared (R2) from predicted and actual values.

        Parameters
        ----------
        predicted_df : pd.DataFrame
            DataFrame containing the predicted values.
        actual_df : pd.DataFrame
            DataFrame containing the actual values.

        Returns
        -------
        tuple
        mse : float
            Mean Squared Error between predicted and actual values.
        mae : float
            Mean Absolute Error between predicted and actual values.
        r2 : float
            R-squared value indicating the proportion of variance explained by the model.
    """

    assert (
        predicted_df.shape == actual_df.shape
    ), "Predicted and actual DataFrames must have the same shape."
    # if predicted_df.isinstance(pd.DataFrame):
    #     predicted = predicted_df.values
    # if isinstance(predicted_df, pd.Series):
    #     actual = actual_df.values

    mse = mean_squared_error(actual_df, predicted_df)
    mae = mean_absolute_error(actual_df, predicted_df)
    r2 = r2_score(actual_df, predicted_df)

    return mse, mae, r2


# In[3]:


# load the training data
profile_file_dir = pathlib.Path("../data_splits/test.parquet").resolve(strict=True)

models_path = pathlib.Path("../models").resolve(strict=True)

terminal_column_names = pathlib.Path("../results/terminal_columns.txt").resolve(
    strict=True
)
predictions_save_path = pathlib.Path("../results/model_performances.parquet").resolve()
terminal_column_names = [
    line.strip() for line in terminal_column_names.read_text().splitlines()
]
results_dir = pathlib.Path("../results/").resolve()
results_dir.mkdir(parents=True, exist_ok=True)

profile_df = pd.read_parquet(profile_file_dir)
print(profile_df.shape)
profile_df.head()


# In[4]:


single_feature = "Terminal_Cytoplasm_Intensity_UpperQuartileIntensity_AnnexinV"


# In[5]:


terminal_df = profile_df[terminal_column_names]
single_feature_df = profile_df[[single_feature]]
profile_df = profile_df.drop(columns=terminal_column_names)


# In[ ]:


models = pathlib.Path(models_path).glob("*.joblib")
models_dict = {
    "model_name": [],
    "model_path": [],
    "shuffled": [],
    "feature": [],
}

for model_path in models:
    models_dict["model_name"].append(model_path.name)
    models_dict["model_path"].append(model_path)
    models_dict["shuffled"].append(
        "shuffled" if "shuffled" in model_path.name else "not_shuffled"
    )
    models_dict["feature"].append(
        "Terminal_Cytoplasm_Intensity_UpperQuartileIntensity_AnnexinV"
        if "terminal_feature" in model_path.name
        else "all_terminal_features"
    )


# In[7]:


results_dict = {
    "model_name": [],
    "shuffled": [],
    "feature": [],
    "mse": [],
    "mae": [],
    "r2": [],
}


# In[8]:


metadata_columns = [x for x in profile_df.columns if "metadata" in x.lower()]
terminal_column_names = [x for x in profile_df.columns if "Terminal" in x]
features_df = profile_df.drop(columns=metadata_columns, errors="ignore")


# In[9]:


models_dict


# In[10]:


for i, model_name in enumerate(models_dict["model_name"]):
    print(f"Processing model {i + 1}/{len(models_dict['model_name'])}: {model_name}")
    print(models_dict["feature"][i])
    model = joblib.load(models_dict["model_path"][i])
    if "all_terminal_features" not in str(models_dict["feature"][i]):
        predictions = model.predict(features_df[model.feature_names_in_])
        print(predictions.shape, single_feature_df[single_feature].shape)

        mse, mae, r2 = model_stats_grab(predictions, single_feature_df[single_feature])
        results_dict["shuffled"].append(models_dict["shuffled"][i])
        results_dict["feature"].append(models_dict["feature"][i])
    else:
        predictions = model.predict(features_df[model.feature_names_in_])
        print(predictions.shape, terminal_df.shape)
        mse, mae, r2 = model_stats_grab(predictions, terminal_df)
        results_dict["shuffled"].append(models_dict["shuffled"][i])
        results_dict["feature"].append("all_terminal_features")

    results_dict["model_name"].append(model_name)
    results_dict["mse"].append(mse)
    results_dict["mae"].append(mae)
    results_dict["r2"].append(r2)
results_df = pd.DataFrame(results_dict)
results_df.to_parquet(predictions_save_path, index=False)
results_df.head()


# In[11]:


# plot the performance of the models

sns.set(style="whitegrid")
plt.figure(figsize=(12, 6))
sns.barplot(x="feature", y="mse", hue="shuffled", data=results_df, palette="viridis")
plt.xticks(rotation=45, ha="right")
plt.title("Model Performance: Mean Squared Error (MSE)")
plt.tight_layout()
plt.figure(figsize=(12, 6))
sns.barplot(x="feature", y="mae", hue="shuffled", data=results_df, palette="viridis")
plt.xticks(rotation=45, ha="right")
plt.title("Model Performance: Mean Absolute Error (MAE)")
plt.tight_layout()
plt.figure(figsize=(12, 6))
sns.barplot(x="feature", y="r2", hue="shuffled", data=results_df, palette="viridis")
plt.xticks(rotation=45, ha="right")
plt.title("Model Performance: R-squared (R2)")
plt.tight_layout()
