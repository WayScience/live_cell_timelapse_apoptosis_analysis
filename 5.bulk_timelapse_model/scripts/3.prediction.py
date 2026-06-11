#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import joblib
import numpy as np
import pandas as pd

# In[2]:


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
models_path = pathlib.Path("../models").resolve()
data_split_df = pd.read_parquet(train_test_wells_path)
df = pd.read_parquet(profile_data_path)


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


# In[4]:


# map the train/test wells to the aggregate data
df["Metadata_data_split"] = df["Metadata_Well"].map(
    data_split_df.set_index("Metadata_Well")["data_split"]
)
data_split = df.pop("Metadata_data_split")
df.insert(0, "Metadata_data_split", data_split)
df["Metadata_Time"] = df["Metadata_Time"].astype(float)
# drop NaN values in the terminal columns
df = df.dropna(subset="Metadata_data_split")
df["Metadata_data_split"].unique()


# In[5]:


# if the data_split is train and the time is not 12 then set to non_trained_pair
# where 12 is the last time point
df["Metadata_data_split"] = df.apply(
    lambda x: (
        "non_trained_pair"
        if (x["Metadata_data_split"] == "train" and x["Metadata_Time"] != 12.0)
        else x["Metadata_data_split"]
    ),
    axis=1,
)


# In[6]:


metadata_columns = [x for x in df.columns if "metadata" in x.lower()]
aggregate_features_df = df.drop(columns=metadata_columns, errors="ignore")


# In[7]:


results_dict = {}
for i, model_name in enumerate(models_dict["feature"]):
    model = joblib.load(models_dict["model_path"][i])
    if models_dict["feature"][i] != "all_terminal_features":
        print(models_dict["feature"][i])
        predicted_df = pd.DataFrame(
            model.predict(aggregate_features_df),
            columns=[models_dict["feature"][i]],
        )
    else:
        print("all_terminal_features")
        predicted_df = pd.DataFrame(
            model.predict(aggregate_features_df),
            columns=terminal_column_names,
        )
    predicted_df[metadata_columns] = df[metadata_columns]
    predicted_df["shuffled"] = models_dict["shuffled"][i]
    # drop nan value
    predicted_df = predicted_df.dropna()

    # check if a key for the feature already exists in results_dict
    if f"{models_dict['feature'][i]}" in results_dict:
        temporary_df = pd.concat(
            [results_dict[f"{models_dict['feature'][i]}"], predicted_df],
            ignore_index=True,
            sort=False,
        )
        results_dict[f"{models_dict['feature'][i]}"] = temporary_df
    else:
        results_dict[f"{models_dict['feature'][i]}"] = predicted_df

    print(results_dict[f"{models_dict['feature'][i]}"].shape)


# In[8]:


for model in results_dict.keys():
    save_path = pathlib.Path(f"../results/{model}.parquet").resolve()
    save_path.parent.mkdir(parents=True, exist_ok=True)
    results_dict[model].to_parquet(save_path, index=False)
