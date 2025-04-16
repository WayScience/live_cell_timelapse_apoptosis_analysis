#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tqdm
import umap

# In[2]:


def fit_umap_to_the_first_timepoint(
    df: pd.DataFrame,
    timepoint_column: str = "Metadata_Time",
    metadata_columns: list = None,
    feature_columns: list = None,
    umap_model: umap.UMAP = None,
) -> pd.DataFrame:
    """
    This function fits a UMAP model to the first timepoint of the data and then applies the model to the rest of the data.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing all feature, metadata, and timepoint columns.
    timepoint_column : str, optional
        The name of the column containing the timepoint information, by default "Metadata_Time"
    metadata_columns : list, optional
        The names of the columns containing the metadata information, by default None
    feature_columns : list, optional
        The names of the columns containing the feature information, by default None
    umap_model : umap.UMAP, optional
        The UMAP model to use, by default None. If None, a new UMAP model will be created with default parameters.

    Returns
    -------
    pd.DataFrame
        The UMAP embeddings for the data, with the metadata columns included.
    """

    df = df.copy()
    metadata_df = df[metadata_columns]

    # get the first timepoint and the subset of the data for that timepoint
    first_time = df[timepoint_column].min()
    first_timepoint_subset_df = df[df[timepoint_column] == first_time]

    # Prepare the first timepoint subset by dropping metadata columns, selecting feature columns, and removing rows with missing values
    first_timepoint_subset_df = first_timepoint_subset_df.drop(metadata_columns, axis=1)
    first_timepoint_subset_df = first_timepoint_subset_df[feature_columns]
    first_timepoint_subset_df = first_timepoint_subset_df.dropna(axis=0)
    # fit the model to the first timepoint
    _ = umap_model.fit_transform(first_timepoint_subset_df)

    # get the rest of the data fo transformation
    df = df.drop(metadata_columns, axis=1)
    df = df[feature_columns]
    df.dropna(axis=0, inplace=True)
    metadata_df = metadata_df.loc[df.index]
    df.reset_index(drop=True, inplace=True)
    metadata_df.reset_index(drop=True, inplace=True)

    # apply the model to the rest of the data
    umap_embeddings = umap_model.transform(df)
    # create a dataframe with the umap fit and the metadata
    umap_df = pd.DataFrame(umap_embeddings, columns=["UMAP_0", "UMAP_1"])
    # add the metadata to the dataframe
    umap_df = pd.concat([umap_df, metadata_df], axis=1)

    return umap_df


# In[3]:


CP_scDINO_profile_file_path = pathlib.Path(
    "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs.parquet"
).resolve(strict=True)
df = pd.read_parquet(CP_scDINO_profile_file_path)
df.head()


# In[4]:


metadata_columns = [x for x in df.columns if "Metadata" in x]
scDINO_columns = [x for x in df.columns if "scDINO" in x]
CP_columns = df.drop(columns=metadata_columns + scDINO_columns).columns
CP_scDINO_columns = df.drop(metadata_columns, axis=1).columns

feature_set_dict = {
    "scDINO": scDINO_columns,
    "CP": CP_columns,
    "CP_scDINO": CP_scDINO_columns,
}


# In[5]:


umap_model = umap.UMAP(
    n_neighbors=15,  # higher number focuses on global structure
    n_components=2,
    metric="euclidean",
    random_state=0,
    min_dist=0.1,  # lower number means tighter points
    spread=0.5,  #
)


# In[6]:


for feature_set_name, feature_set in tqdm.tqdm(feature_set_dict.items()):
    umap_df = fit_umap_to_the_first_timepoint(
        df,
        timepoint_column="Metadata_Time",
        metadata_columns=metadata_columns,
        feature_columns=feature_set,
        umap_model=umap_model,
    )
    # set the save path of the umap data
    umap_save_path = pathlib.Path(
        f"../results/UMAP/{feature_set_name}_umap.parquet"
    ).resolve()
    umap_save_path.parent.mkdir(parents=True, exist_ok=True)
    # save the umap data
    umap_df.to_parquet(umap_save_path, index=False)
