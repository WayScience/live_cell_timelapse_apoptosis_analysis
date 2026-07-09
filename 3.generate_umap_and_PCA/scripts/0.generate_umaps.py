#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
import umap

# In[2]:


def fit_umap_to_the_first_timepoint(
    df: pd.DataFrame,
    timepoint_column: str,
    metadata_columns: list,
    feature_columns: list,
    umap_model: umap.UMAP,
) -> pd.DataFrame:
    """
    This function fits a UMAP model to the first timepoint of the data and then applies the model to the rest of the data.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing all feature, metadata, and timepoint columns.
    timepoint_column : str
        The name of the column containing the timepoint information
    metadata_columns : list
        The names of the columns containing the metadata information
    feature_columns : list
        The names of the columns containing the feature information
    umap_model : umap.UMAP
        The UMAP model to use.

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


def fit_umap_to_all_timepoints(
    df: pd.DataFrame,
    metadata_columns: list,
    feature_columns: list,
    umap_model: umap.UMAP,
) -> pd.DataFrame:
    """
    This function fits a UMAP model to all timepoints of the data.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing all feature, metadata, and timepoint columns.
    timepoint_column : str
        The name of the column containing the timepoint information
    metadata_columns : list
        The names of the columns containing the metadata information
    feature_columns : list
        The names of the columns containing the feature information
    umap_model : umap.UMAP
        The UMAP model to use.

    Returns
    -------
    pd.DataFrame
        The UMAP embeddings for the data, with the metadata columns included.
    """

    df = df.copy()
    metadata_df = df[metadata_columns]

    # Prepare the data by dropping metadata columns, selecting feature columns, and removing rows with missing values
    df = df.drop(metadata_columns, axis=1)
    df = df[feature_columns]
    df.dropna(axis=0, inplace=True)
    metadata_df = metadata_df.loc[df.index]
    df.reset_index(drop=True, inplace=True)
    metadata_df.reset_index(drop=True, inplace=True)

    # fit the model to all timepoints
    umap_embeddings = umap_model.fit_transform(df)
    # create a dataframe with the umap fit and the metadata
    umap_df = pd.DataFrame(umap_embeddings, columns=["UMAP_0", "UMAP_1"])
    # add the metadata to the dataframe
    umap_df = pd.concat([umap_df, metadata_df], axis=1)

    return umap_df


# In[3]:


dictionary_of_feature_sets = {
    "single-cell_profiles": {
        "CP": pathlib.Path(
            "../../data/CP_feature_select/profiles/features_selected_profile.parquet"
        ).resolve(strict=True),
        "scDINO": pathlib.Path(
            "../../data/scDINO/CLS_features_annotated_normalized_feature_selected.parquet"
        ).resolve(strict=True),
        "CP_scDINO": pathlib.Path(
            "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs.parquet"
        ).resolve(strict=True),
    },
    "bulk_profiles": {
        "CP": pathlib.Path(
            "../../data/CP_aggregated/profiles/aggregated_profile.parquet"
        ).resolve(strict=True),
        "scDINO": pathlib.Path(
            "../../data/scDINO/CLS_features_annotated_normalized_feature_selected_aggregated.parquet"
        ).resolve(strict=True),
        "CP_scDINO": pathlib.Path(
            "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs_aggregated.parquet"
        ).resolve(strict=True),
    },
}


# In[4]:


umap_model = umap.UMAP(
    n_neighbors=15,  # higher number focuses on global structure
    n_components=2,
    metric="euclidean",
    random_state=0,
    min_dist=0.1,  # lower number means tighter points
    spread=0.5,  #
)


# In[5]:


for profile_level in dictionary_of_feature_sets.keys():
    for profile in dictionary_of_feature_sets[profile_level].keys():
        profile_df = pd.read_parquet(dictionary_of_feature_sets[profile_level][profile])
        metadata_columns = [x for x in profile_df.columns if "Metadata_" in x]
        feature_columns = [x for x in profile_df.columns if "Metadata_" not in x]
        umap_df = fit_umap_to_the_first_timepoint(
            profile_df,
            timepoint_column="Metadata_Time",
            metadata_columns=metadata_columns,
            feature_columns=feature_columns,
            umap_model=umap_model,
        )
        # set the save path of the umap data
        umap_save_path = pathlib.Path(
            f"../results/UMAP/{profile_level}_{profile}_umap.parquet"
        ).resolve()
        umap_save_path.parent.mkdir(parents=True, exist_ok=True)
        # save the umap data
        umap_df.to_parquet(umap_save_path, index=False)
        all_timepoints_umap_df = fit_umap_to_all_timepoints(
            profile_df,
            metadata_columns=metadata_columns,
            feature_columns=feature_columns,
            umap_model=umap_model,
        )
        # set the save path of the umap data for all timepoints
        all_timepoints_umap_save_path = pathlib.Path(
            f"../results/UMAP/{profile_level}_{profile}_all_timepoints_umap.parquet"
        ).resolve()
        all_timepoints_umap_save_path.parent.mkdir(parents=True, exist_ok=True)
        # save the umap data for all timepoints
        all_timepoints_umap_df.to_parquet(all_timepoints_umap_save_path, index=False)
