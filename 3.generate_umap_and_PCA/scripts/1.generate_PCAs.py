#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
from sklearn.decomposition import PCA

# In[2]:


def fit_pca_to_the_first_timepoint(
    df: pd.DataFrame,
    timepoint_column: str,
    metadata_columns: list,
    feature_columns: list,
    pca_model: PCA,
) -> pd.DataFrame:
    """
    This function fits a pca model to the first timepoint of the data and then applies the model to the rest of the data.

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
    pca_model : pca.pca
        The pca model to use

    Returns
    -------
    pd.DataFrame
        The pca embeddings for the data, with the metadata columns included.
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
    _ = pca_model.fit_transform(first_timepoint_subset_df)

    # get the rest of the data fo transformation
    df = df.drop(metadata_columns, axis=1)
    df = df[feature_columns]
    df.dropna(axis=0, inplace=True)
    metadata_df = metadata_df.loc[df.index]
    df.reset_index(drop=True, inplace=True)
    metadata_df.reset_index(drop=True, inplace=True)

    # apply the model to the rest of the data
    pca_embeddings = pca_model.transform(df)
    # create a dataframe with the pca fit and the metadata
    pca_df = pd.DataFrame(pca_embeddings, columns=["PCA0", "PCA1"])
    # add the metadata to the dataframe
    pca_df = pd.concat([pca_df, metadata_df], axis=1)

    return pca_df


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


pca_model = PCA(n_components=2)


# In[5]:


for profile_level in dictionary_of_feature_sets.keys():
    for profile in dictionary_of_feature_sets[profile_level].keys():
        profile_df = pd.read_parquet(dictionary_of_feature_sets[profile_level][profile])
        metadata_columns = [x for x in profile_df.columns if "Metadata_" in x]
        feature_columns = [x for x in profile_df.columns if "Metadata_" not in x]
        pca_df = fit_pca_to_the_first_timepoint(
            profile_df,
            timepoint_column="Metadata_Time",
            metadata_columns=metadata_columns,
            feature_columns=feature_columns,
            pca_model=pca_model,
        )
        # set the save path of the pca data
        pca_save_path = pathlib.Path(
            f"../results/pca/{profile_level}_{profile}_pca.parquet"
        ).resolve()
        pca_save_path.parent.mkdir(parents=True, exist_ok=True)
        # save the pca data
        pca_df.to_parquet(pca_save_path, index=False)
