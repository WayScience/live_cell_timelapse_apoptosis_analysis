#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import random
import warnings

import numpy as np
import pandas as pd
from copairs import map
from copairs.matching import assign_reference_index

# suppress warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# In[2]:


well_to_dose_path = pathlib.Path(
    "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs.parquet"
).resolve(strict=True)
well_to_dose_df = pd.read_parquet(
    well_to_dose_path, columns=["Metadata_Well", "Metadata_dose"]
)
# drop duplicates
well_to_dose_df = well_to_dose_df.drop_duplicates(
    subset=["Metadata_Well", "Metadata_dose"]
)
# make dict that maps well to dose
well_to_dose = well_to_dose_df.set_index("Metadata_Well")["Metadata_dose"].to_dict()


# In[3]:


CP_scDINO_profile_file_path = pathlib.Path(
    "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs_aggregated.parquet"
).resolve(strict=True)
df = pd.read_parquet(CP_scDINO_profile_file_path)
# map dose to well from well_to_dose dict
df["dose"] = df["Metadata_Well"].map(well_to_dose)
dose = df.pop("dose")
df.insert(0, "Metadata_dose", dose)
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


df["Metadata_Time"] = df["Metadata_Time"].astype(float)
timepoints = df["Metadata_Time"].unique().tolist()
timepoints.sort()
shuffle_options = [False, True]


# In[6]:


# This loops through the CP, scDINO, and CP_scDINO feature sets and calculates the mean average precision for each timepoint
for shuffle in shuffle_options:

    for feature_set_name, feature_set in feature_set_dict.items():
        feature_subset_df = pd.concat([df[feature_set], df[metadata_columns]], axis=1)

        list_of_maps = []
        # calculate the mean average precision for each timepoint
        # and each shuffle option
        for timepoint in timepoints:
            if shuffle:
                random.seed(0)
                # permutate the data
                for col in feature_subset_df.columns:
                    if col != "Metadata_dose":
                        feature_subset_df[col] = np.random.permutation(
                            feature_subset_df[col]
                        )

            time_subset_df = feature_subset_df[
                feature_subset_df["Metadata_Time"] == timepoint
            ]
            reference_col = "Metadata_reference_index"
            df_activity = assign_reference_index(
                time_subset_df,
                "Metadata_dose == '0.0'",
                reference_col=reference_col,
                default_value=-1,
            )
            pos_sameby = ["Metadata_dose", reference_col]
            pos_diffby = []
            neg_sameby = []
            neg_diffby = ["Metadata_dose", reference_col]
            metadata = df_activity.filter(regex="Metadata")
            profiles = df_activity.filter(regex="^(?!Metadata)")
            # drop the nans from the profiles
            profiles = profiles.dropna(axis=0)
            metadata = metadata.loc[profiles.index]
            profiles = profiles.values
            # suppress the warning
            activity_ap = map.average_precision(
                metadata, profiles, pos_sameby, pos_diffby, neg_sameby, neg_diffby
            )

            activity_ap = activity_ap.query("Metadata_dose != '0.0'")
            activity_map = map.mean_average_precision(
                activity_ap, pos_sameby, null_size=1000000, threshold=0.05, seed=0
            )
            activity_map.reset_index(drop=True, inplace=True)
            activity_map.insert(
                0,
                "Metadata_Time",
                timepoint,
            )
            if shuffle:
                activity_map.insert(
                    1,
                    "Shuffle",
                    "True",
                )
            else:
                activity_map.insert(
                    1,
                    "Shuffle",
                    "False",
                )
            list_of_maps.append(activity_map)

    final_mAP_df = pd.concat(list_of_maps, ignore_index=True)
    # save the final dataframe to a parquet file
    final_mAP_path = pathlib.Path(
        f"../data/mAP/mAP_scores_{feature_set_name}.parquet"
    ).resolve()
    final_mAP_path.parent.mkdir(parents=True, exist_ok=True)
    final_mAP_df.to_parquet(
        final_mAP_path,
        index=False,
    )
