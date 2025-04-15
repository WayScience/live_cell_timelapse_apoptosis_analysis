import numpy as np
import pandas as pd
from copairs import map
from copairs.matching import assign_reference_index


def run_mAP_across_time(
    df: pd.DataFrame,
    seed: int = 0,
    time_column: str = "Metadata_Time",
    reference_column_name: str = "Metadata_dose",
    reference_group: str = "0.0",
) -> dict[pd.DataFrame]:
    """
    Run mAP across timepoints specifies and hardcoded columns for this data

    Parameters
    ----------
    df : pd.DataFrame
        An aggregated dataframe with metadata and features and temporal information
    seed : int, optional
        Random seed for reproducibility, by default 0
    time_column : str, optional
        The column name for timepoints, by default "Metadata_Time"
    reference_column_name : str, optional
        The column name for grouping, by default "Metadata_dose"
    reference_group : str, optional
        The reference group for the analysis, by default "DMSO CTL"

    Returns
    -------
    dict
        A dictionary of dataframes with the mAP results for each
        timepoint.
    """

    unique_timepoints = df[time_column].unique()
    dict_of_map_dfs = {}
    for timepoint in unique_timepoints:
        single_time_df = df.loc[df[time_column] == timepoint]
        reference_col = "Metadata_reference_index"
        df_activity = assign_reference_index(
            single_time_df,
            f"{reference_column_name} == {reference_group}",
            reference_col=reference_col,
            default_value=-1,
        )
        pos_sameby = [reference_column_name, reference_col]
        pos_diffby = []
        neg_sameby = []
        neg_diffby = [reference_column_name, reference_col]
        metadata = df_activity.filter(regex="Metadata")
        profiles = df_activity.filter(regex="^(?!Metadata)")
        profile_shape = profiles.shape
        profiles = profiles.dropna(axis=0)
        metadata = metadata.loc[profiles.index]
        if profile_shape[0] != profiles.shape[0]:
            print(
                f"Warning: Dropped {profile_shape[0] - profiles.shape[0]} rows with NaN values"
            )
        profiles = profiles.reset_index(drop=True)
        metadata = metadata.reset_index(drop=True)
        profiles = profiles.values

        activity_ap = map.average_precision(
            metadata, profiles, pos_sameby, pos_diffby, neg_sameby, neg_diffby
        )

        activity_ap = activity_ap.query(f"{reference_column_name} != {reference_group}")
        activity_map = map.mean_average_precision(
            activity_ap, pos_sameby, null_size=1000000, threshold=0.05, seed=seed
        )
        activity_map["-log10(p-value)"] = -activity_map["corrected_p_value"].apply(
            np.log10
        )
        # flatten the multi-index columns to make it easier to work with
        dict_of_map_dfs[timepoint] = activity_map
    return dict_of_map_dfs
