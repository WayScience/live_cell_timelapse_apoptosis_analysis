#!/usr/bin/env python
# coding: utf-8

# This notebook performs a linear model regression on timelapse data to understand the temporal contribution to each morphology feature.

# Note that I also separate the disctinction between Cellprofiler (CP) and scDINO features.

# ## The linear model will be constructed as follows:
# ### $Y_{feature} = \beta_t X_t + \beta_{cell count} X_{cell count} + \beta_{Stuarosporine \space dose} X_{Stuarosporine \space dose} + \beta_{interaction_{cell \space count + time}}(X_{cell \space count}X_t) + \beta_{interaction_{dose + time}}(X_{dose}X_t) +  \beta_0$
#
#
#

# ### The model explained
# This model fits the factors to every individual morphology feature.
# This ulitmatly allows us to understand and interpret the contribution of each factor to the morphology feature.
# We also add interaction terms to further understand the relationship of multiple factors with eachother their contribution to each morphology feature.
# ### We define and interpret each term as follows:
# - $Y_{feature}$: the morphology feature we are trying to predict
#     - This is a morphology feature extracted from either CellProfiler or scDINO.
# - $\beta_t$: the coefficient for the time variable
#     - This coefficient represents the contribution of time to the morphology feature.
# - $X_t$: the time variable
#     - This variable represents the time elapsed from the start of imaging.
#     - Note that these are timepoints starting at 0 and ending at 12 with increments of 1.
# - $\beta_{cell count}$: the coefficient for the cell count variable
#     - This coefficient represents the contribution of cell count/well to the morphology feature.
# - $X_{cell count}$: the cell count variable
#     - This variable represents the number of cells in each well.
# - $\beta_{Stuarosporine \space dose}$: the coefficient for the Stuarosporine dose variable
#     - This coefficient represents the contribution of Stuarosporine dose to the morphology feature.
# - $X_{Stuarosporine \space dose}$: the Stuarosporine dose variable
#     - This variable represents the Stuarosporine dose in each well.
#     - note that the number input here is on a continuous scale with the attached units of nM.
# #### The interaction terms
# - $\beta_{interaction}(X_{cell \space count}X_t)$: the coefficient for the interaction term
#     - This coefficient represents the contribution of the interaction between cell count and time to the morphology feature.
# - $X_{cell \space count}X_t$: the interaction term of cell count and time
#     - This variable represents the interaction between cell count and time.
# - $\beta_{interaction}(X_{dose}X_t)$: the coefficient for the interaction term
#     - This coefficient represents the contribution of the interaction between Stuarosporine dose and time to the morphology feature.
# - $X_{dose}X_t$: the interaction term of Stuarosporine dose and time
#     - This variable represents the interaction between Stuarosporine dose and time.
# - $\beta_0$: the intercept
#
# #### Hypothesis:
# The null hypothesis is that all factors for every feature will have a beta coeffienct of 0. This would imply that the factors do not contribute to the morphology feature.
# The alternative hypothesis is that at least one factor for individual features will have a beta coeffienct > 0. This would imply that the factors do contribute to the morphology feature.

# In[1]:


import pathlib
import warnings

import joblib
import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm
import tqdm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


# In[2]:


def fit_linear_model(
    X: np.ndarray,
    y: np.ndarray,
    feature: str,
    write: bool = False,
) -> statsmodels.regression.linear_model.RegressionResultsWrapper:
    """
    Fit a linear model to the data and save the model to a file.

    Parameters
    ----------
    X : np.ndarray
        The input data to fit on. Should be numeric.
    y : np.ndarray
        The target data to fit on. Should be numeric.
    feature : str
        The feature name that is used from y to fit the model.
    shuffle : bool
        Whether to shuffle the data before fitting the model.

    Returns
    -------
    model : statsmodels.regression.linear_model.RegressionResultsWrapper
        The fitted model.
    """
    # Ensure X and y are numeric
    X = X.apply(pd.to_numeric, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")

    # Drop rows with missing values and match the indices in the y series
    X = X.dropna()
    y = y.loc[X.index]

    # Add a constant for the intercept
    X = sm.add_constant(X)

    # Fit the model
    model = sm.OLS(y, X).fit()
    if write:
        # write the model to a file joblib
        joblib_path = pathlib.Path(f"../linear_models/lm_{feature}.joblib").resolve()
        joblib_path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(model, joblib_path)

    return model


# In[3]:


agg_profile_file_path = pathlib.Path(
    "../../data/CP_scDINO_features/combined_CP_scDINO_norm_fs_aggregated.parquet"
).resolve(strict=True)

all_features_beta_df_path = pathlib.Path(
    "../results/all_features_beta_df.parquet"
).resolve()
all_features_beta_df_path.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_parquet(agg_profile_file_path)
df.head()


# In[4]:


# get the metadata features
metadata_columns = [x for x in df.columns if "Metadata" in x]
# get the features
feature_columns = [x for x in df.columns if "Metadata" not in x]
time_column = "Metadata_Time"
single_cells_count_column = "Metadata_number_of_singlecells"
dose_column = "Metadata_dose"
interaction_column1 = "Metadata_interaction1"
interaction_column2 = "Metadata_interaction2"

# ensure that the interaction terms are both numeric
df["Metadata_number_of_singlecells"] = pd.to_numeric(
    df["Metadata_number_of_singlecells"], errors="coerce"
)
df["Metadata_Time"] = pd.to_numeric(df["Metadata_Time"], errors="coerce")
df["Metadata_dose"] = pd.to_numeric(df["Metadata_dose"], errors="coerce")
df[interaction_column1] = df["Metadata_number_of_singlecells"] * df["Metadata_Time"]
df[interaction_column2] = df["Metadata_Time"] * df["Metadata_dose"]


# In[5]:


coefficient_names = {
    "beta": [],
    "p_value": [],
    "variate": [],
    "r2": [],
    "feature": [],
}


# In[6]:


model_dict = {}


# In[7]:


X = df[
    [
        time_column,
        single_cells_count_column,
        dose_column,
        interaction_column1,
        interaction_column2,
    ]
]
for feature in tqdm.tqdm(feature_columns):
    y = df[feature]
    model = fit_linear_model(
        X,
        y,
        feature,
    )
    # get the model coefficients and p-values
    for variate in model.params.keys():
        coefficient_names["beta"].append(model.params[variate])
        coefficient_names["variate"].append(variate)
    for pval in model.pvalues.keys():
        coefficient_names["p_value"].append(model.pvalues[pval])
        coefficient_names["r2"].append(model.rsquared)
        coefficient_names["feature"].append(feature)
    model_dict[feature] = model


# In[8]:


# write the model to a file joblib
joblib_path = pathlib.Path("../linear_models/lm_all_features.joblib").resolve()
# write the model to a file
joblib_path.parent.mkdir(parents=True, exist_ok=True)
joblib.dump(model_dict, joblib_path)


# In[9]:


all_features_beta_df = pd.DataFrame.from_dict(coefficient_names)
# remove any "Metadata_" string from the feature names
all_features_beta_df["variate"] = all_features_beta_df["variate"].str.replace(
    "Metadata_", ""
)
all_features_beta_df["variate"] = all_features_beta_df["variate"].str.replace(
    "number_of_singlecells", "Cell count"
)
all_features_beta_df.head()


# Extract feature information and save the dataframe

# In[10]:


# split the df into two dfs, one with CP and with scDINO
cp_df = all_features_beta_df[all_features_beta_df["feature"].str.contains("CP")]

cp_df = cp_df.copy()
cp_df[
    [
        "Compartment",
        "Feature_type",
        "Measurement",
        "Channel",
        "extra1",
        "extra2",
        "extra3",
        "extra4",
        "extra5",
        "extra6",
    ]
] = cp_df["feature"].str.split("_", expand=True)
cp_df["featurizer_id"] = "CP"


# realign the channel and feature types - Cellprofiler outputs an unaligned df of feature names

# make the areashape channel None
cp_df.loc[cp_df["Feature_type"].str.contains("AreaShape"), "Channel"] = "None"

cp_df.loc[cp_df["Channel"].str.contains("CL", na=False), "Channel"] = (
    cp_df["Channel"].fillna("") + "_" + cp_df["extra1"].fillna("")
)
# merge the 488 with which emmission specrtra 1 or 2
cp_df.loc[cp_df["Channel"].str.contains("488", na=False), "Channel"] = (
    cp_df["Channel"].fillna("") + "_" + cp_df["extra2"].fillna("")
)

cp_df.loc[cp_df["extra1"].str.contains("CL", na=False), "extra1"] = (
    cp_df["extra1"].fillna("") + "_" + cp_df["extra2"].fillna("")
)
cp_df.loc[cp_df["extra1"].str.contains("488", na=False), "extra1"] = (
    cp_df["extra1"].fillna("") + "_" + cp_df["extra3"].fillna("")
)

# specify that it is ChomaLive (CL)
cp_df.loc[cp_df["extra2"].str.contains("CL", na=False), "extra2"] = (
    cp_df["extra2"].fillna("") + "_" + cp_df["extra3"].fillna("")
)
cp_df.loc[cp_df["extra2"].str.contains("488", na=False), "extra2"] = (
    cp_df["extra2"].fillna("") + "_" + cp_df["extra4"].fillna("")
)

cp_df.loc[cp_df["extra3"].str.contains("CL", na=False), "extra3"] = (
    cp_df["extra3"].fillna("") + "_" + cp_df["extra4"].fillna("")
)
cp_df.loc[cp_df["extra3"].str.contains("488", na=False), "extra3"] = (
    cp_df["extra3"].fillna("") + "_" + cp_df["extra5"].fillna("")
)

cp_df.loc[cp_df["extra4"].str.contains("CL", na=False), "extra4"] = (
    cp_df["extra4"].fillna("") + "_" + cp_df["extra5"].fillna("")
)
cp_df.loc[cp_df["extra4"].str.contains("488", na=False), "extra4"] = (
    cp_df["extra4"].fillna("") + "_" + cp_df["extra6"].fillna("")
)

cp_df.rename(columns={"extra1": "Channel2"}, inplace=True)
# remove the extra columns to retain the feature types and channels
cp_df.drop(
    columns=["Measurement", "extra2", "extra3", "extra4", "extra5", "extra6"],
    inplace=True,
)
# make channel2 None if feature is not correlation
cp_df.loc[~cp_df["Feature_type"].str.contains("Correlation"), "Channel2"] = "None"


# In[11]:


scdino_df = all_features_beta_df[all_features_beta_df["feature"].str.contains("scDINO")]
scdino_df = scdino_df.copy()
scdino_df["feature"] = scdino_df["feature"].str.replace("channel_", "")
scdino_df["feature"] = scdino_df["feature"].str.replace("channel", "")
scdino_df[["Channel", "remove", "feature", "feature_number", "featurizer_id"]] = (
    scdino_df["feature"].str.split("_", expand=True)
)
scdino_df.drop(columns=["remove"], inplace=True)
# set scDINO to be the feature type, Compartment and measurement
# as scDINO does not have these features only a channel and feature number
scdino_df[["Compartment", "Feature_type", "Measurement"]] = "scDINO"


# In[12]:


final_df = pd.concat([cp_df, scdino_df], axis=0)
final_df["Channel"] = final_df["Channel"].str.replace("Adjacent", "None")

final_df["Channel"] = final_df["Channel"].str.replace("CL_488_1", "488-1")
final_df["Channel"] = final_df["Channel"].str.replace("CL_488_2", "488-2")
final_df["Channel"] = final_df["Channel"].str.replace("CL_561", "561")
final_df["Channel"] = final_df["Channel"].str.replace("488-1", "CL 488-1")
final_df["Channel"] = final_df["Channel"].str.replace("488-2", "CL 488-2")
final_df["Channel"] = final_df["Channel"].str.replace("561", "CL 561")
final_df["Channel"].unique()


# In[13]:


final_df["variate"] = final_df["variate"].str.replace("dose", "Dose")
final_df["variate"] = final_df["variate"].str.replace(
    "interaction1", "Time x \nCell count"
)
final_df["variate"] = final_df["variate"].str.replace("interaction2", "Time x \nDose")


# In[14]:


# correct the p-values using the Benjamini/Hochberg method
final_df["p_value_corrected"] = multipletests(final_df["p_value"], method="fdr_bh")[1]


# In[15]:


# save the final df to a file
final_df.to_parquet(all_features_beta_df_path, index=False)
final_df.head()
