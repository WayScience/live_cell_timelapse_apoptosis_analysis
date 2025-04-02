#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
import statsmodels.stats
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# In[2]:


endpoint_path = pathlib.Path(
    "../../data/CP_feature_select/endpoint_whole_image/feature_selected_whole_image.parquet"
).resolve(strict=True)
intensity_feature_path = pathlib.Path(
    "../data/0.ground_truth/annexinv_intensity_features_df.parquet"
).resolve()
intensity_feature_path.parent.mkdir(parents=True, exist_ok=True)
tukey_results_path = pathlib.Path(
    "../data/0.ground_truth/tukey_results.parquet"
).resolve()
tukey_results_path.parent.mkdir(parents=True, exist_ok=True)

endpoint_df = pd.read_parquet(endpoint_path)
endpoint_df.head()


# In[3]:


metadata_columns = [x for x in endpoint_df.columns if "Metadata_dose" in x]
# get the annexinV columns
annexinV_columns = [x for x in endpoint_df.columns if "Intensity" in x]
annexinv_df = endpoint_df[metadata_columns + annexinV_columns]

annexinv_df.head()
# save the intensity feature df

annexinv_df.to_parquet(intensity_feature_path)


# Interesting result here - should be faceted by the channel.
# I am interested in determining the key dose that is the most effective

# In[ ]:


# perform ANOVA for each intensity column for each dose
list_of_anova_results = []
for column in annexinv_df.columns:
    if column == "Metadata_dose":
        continue
    model = ols(f"{column} ~ C(Metadata_dose)", data=annexinv_df).fit()
    anova_results = anova_lm(model, typ=2)
    anova_results.reset_index(inplace=True)
    anova_results["feature"] = column
    # post hoc test
    tukey = pairwise_tukeyhsd(
        endog=annexinv_df[column], groups=annexinv_df["Metadata_dose"], alpha=0.05
    )
    tukey_results = pd.DataFrame(
        data=tukey._results_table.data[1:], columns=tukey._results_table.data[0]
    )
    tukey_results["feature"] = column
    list_of_anova_results.append(tukey_results)
df = pd.concat(list_of_anova_results)
# correct for multiple testing
df["p-adj_bh"] = statsmodels.stats.multitest.multipletests(
    df["p-adj"], method="fdr_bh"
)[1]

df.to_parquet(tukey_results_path)
df.head()
