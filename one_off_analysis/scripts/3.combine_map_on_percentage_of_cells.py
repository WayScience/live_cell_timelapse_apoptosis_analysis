#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pathlib

import numpy as np
import pandas as pd
from copairs import map
from copairs.matching import assign_reference_index

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False
import warnings

import pycytominer.aggregate
import tqdm

# Suppress all RuntimeWarnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

#


# In[ ]:


input_dir = pathlib.Path("../results/mAP_cell_percentages/").resolve(strict=True)
output_file = pathlib.Path("../results/mAP_cell_percentages.parquet").resolve()
# get a list of all files in the directory
df = pd.concat([pd.read_parquet(file) for file in list(input_dir.glob("*"))])
df.to_parquet(output_file)
df.head()
