#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import random

import numpy as np

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False
import itertools

#


# In[2]:


seeds_path = pathlib.Path("../combinations/seeds.txt").resolve()
percentage_path = pathlib.Path("../combinations/percentage.txt").resolve()

seeds_path.parent.mkdir(parents=True, exist_ok=True)


# In[3]:


random.seed(0)
num_of_iterations = 100
seed_space = np.random.randint(0, 1000000, num_of_iterations).tolist()  # min, max
percentages_to_run = np.linspace(0.1, 1, 10).tolist()  # start, stop, num
percentages_to_run = [round(x, 2) for x in percentages_to_run]
# write each combination to a line in the file
with open(seeds_path, "w") as f:
    for seed in seed_space:
        f.write(f"{seed}\n")

with open(percentage_path, "w") as f:
    for percentage in percentages_to_run:
        f.write(f"{percentage}\n")
