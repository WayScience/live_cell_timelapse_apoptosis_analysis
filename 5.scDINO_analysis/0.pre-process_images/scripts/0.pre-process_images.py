#!/usr/bin/env python
# coding: utf-8

# This notebook converts single channel grayscale images to 5 channel images by adding blank channels.
# This is done to make the images compatible with the pre-trained models that expect 5 channel images.
# The code in this notebook will need to change to match a unique dataset, regretfully.
#
# Note that the data used here has four channels, but the model needs 5 channels input.

# In[1]:


import multiprocessing as mp
import pathlib
from typing import Optional

import cv2

# show the image
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tiff
import tqdm

# ## Import paths

# In[2]:


# set the path to the data directory
data_file_dir = pathlib.Path(
    "../../../4.process_features/data/annotated_data/run_20231017ChromaLive_6hr_4ch_MaxIP_sc.parquet"
).resolve(strict=True)

# read in the data
cp_feature_data = pd.read_parquet(data_file_dir)
# print the data
print(cp_feature_data.shape)
cp_feature_data.head()


# In[3]:


# get columns that contain Metadata
metadata_columns = [col for col in cp_feature_data.columns if "Metadata" in col]
metadata_df = cp_feature_data[metadata_columns]
# get columns that contain Features
feature_df = cp_feature_data.drop(columns=metadata_columns)
metadata_df.head()


# This cell is not run as it takes a long time to run...

# In[4]:


# define the function to get the metadata and crop the image


def crop_image(
    i: int,
    image_path: str,
    radius: int = 50,
    add_channels: Optional[bool] = False,
    total_channels: int = 5,
) -> None:
    """
    Crop the image based on the metadata and save the cropped image to disk
    Also output extracted metadata for the cropped image

    Parameters
    ----------
    i : int
        This is the iterator index for the metadata_df
    image_path : str
        Path to the image directory
    radius : int, optional
        The radius to crop the image by, by default 50
    add_channels : Optional[bool], optional
        This is a bool argument if set True will add extra channels to add up to 5 total , by default False

    Returns
    -------
    None
    """
    image_information_df = metadata_df.copy().iloc[i]

    center_y = image_information_df["Metadata_Nuclei_Location_Center_Y"].astype(int)
    center_x = image_information_df["Metadata_Nuclei_Location_Center_X"].astype(int)
    # DNA
    image_name_DNA = pathlib.Path(image_information_df["Metadata_Image_FileName_DNA"])
    image_path_DNA = pathlib.Path(image_path / image_name_DNA).resolve(strict=True)
    # 488_1
    image_name_488_1 = pathlib.Path(
        image_information_df["Metadata_Image_FileName_488_1"]
    )
    image_path_488_1 = pathlib.Path(image_path / image_name_488_1).resolve(strict=True)
    # 488_2
    image_name_488_2 = pathlib.Path(
        image_information_df["Metadata_Image_FileName_488_2"]
    )
    image_path_488_2 = pathlib.Path(image_path / image_name_488_2).resolve(strict=True)
    # 561
    image_name_561 = pathlib.Path(image_information_df["Metadata_Image_FileName_561"])
    image_path_561 = pathlib.Path(image_path / image_name_561).resolve(strict=True)
    image_DNA = tiff.imread(image_path_DNA)
    image_488_1 = tiff.imread(image_path_488_1)
    image_488_2 = tiff.imread(image_path_488_2)
    image_561 = tiff.imread(image_path_561)

    image_DNA_crop = image_DNA[
        center_y - radius : center_y + radius, center_x - radius : center_x + radius
    ]
    image_488_1_crop = image_488_1[
        center_y - radius : center_y + radius, center_x - radius : center_x + radius
    ]
    image_488_2_crop = image_488_2[
        center_y - radius : center_y + radius, center_x - radius : center_x + radius
    ]
    image_561_crop = image_561[
        center_y - radius : center_y + radius, center_x - radius : center_x + radius
    ]

    # check if crop is an edge case
    # Where edge case is cells that are too close to the edge of the image to crop
    # This ensures that all crops are the same dimensions and can be used in the model
    if image_DNA_crop.shape[0] < radius * 2 or image_DNA_crop.shape[1] < radius * 2:
        print(
            f"Image {image_information_df['Metadata_Image_FileName_DNA']} is an edge case. Ommitting..."
        )
        ommited_couunt = 1
        return None
    # merge the channels to a single image
    image_merge = np.stack(
        [image_DNA_crop, image_488_1_crop, image_488_2_crop, image_561_crop], axis=-1
    )
    if add_channels:
        if image_merge.shape[-1] < total_channels:
            channels_to_add = total_channels - image_merge.shape[-1]
            for _ in range(channels_to_add):
                # add a new channel of all zeros
                new_channels = np.zeros((image_merge.shape[0], image_merge.shape[1], 1))
                image_merge = np.concatenate((image_merge, new_channels), axis=-1)
    # save images to disk
    image_save_path = pathlib.Path(
        f"../../../data/processed_images/crops/{image_information_df['Metadata_Well']}/image_number_{image_information_df['Metadata_ImageNumber']}_cell_number_{image_information_df['Metadata_Nuclei_Number_Object_Number']}/"
    )
    image_save_path.mkdir(parents=True, exist_ok=True)
    file_name = image_information_df["Metadata_Image_FileName_DNA"].replace(
        ".tiff",
        f'image_number_{image_information_df["Metadata_ImageNumber"]}_{image_information_df["Metadata_Nuclei_Number_Object_Number"]}_crop.tiff',
    )
    image_save_path = pathlib.Path(image_save_path / file_name)
    tiff.imwrite(image_save_path, image_merge)
    return None


# In[5]:


image_path = pathlib.Path(
    "../../../2.cellprofiler_ic_processing/illum_directory/20231017ChromaLive_6hr_4ch_MaxIP/"
)
radius = 50


# In[6]:


# set the number of processes to use
num_processes = mp.cpu_count() - 2
print(f"Number of processes: {num_processes}")
# get the total number of rows in the metadata_df
total_crops = metadata_df.shape[0]
print(f"There are {total_crops} images to crop")

process_list = [
    mp.Process(target=crop_image, args=(i, image_path, radius, False))
    for i in range(total_crops)
]
print(f"There are {len(process_list)} processes to run")
pool = mp.Pool(num_processes)
# run the processes in the pool with multiple args
pool.starmap_async(
    crop_image, [(i, image_path, radius, False) for i in range(total_crops)]
)
pool.close()
pool.join()
pool.terminate()