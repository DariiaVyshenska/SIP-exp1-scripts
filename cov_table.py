#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 09:13:50 2019

@author: DVyshenska
"""

import pandas as pd
#import numpy as np
from os import chdir, getcwd, listdir
wd=getcwd()
chdir(wd)

#%%
def dispatch1 (vec, action):
    if action == "mean":
        return vec.mean()
    elif action == "median":
        return vec.median()

def files_list_f(pat_match, folder_path):
    return [f for f in listdir(folder_path) if pat_match in f]

def file_id_f(f_name, pat_name, pat_ext):
    return f_name.replace(pat_name, "").replace(pat_ext, "")

def import_data(sum_table, basecov_path, microbe, stat):
    col_name = microbe + "_basecov"
    sum_table[col_name] = "NA"
    for file_name in files_list_f("basecov_", basecov_path):
        file_id = file_id_f(file_name, "basecov_", ".txt")
        print (file_name)
        t = pd.read_csv(basecov_path + file_name, sep="\t")
        sum_table.loc[file_id,col_name] = dispatch1(t.Coverage, stat)
    print("DONE!")

#%%

basecov_path = "./map_basecov_files/ecoli/"
stat = "mean"
microbe = "ecoli"

#%%
sum_table = pd.read_csv("./summary_files/summary_out_wd.csv")
sum_table.set_index('RQC_Seq_Unit_Name', inplace=True)

sum_table_mean=sum_table
import_data(sum_table_mean, "./map_basecov_files/ecoli/", "ecoli", "mean")
import_data(sum_table_mean, "./map_basecov_files/pputida/", "pputida", "mean")
sum_table_mean.to_csv("summary_ep_meanBasecov.csv")

sum_table_median=sum_table
import_data(sum_table_median, "./map_basecov_files/ecoli/", "ecoli", "median")
import_data(sum_table_median, "./map_basecov_files/pputida/", "pputida", "median")
sum_table_median.to_csv("summary_ep_medianBasecov.csv")

