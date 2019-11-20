#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:25:38 2019

@author: DVyshenska
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
#import numpy as np
from os import chdir, getcwd, listdir
#import math
wd=getcwd()
chdir(wd)

#%% FUNCTIONS

def files_list_f(pat_match, folder_path):
    return [f for f in listdir(folder_path) if pat_match in f]

def file_id_f(f_name, pat_name, pat_ext):
    return f_name.replace(pat_name, "").replace(pat_ext, "")
#%% VARIABLES

spikin_path = "./spikin_stats_files/"
summ_file_path = "./summary_files/Summary.Lib.Info_updated.csv"
map_path = "./map_stats_files/"

#%%

lib_info = pd.read_csv(summ_file_path)
lib_info['RQC_Seq_Unit_Name'].replace(regex=True,inplace=True,
        to_replace=r'.fastq.gz', value=r'')


spkn_col_names=["total_all_reads",
"total_spkn_reads",
"Batch146B_p061_GC-50.8_13Clabel-57.2",
"Batch146B_p060_GC-50.5_13Clabel-42.9",
"Batch146B_p059_GC-48.9_13Clabel-14.3",
"Batch146B_p063_GC-51.1_13Clabel-71.5",
"Batch146B_p058_GC-49.7_13Clabel-28.6",
"Batch146B_p056_GC-52.3_13Clabel-100",
"Batch146B_p057_GC-47.7_13Clabel-0",
"Batch146B_p064_GC-51.8_13Clabel-85.8",
"pSET152_GC-63_13Clabel-100",
"pW5Y-AprR_GC-37_13Clabel-0",
"pputida_assigned_reads",
"ecoli_assigned_reads"]


lib_info_updated = lib_info.reindex(columns = lib_info.columns.tolist() + \
                                    spkn_col_names)
lib_info_updated.set_index('RQC_Seq_Unit_Name', inplace=True)

for file_name in files_list_f("stats_", spikin_path):
    file_id = file_id_f(file_name, "stats_", ".txt")
    
    fileH = open(spikin_path + file_name, "r")
    for l in fileH:
        if "Total" in l:
            lib_info_updated.loc[file_id,"total_all_reads"]=int(l.split('\t')[1])
        elif "Matched" in l:
            lib_info_updated.loc[file_id,"total_spkn_reads"]=int(l.split('\t')[1])
        for s in spkn_col_names[2:]:
            if s in l:
                lib_info_updated.loc[file_id,s]=int(l.split('\t')[1])
    fileH.close()

#lib_info_updated.to_csv("summary_out_SP.csv")

lib_info_updated.fillna(0,inplace=True)
for file_name_m in files_list_f("mapstats_", map_path):
    file_id_m = file_id_f(file_name_m, "mapstats_", ".txt")
    temp_t = pd.read_csv(map_path + file_name_m, sep="\t")
    temp_t.set_index('#name', inplace=True)
    lib_info_updated.loc[file_id_m,"pputida_assigned_reads"] = int(temp_t.loc['pputida', 'assignedReads'])
    lib_info_updated.loc[file_id_m,"ecoli_assigned_reads"] = int(temp_t.loc['ecoli', 'assignedReads'])

lib_info_updated.to_csv("summary_out_wd.csv")
print("DONE!\n")
