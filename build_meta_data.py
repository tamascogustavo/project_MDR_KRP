#!/user/bin/envpython3
"""
Author:GustavoTamasco
Scripttorun:

Thescripttakes:


"""

#importstatements

import os.path
import subprocess
import os
from os import listdir
from os.path import isfile,join
import itertools
import shutil
import pandas as pd
import numpy as np

#functions and classes
def move_file(file_path, new_dir):
    '''
    This function will move the fna files to the created folder

    :param bof:
    :param new_dir:
    :return: none
    '''
    os.chdir(new_dir)
    out_name = file_path.strip().split("/")[-1]
    if os.path.exists(out_name):
        print("{} already in the dir".format(out_name))
    else:
        shutil.copy(file_path, new_dir)
        print("{} was moved to the dir".format(file_path))

def create_dir(path):
    '''
    This function creates a new directory in the local machine were the fna files will
    be stored
    :param path: path to the new dir
    :return: None
    '''
    if os.path.exists(path):
        print("Dir {} was already created".format(path))
    else:
        os.makedirs(path)
        print("{} was created".format(path))

def list_files(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = []
    files_ = [file for file in listdir(path) if isfile(join(path,file))]
    for file in files_:
        if ".snps" in file:
            files.append(file)
    return (files)

def parse_positions(snp_file, snp_positions):
    for line in snp_file:
        if "NODE" in line:
            position = line.strip().split()[0]
            snp_positions.add(position)
    return snp_positions


def get_specific_position(file, database, name, dict_info):
    all_positions = []
    database = list(database)
    dict_info.setdefault(name, [])
    for line in file:
        if "NODE" in line:
            position = line.strip().split()[0]
            all_positions.append(position)

    for i in range(len(database)):
        pos = database[i]
        if pos in all_positions:
            pos =1
            dict_info[name].append(pos)
        else:
            pos =0
            dict_info[name].append(pos)
    return dict_info













def main():
    """Main code of the script"""
    all_snp_positions = set()
    dict_info = {}
    #Getthefiles
    path_to_genomes = "/Users/gustavotamasco/mdrkrp/project_MDR_KRP_snps_info"
    dir_path = os.getcwd()
    data_dir = "{}_snps_info".format(dir_path)
    os.chdir(path_to_genomes)
    all_files = list_files(path_to_genomes)
    for file in all_files:
        with open(file) as snp_info:
            parse_positions(snp_info, all_snp_positions)

    for f in all_files:
        with open(f) as snp_info_2:
            metadata = get_specific_position(snp_info_2, all_snp_positions, f, dict_info)

    lol = [len(x) for x in metadata.values()]
    print(lol)
    df = pd.DataFrame.from_dict(metadata, orient='index')
    df = df.transpose()
    df.to_csv('metadata_df', sep='\t', encoding='utf-8')
    #print(df)












#main
if __name__=='__main__':
    main()
