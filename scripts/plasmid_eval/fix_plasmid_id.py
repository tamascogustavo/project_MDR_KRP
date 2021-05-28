#!/user/bin/envpython3
"""
Author:GustavoTamasco
Scripttorun:

Thescripttakes:

Runthecodeby:ex_p5.py<reference_fasta_file><related_fasta_file>
"""

#importstatements
from sys import argv
import os.path
import subprocess
import os
from os import listdir
from os.path import isfile,join
import shutil
from collections import Counter
import pandas as pd
import matplotlib


#functionsandclasses

def list_directories(path):
    '''
    This function lists all file in a the path
    :param path:path of a dir
    :return: allfiles path in a list
    '''
    all_dirs = listdir(path)
    return all_dirs

def list_files_simple(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = []
    files_ = [file for file in listdir(path) if isfile(join(path,file))]
    for f in files_:
        if ".tsv" in f:
            files.append(f)
    return (files)


def list_files(all_paths, path,organism):
    '''
    This function lists all file in a the path

    :param path:path of adir
    :return:all files path in a list

    '''
    complete_path = "{}/{}/ABRICATE".format(path,organism)
    files = list_directories(complete_path)
    for file in files:
        if "_argannot.tsv" in file:
            file_path = "{}/{}".format(complete_path, file)
            all_paths.append(file_path)
    return all_paths


def move_file(file_path, new_dir):
    '''
    This function will move the fna files to the created folder

    :param bof:
    :param new_dir:
    :return: none
    '''
    out_name = file_path.strip().split("/")[-1]
    if os.path.exists(out_name):
        print("{} already in the dir".format(out_name))
    else:
        shutil.copy(file_path, new_dir)
        print("{} was moved to the dir".format(file_path))

def list_files_new_source(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path,file))]
    return (files)

def parse_fasta(organism, file, all_plasmid_data):
    name = organism.strip().split("_plasflow_plasmids.fasta")[0]
    for line in file:
        if line.startswith(">"):
            seq_id = line.strip()[1:]
            seq_id = "{}_{}".format(name, seq_id)
            all_plasmid_data[seq_id] = ""
        else:
            all_plasmid_data[seq_id] += line.strip()



def main():
    """Main code of the script"""

    #Getthefiles
    all_plasmid_data = {}
    path_to_all_info = '/Users/gustavotamasco/mdrkrp/plasmids/plasmid_seqs'
    dirpath=os.getcwd()
    os.chdir(path_to_all_info)
    all_files = list_files_new_source(path_to_all_info)
    for organism in all_files:
        with open(organism) as file:
            parse_fasta(organism, file, all_plasmid_data)
    for k,v in all_plasmid_data.items():
        print(">{}\n{}".format(k,v))



#main
if __name__=='__main__':
    main()
