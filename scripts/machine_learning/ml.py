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
import re
import pandas as pd

#functionsandclasses
def list_directories(path):
    '''
    This function lists all file in a the path
    :param path:path of a dir
    :return: allfiles path in a list
    '''
    all_dirs = listdir(path)
    return all_dirs

def create_genomes_dir(path):
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


def run_dsk(dir_path,list_outs):
    '''
    This function runs parsnp
    :param dir_path: is the location in the start of the section
    :param genomes_dir: path to the fna files
    :return:
    '''

    out_name = dir_path.split("/")[-1][0:-4]
    out_name = "{}_dsk_out.h5".format(out_name)
    list_outs.append(out_name)
    if os.path.exists(out_name):
        print("dsk was already executed for this genome, {} already exist".format(out_name))
    else:
        cmd_dsk = "dsk -file {} -out {} -kmer-size 31 ".format(dir_path, out_name)
        print(cmd_dsk)
        exit_message = subprocess.check_call(cmd_dsk, shell=True)
        print("Exit status: {0}".format(exit_message))
        print("{0} dsk was executed".format(out_name))

def run_dsk_translator(file, kmer_list):
    out_name = file.split("_dsk_out.h5")[0]
    out_name = "{}_kmers_dataset.txt".format(out_name)
    kmer_list.append(out_name)
    if os.path.exists(out_name):
        pass
    else:
        cmd_trn = "dsk2ascii -file {} -out {}".format(file, out_name)
        exit_message = subprocess.check_call(cmd_trn, shell=True)
        print("Exist status: {}".format(exit_message))
        print("{} was created".format(out_name))

def parse_all_info(file, organism, organism_dict):
    '''
    This function creates a dictionary for organisms and their kmers

    The filter using count >=2 is the same used by the authors of the original
    paper
    :param file: is a openned file of kmers count
    :param organism: is the name of the file
    :param organism_dict: is a meta dict containing all info
    :return: updated dict
    '''
    organism_name = organism.split("_kmers_dataset.txt")
    organism_dict[organism] = []
    for line in file:
        count = int(line.strip().split()[-1])
        if count >=2:
            info = line.strip().split()[0]
            organism_dict[organism].append(info)
    return organism_dict


def list_files(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path,file))]
    return (files)

def build_database(organis_dict, all_kmers, metadata):
    if os.path.exists("all_dataset_kmers.csv"):
        print("all_dataset_kmers.csv already exists")
    else:
        all_kmers = list(all_kmers)
        names = [x for x in organis_dict.keys()]
        for name in names:
            final_name = name.split("_kmers_dataset.txt")[0]
            metadata.setdefault(final_name, [])

            for i in range(len(all_kmers)):
                pos = all_kmers[i]
                if pos in organis_dict[name]:
                    pos =1
                    metadata[final_name].append(pos)
                else:
                    pos = 0
                    metadata[final_name].append(pos)
    return metadata






def main():
    """Main code of the script"""

    #Getthefiles
    all_dsk_out = []
    all_kmers_isolates = []
    organism_dict = {}
    metadata = {}
    path_to_all_info = '/Users/gustavotamasco/mdrkrp/genomes_final/genomes'
    genomes = list_files(path_to_all_info)
    genomes = ["{}/{}".format(path_to_all_info,x) for x in genomes]
    #print(genomes)
    #generating kmers
    for genome in genomes:
        run_dsk(genome, all_dsk_out)
    #creating individual db
    for kmer in all_dsk_out:
        run_dsk_translator(kmer, all_kmers_isolates)
    #Building full metadata
    for isolate_kmer in all_kmers_isolates:
        with open(isolate_kmer) as file:
            organism_metadata = parse_all_info(file, isolate_kmer, organism_dict)
    #all kmers set
    all_kmer_set = [x for x in organism_dict.values()]
    flat_list_kmer = set([item for sublist in all_kmer_set for item in sublist])

    #Building dataframe
    '''
    columns = flat_list_kmer 
    rows = is the organism name followed by presence absence or each kmer
    
    '''
    build_database(organism_dict, flat_list_kmer, metadata)
    if os.path.exists("all_dataset_kmers.csv"):
        print("done")
    else:
        df = pd.DataFrame.from_dict(metadata, orient='index')
        #df = df.transpose()
        df.to_csv('all_dataset_kmers.csv', sep='\t', encoding='utf-8')

    #rename colnames
    if os.path.exists('mdkrp_dataset_kmers.csv'):
        print("all right for running ")
    else:
        data = pd.read_csv("all_dataset_kmers.csv", sep='\t')
        kmer_col = list(flat_list_kmer)
        col_names = ["isolate"] + kmer_col

        data.columns = list(col_names)

        data.to_csv('mdkrp_dataset_kmers.csv', sep='\t', encoding='utf-8')

    for x in organism_dict.keys():
        print(x)



#main
if __name__=='__main__':
    main()
