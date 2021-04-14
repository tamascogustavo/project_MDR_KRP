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
matplotlib.use("TkAgg")
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

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

def list_genome_files(all_paths, path,organism):
    '''
    This function lists all file in a the path

    :param path:path of adir
    :return:all files path in a list

    '''
    complete_path = "{}/{}/ANNOT/ABRICATE".format(path,organism)
    files = list_directories(complete_path)
    for file in files:
        if "_argannot.tsv" in file:
            file_path = "{}/{}".format(complete_path, file)
            all_paths.append(file_path)
    return all_paths


def print_status(files):
    '''
    This function prints info about the files selected

    :param files: List of fna files
    :return: None
    '''

    message = "A total of {} file containing _argannot.tsv in their name were found." \
              "Next step, building a directory of these files".format(len(files))
    print(message)
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

def list_files_new_source(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path,file))]
    return (files)

def parse_genes(name, arg_file, arg_db):
    all_genes = []
    name = name.split("_vfdb")[0]
    for line in arg_file:
        if line.strip().split()[4] != "GENE":
            gene = line.strip().split()[4]
            all_genes.append(gene)
    arg_db[name] = all_genes
    return arg_db

def parse_genes_v2(name, arg_file, arg_db):
    all_genes = []
    new_genes = []
    name = name.split("_vfdb")[0]
    if name not in arg_db:
        for line in arg_file:
            if line.strip().split()[4] != "GENE":
                gene = line.strip().split()[4]
                all_genes.append(gene)
        arg_db[name] = all_genes
    else:
        for new_line in arg_file:
            if new_line.strip().split()[4] != "GENE":
                genes_to_add = new_line.strip().split()[4]
                new_genes.append(genes_to_add)
        arg_db[name]+=new_genes


def get_all_genes(metadata):
    all_genes = []
    for k, v in metadata.items():
        for gene in v:
            all_genes.append(gene)
    return all_genes

def build_class_df(df_info, classes, metadata):
    for organism, genes in metadata.items():
        abundance = Counter(genes)
        #print(abundance)
        array = get_array(abundance, classes)
        df_info[organism] = array
    return df_info

def get_array(genes_dict, classes):
    array = []
    for x in classes:
        if x in genes_dict.keys():
            array.append(genes_dict[x])
        else:
            value = 0
            array.append(value)
    return (array)



def main():
    """Main code of the script"""

    #Getthefiles
    all_plasmid_path = []
    path_to_all_info = '/Users/gustavotamasco/mdrkrp/project_MDR_KRP/all_vir_files'
    dirpath=os.getcwd()
    os.chdir(path_to_all_info)
    files = list_files_simple(path_to_all_info)

    '''Building metadata'''
    #All genes to each genome
    metadata = {}
    for file in files:
        with open(file) as vir_info:
            parse_genes_v2(file, vir_info, metadata)


    #All genes that occured
    all_genes = sorted(set(get_all_genes(metadata)))
    print(all_genes)

    #All vir classess


    '''Build dataframe for the classes plot'''
    df_info = {}
    df_major_classes = build_class_df(df_info, all_genes, metadata)
    df = pd.DataFrame.from_dict(df_major_classes, orient='index', columns=['clpC', 'entA', 'entB', 'entE', 'entS', 'fepA', 'fepB', 'fepC', 'fepD', 'fepG', 'fimA', 'fimE', 'fyuA', 'irp1', 'irp2', 'mgtB', 'mgtC', 'ompA', 'pilA', 'xcpA/pilD', 'xcpR', 'yagV/ecpE', 'yagW/ecpD', 'yagX/ecpC', 'yagY/ecpB', 'yagZ/ecpA', 'ybtA', 'ybtE', 'ybtP', 'ybtQ', 'ybtS', 'ybtT', 'ybtU', 'ybtX', 'ykgK/ecpR'])
    #df = df.transpose()
    #df.to_csv('arg_genes.csv', sep='\t', encoding='utf-8')
    #sns.set(font_scale=0.65)
    #Need both
    #not_full = sns.clustermap(df, label='small', cmap="vlag", standard_scale=1, linewidths=0)
    full_plot = sns.clustermap(df, label='small', cmap="vlag", linewidths=0)
    #plt.title('Antibiotic resistance genes across 34 organism', fontsize=15)
    #sns.set(font_scale=1)
    plt.show()
    full_plot.savefig("genome_plasmid_vir.pdf", bbox_inches='tight')
    #not_full.savefig("genome_plasmid_vir_not_scalled.pdf", bbox_inches='tight')




#main
if __name__=='__main__':
    main()
