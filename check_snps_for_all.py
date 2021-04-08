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
        if ".fna" in file:
            files.append(file)
    return (files)

def run_nucmer(comb):
    '''
    This function creates the delta files
    :param query: is a fasta file from a query genome
    :param ref: is a fasta file from a ref genome
    :return:
    '''
    out_name = "{}_{}".format(comb[0].split(".")[0], comb[1].split(".")[0])
    out_name_nucmer = "{}.delta".format(out_name)
    if os.path.exists(out_name_nucmer):
        print("{} already exists".format(out_name_nucmer))
    else:
        cmd_nucmer = "nucmer --prefix={} {} {}".format(out_name, comb[0], comb[1])
        exit_message = subprocess.check_call(cmd_nucmer, shell=True)
        print("Exit status: {0}".format(exit_message))
        print("{0} nucmer was executed".format(out_name))
    return out_name_nucmer

def run_snp_detection(delta_file):
    '''
    This function generates a file with the position where snp happens
    :param delta_file: file name
    :return:
    '''
    out_name = delta_file.split(".")[0]
    out_name_snps = "{}.snps".format(out_name)
    if os.path.exists(out_name_snps):
        print("{} already exists".format(out_name_snps))
    else:
        cmd_snp_view = "show-snps -Clr {} > {}".format(delta_file, out_name_snps)
        exit_message = subprocess.check_call(cmd_snp_view, shell=True)
        print("Exit status: {0}".format(exit_message))
        print("Check the file {0} for snps info".format(out_name))
    return out_name_snps


def generate_combinations(files):

    combinations = itertools.combinations(files,2)
    return (list(combinations))


def main():
    """Main code of the script"""

    #Getthefiles
    path_to_genomes = "/Users/gustavotamasco/mdrkrp/genomes"
    dir_path = os.getcwd()
    data_dir = "{}_snps_info".format(dir_path)
    os.chdir(path_to_genomes)
    all_files = list_files(path_to_genomes)
    files_cluster_closeorgs = list_files(path_to_genomes)
    files_cluster_3 = ["genome24.fna","Hemo_128.fna","genome48.fna"]
    files_cluster_mid = ["URO_1219.fna", "MI_432.fna",
                         "MI_45.fna", "URO_353.fna",
                         "MIII_62.fna", "MI_209.fna"]
    for elements1 in files_cluster_3:
        if elements1 in files_cluster_closeorgs:
            files_cluster_closeorgs.remove(elements1)

    for elements2 in files_cluster_mid:
        if elements2 in files_cluster_closeorgs:
            files_cluster_closeorgs.remove(elements2)

    '''
    Create all possible combinations
    '''


    #For the cluster with 3
    all_names_1 = []
    new_path_1 = "{}/cluster_of_three".format(os.getcwd())
    create_dir(new_path_1)

    combination_of_tree_org = generate_combinations(files_cluster_3)
    for comb1 in combination_of_tree_org:
        nuc_file = run_nucmer(comb1)
        all_names_1.append(nuc_file)
        snp_file = run_snp_detection(nuc_file)
        all_names_1.append(snp_file)
    for n1 in all_names_1:
        n1 = "{}/{}".format(path_to_genomes,n1)
        move_file(n1, new_path_1)

    #For the mid cluster
    all_names_2 = []
    new_path_2 = "{}/mid_cluster".format(path_to_genomes)
    create_dir(new_path_2)

    os.chdir(path_to_genomes)
    combination_mid_org = generate_combinations(files_cluster_mid)
    for comb2 in combination_mid_org:
        nuc_file_2 = run_nucmer(comb2)
        all_names_2.append(nuc_file_2)
        snp_file_2 = run_snp_detection(nuc_file_2)
        all_names_2.append(snp_file_2)
    for n2 in all_names_2:
        n2 = "{}/{}".format(path_to_genomes,n2)
        move_file(n2, new_path_2)

    #For the closest organism in the tree
    os.chdir(path_to_genomes)
    all_names_3 = []
    new_path_3 = "{}/close_org_cluster".format(path_to_genomes)
    create_dir(new_path_3)
    combination_close_org = generate_combinations(files_cluster_closeorgs)
    for comb3 in combination_close_org:
        nuc_file_3 = run_nucmer(comb3)
        all_names_3.append(nuc_file_3)
        snp_file_3 = run_snp_detection(nuc_file_3)
        all_names_3.append(snp_file_3)
    for n3 in all_names_3:
        n3 = "{}/{}".format(path_to_genomes, n3)
        move_file(n3, new_path_3)

    '''
    
    '''








#main
if __name__=='__main__':
    main()
