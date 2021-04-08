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
import math
import matplotlib
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
from matplotlib import pyplot as plt




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

    combinations = itertools.combinations(files,6)
    return (list(combinations))

def get_total_snp_value(name, file, snp_dict):
    size = 0
    for line in file:
        if "NODE" in line:
            size+=1
    snp_dict[name] = size
    return snp_dict

def print_total_snps(snp_data):
    for k,v in snp_data.items():
        name = k.split(".")[0]
        snp_total = v
        print("The comparison {} has {} as total number of snps ".format(name,snp_total))

def parse_info(name, file, coord_dict):
    out_name = "{}_coords".format(name.split(".")[0])
    coords_list = []
    for line in file:
        if "NODE" in line:
            position_ref = float(line.strip().split()[0])
            position_query = float(line.strip().split()[3])
            position_value = str(math.log10(position_ref+position_query))
            coords_list.append(position_value)
            coord_dict[out_name] = coords_list
        if len(coords_list) == 0:
            coord_dict[out_name] = 1
    return (out_name)


def plot_ven(list_names, data):

    list1 = data[list_names[0]]
    list2 = data[list_names[1]]
    list3 = data[list_names[2]]

    venn3([set(list1), set(list2), set(list3)], set_labels=("lol", "show", "hi"))
    plt.show()



def main():
    """Main code of the script"""

    snp_total = {}
    #Getthefiles
    path_to_snps = "/Users/gustavotamasco/mdrkrp/project_MDR_KRP_snps_info"
    dir_path = os.getcwd()
    data_dir = "{}_snps_stats".format(dir_path)
    os.chdir(path_to_snps)
    snp_files = list_files(path_to_snps)
    '''
    Total snps
    '''
    for file in snp_files:
        with open(file) as snp_data:
            get_total_snp_value(file,snp_data, snp_total)
    print_total_snps(snp_total)

    '''
    Deep stats for ven Diagram
    '''
    all_names = []
    all_coords = {}
    snp_files_ = list_files(path_to_snps)
    for file_ in snp_files_:
        with open(file_) as snp_ven:
            name = parse_info(file_,snp_ven, all_coords)
            all_names.append(name)
    '''
    The comparison 432 vs 432 is removed, because its comparing to itself
    # To check unmark the following code
    '''
    for n in all_names:
        if n not in all_coords.keys():
            print(n)

    '''
    Creating the venn diagram 
    
    '''
    for i in range(0,len(all_names),3):
        group_for_ven = all_names[i:i+3]
        plot_ven(group_for_ven, all_coords)











#main
if __name__=='__main__':
    main()
