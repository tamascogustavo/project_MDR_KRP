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
matplotlib.use("TkAgg")
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


def get_total_snp_value(name, file, snp_dict):
    '''
    This function calculates the total of snps in a genome comparison
    :param name: name of the file
    :param file: is the file
    :param snp_dict: a dict that will contain the info
    :return: an updated info
        :key is the name
        :value is the total number of snps
    '''
    size = 0
    for line in file:
        if "NODE" in line:
            size+=1
    snp_dict[name] = size
    return snp_dict

def print_total_snps(snp_data):
    '''
    This function prints the info of the dict of snps

    :param snp_data: dict with snp data
    :return: none
    '''
    for k,v in snp_data.items():
        name = k.split(".")[0]
        snp_total = v
        print("The comparison {} has {} as total number of snps ".format(name,snp_total))

def parse_info(name, file, coord_dict):
    '''
    This function creates a list with name and update a dict with new name and coordinates
    of snps

    Here we decided to use log to generate a unique code for each snp position related to
    the genome reference.
    Appying log10(ref+query)

    :param name: name of the old file
    :param file: is the .snps file
    :param coord_dict: a dict of snp coord associated with name
    :return: new name
    '''
    out_name = "{}_coords".format(name.split(".")[0])
    coords_list = []
    for line in file:
        if "NODE" in line:
            position_ref = float(line.strip().split()[0])
            position_query = float(line.strip().split()[3])
            position_value = position_ref
            coords_list.append(position_value)
            coord_dict[out_name] = coords_list
        if len(coords_list) == 0:
            coord_dict[out_name] = [1]
    return (out_name)


def plot_ven(list_names, data):
    '''
    This function plots a ven diagram

    :param list_names: list containing 3 names
    :param data: a dict build previously
    :return: None
    '''

    list1 = data[list_names[0]]
    list2 = data[list_names[1]]
    list3 = data[list_names[2]]

    venn3([set(list1), set(list2), set(list3)], set_labels=(list_names[0], list_names[1], list_names[2]))
    plt.show()

def get_data(names, info_dict):
    '''
    Get the info for a given name

    :param names: name of the file
    :param info_dict: dict with snps position info
    :return: dict specific for the request
    '''
    data_dir = {}
    for name in names:
        new_name = name.split(".")[0]
        new_name = "{}_MI_432_coords".format(new_name)
        if new_name in info_dict.keys():
            data_dir[new_name] = info_dict[new_name]
    return data_dir

def plot_group(data):
    '''
    Plot the ven diagram for manual selected groups
    :param data: is a dict where k is name and v is the snp position
    :return:
    '''

    all_names = [x for x in data.keys()]
    for i in range(0, len(all_names), 3):
        group_for_ven = all_names[i:i+3]
        if len(group_for_ven) == 3:
            list1 = data[group_for_ven[0]]
            list2 = data[group_for_ven[1]]
            list3 = data[group_for_ven[2]]

            venn3([set(list1), set(list2), set(list3)], set_labels=(group_for_ven[0], group_for_ven[1], group_for_ven[2]))
            plt.show()
        elif len(group_for_ven) ==2:
            list1 = data[group_for_ven[0]]
            list2 = data[group_for_ven[1]]

            venn2([set(list1), set(list2)],
                  set_labels=(group_for_ven[0], group_for_ven[1]))
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
    Check if something is missing
    '''

    for n in all_names:
        if n not in all_coords.keys():
            print(n)

    '''
    Creating the venn diagram full auto
    
    '''
    #for i in range(0,len(all_names),3):
        #group_for_ven = all_names[i:i+3]
        #plot_ven(group_for_ven, all_coords)

    '''
    Creating manual ven for specific groups
    '''
    #1) Three genomes in the top of the tree
    group_1 = ["genome24.fna", "Hemo_128.fna", "genome48.fna"]
    group_1_data = get_data(group_1, all_coords)
    print("----------------------------------------------------------")
    for k, v in group_1_data.items():
        print("The comparison {} has a total of {} snps".format(k, len(v)))
    #plot_group(group_1_data)

    #2)Organism in the middle of the tree
    group_2 = ['genome36.fna', "MI_432.fna","MI_45.fna",
               "URO_353.fna", "MIII_62.fna", "MI_209.fna",
               'MI_569.fna', "URO_1219.fna"]

    group_2_data = get_data(group_2, all_coords)
    print("----------------------------------------------------------")
    for k2, v2 in group_2_data.items():
        print("The comparison {} has a total of {} snps".format(k2, len(v2)))
    #plot_group(group_2_data)

    #3) Organism in the botton

    group_3 = ['Hemo_736.fna','URO_199.fna','URO_425.fna',
               'MI_88.fna', 'genome41.fna', 'genome1.fna',
               'MI_330.fna', 'MI_345.fna', 'URO_110.fna',
               'URO_775.fna', 'MI_17.fna', 'MI_536.fna',
               'MI_119.fna', 'MI_41.fna',  'URO_401.fna',
               'MI_186.fna', 'MI_449.fna', 'MI_306.fna',
               'Hemo_719.fna', 'URO_461.fna', 'MI_91.fna',
               'Hemo_805.fna', 'URO_770.fna', 'Hemo_989.fna',
               'MI_78.fna', 'Hemo_536.fna', 'MI_329.fna',
               'genome36.fna', 'Hemo_825.fna', 'MI_569.fna']

    group_3_data = get_data(group_3, all_coords)
    print("----------------------------------------------------------")
    for k3, v3 in group_3_data.items():
        print("The comparison {} has a total of {} snps".format(k3, len(v3)))
    #plot_group(group_3_data)

    #3.1) As parsnp botton
    group_3_parsnp = ['Hemo_736.fna', 'URO_199.fna', 'URO_425.fna',
                      'MI_88.fna', 'genome41.fna', 'genome1.fna',
                      'MI_330.fna', 'MI_345.fna', 'URO_110.fna',
                      'URO_775.fna', 'MI_17.fna', 'MI_536.fna',
                      'MI_119.fna', 'MI_41.fna', 'URO_401.fna',
                      'MI_186.fna', 'MI_449.fna', 'MI_306.fna',
                      'Hemo_719.fna', 'URO_461.fna', 'MI_91.fna',
                      'Hemo_805.fna', 'URO_770.fna']
    '''
    To check the differences
    '''
    res = [x for x in group_3 + group_3_parsnp if x not in group_3 or x not in group_3_parsnp]
    # print(res)
    group_3_parsnp_data = get_data(group_3_parsnp, all_coords)
    plot_group(group_3_parsnp_data)










#main
if __name__=='__main__':
    main()
