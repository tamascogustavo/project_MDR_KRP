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


#functions and classes



def list_files(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path,file))]
    return (files)

def run_nucmer(query, ref):
    '''
    This function creates the delta files
    :param query: is a fasta file from a query genome
    :param ref: is a fasta file from a ref genome
    :return:
    '''
    out_name = "{}_{}".format(query.split(".")[0],ref.split(".")[0])
    out_name_nucmer = "{}.delta".format(out_name)
    if os.path.exists(out_name_nucmer):
        print("{} already exists".format(out_name_nucmer))
    else:
        cmd_nucmer = "nucmer --prefix={} {} {}".format(out_name, ref, query)
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

def move_delta(files, new_dir):
    '''
    This moves all delta files to a new dir
    :param files: all files name
    :param new_dir: path to the new dir
    :return:
    '''
    for file in files:
        if ".delta" in file:
            cmd = "mv {} {}".format(file,new_dir)
            subprocess.check_call(cmd, shell=True)

def move_snps(files, new_dir):
    '''
    This moves all snps files to a new dir
    :param files: all files
    :param new_dir: path to new dir
    :return:
    '''
    for file in files:
        if ".snps" in file:
            cmd = "mv {} {}".format(file,new_dir)
            subprocess.check_call(cmd, shell=True)

def main():
    """Main code of the script"""

    #Getthefiles
    path_to_genomes = "/Users/gustavotamasco/mdrkrp/project_MDR_KRPgenomes_parsnp_no431620"
    reference_genome = "MI_432.fna"
    dir_path = os.getcwd()
    data_dir = "{}_snps_info".format(dir_path)
    os.chdir(path_to_genomes)
    files = list_files(path_to_genomes)
    print("The file {} is used af ref for all the data manipulation".format(reference_genome))
    for file in files:
        if ".fna" in file:
            genome = file
            nuc_file = run_nucmer(genome, reference_genome)
            run_snp_detection(nuc_file)

    '''
    Move all files to a new dir 
    
    Only set this on if you have all done!
    '''
    #all_files = list_files(path_to_genomes)
    #move_delta(all_files, data_dir)
    #move_snps(all_files, data_dir)








#main
if __name__=='__main__':
    main()
