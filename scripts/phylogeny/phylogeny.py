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

#functionsandclasses
def list_directories(path):
    '''
    This function lists all file in a the path
    :param path:path of a dir
    :return: allfiles path in a list
    '''
    all_dirs = listdir(path)
    return all_dirs


def list_files(all_paths, path,organism):
    '''
    This function lists all file in a the path

    :param path:path of adir
    :return:all files path in a list

    '''
    complete_path = "{}/{}/ANNOT/PROKKA".format(path,organism)
    files = list_directories(complete_path)
    for file in files:
        if ".fna" in file:
            file_path = "{}/{}".format(complete_path, file)
            all_paths.append(file_path)
    return all_paths

def print_status(files):
    '''
    This function prints info about the files selected

    :param files: List of fna files
    :return: None
    '''

    message = "A total of {} file containing .fna in their name were found." \
              "Next step, building a directory of genomes and running parsnp!".format(len(files))
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


def run_parsnp(dir_path,genomes_dir):
    '''
    This function runs parsnp
    :param dir_path: is the location in the start of the section
    :param genomes_dir: path to the fna files
    :return:
    '''

    out_name = dir_path.split("/")[0:-1]
    out_name = "/".join(map(str,out_name))
    out_name = "{}/parsnp_out".format(out_name)
    if os.path.exists(out_name):
        print("Parsnp was already executed, {} already exist".format(out_name))
    else:
        cmd_parsnp = "parsnp -r ! -d {} -c -o {}".format(genomes_dir, out_name)
        exit_message = subprocess.check_call(cmd_parsnp, shell=True)
        print("Exit status: {0}".format(exit_message))
        print("{0} Parsnp was executed".format(out_name))

def list_files_new_source(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path,file))]
    return (files)

def main():
    """Main code of the script"""

    #Getthefiles
    all_fna_file_path = []
    path_to_all_info = '/Users/gustavotamasco/Google Drive/Shared drives/Projeto MDR KRP/Dados_Sequenciamento/'
    dirpath=os.getcwd()
    os.chdir(path_to_all_info)
    directories = list_directories(path_to_all_info)

    '''Genomes'''
    genomes_path = "{}{}".format(path_to_all_info,directories[0])
    os.chdir(genomes_path)
    genome_dir = list_directories(genomes_path)
    for organism in genome_dir:
        fna_files = list_files(all_fna_file_path,genomes_path,organism)
    print_status(fna_files)

    '''Building a dir of fna files'''
    genomes_fna_path = "{}genomes_parsnp".format(dirpath)
    create_genomes_dir(genomes_fna_path)
    os.chdir(genomes_fna_path)
    for file in fna_files:
        move_file(file, genomes_fna_path)

    '''Adding extra organism from a different source'''
    klebs = "/Users/gustavotamasco/mdrkrp/klebs"
    k_files = list_files_new_source(klebs)
    for k_file in k_files:
        if ".fna" in k_file:
            final_k_file = "{}/{}".format(klebs,k_file)
            move_file(final_k_file, genomes_fna_path)


    '''Run parsnp'''
    run_parsnp(dirpath, genomes_fna_path)




#main
if __name__=='__main__':
    main()
