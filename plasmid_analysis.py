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
        if "fna" in file:
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

    shutil.copy(file_path, new_dir)
    print("{} was moved to the dir".format(file_path))


def run_plasflow(genome_dir):
    '''
    This function runs parsnp
    :param dir_path: is the location in the start of the section
    :param genomes_dir: path to the fna files
    :return:
    '''

    out_name = genome_dir.split(".")[0]
    out_name = "{}_plasflow".format(out_name)

    if os.path.exists(out_name):
        print("Plasflow was already executed, {} already exist".format(out_name))
    else:
        cmd_plasflow = "PlasFlow.py --input {} --out {}".format(genome_dir, out_name)
        exit_message = subprocess.check_call(cmd_plasflow, shell=True)
        print("Exit status: {0}".format(exit_message))
        print("{0} Parsnp was executed".format(out_name))

def parse_fasta(file, name, dict):
    new_name = "{}_contigs.fasta".format(name.split('/')[-1].split(".")[0])
    all_seq = []
    for line in file:
        if line.startswith(">"):
            id = line.strip().split()
            id = "{}_{}".format(new_name,id)
        else:
            seq = line.strip()
            all_seq.append(seq)
    sequence = "".join(map(str, all_seq))
    dict[new_name] = sequence



def parse_info(file):
    swag = True
    seqs = {}
    all_info = []
    data = file.readlines()
    print('Name\tPlasmid len\tBit score and E-value\tIdent\tGaps')
    for i, item in enumerate(data):
        if ">" in item:
            name = [item.split(".")[0][1:]]
            info = data[i+1:i+5]
            info = [x.strip() for x in info]
            final = name+info
            print('{}\t{}\t{}\t{}\t{}'.format(final[0],final[1],final[2],final[3],final[4]))


def run_plasclass(organism):
    out_name = organism.split("_plasflow_plasmids")[0]
    out_name = "{}_statistics".format(out_name)
    if os.path.exists(out_name):
        print("{} already exist".format(out_name))
    else:
        cmd = "classify_fasta.py -f {} -o {}".format(organism, out_name)
        subprocess.check_call(cmd, shell=True)
        print("{} was generated".format(out_name))

def run_blastn(organism):
    out_name = organism.split("_plasflow_plasmids")[0]
    out_name = "{}_blast_result".format(out_name)
    if os.path.exists(out_name):
        print("{} already exist".format(out_name))
    else:
        cmd = "blastn –db nt –query {} –out {} -remote ".format(organism, out_name)
        subprocess.check_call(cmd, shell=True)
        print("{} was generated".format(out_name))



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
    path_to_all_info = '/Users/gustavotamasco/mdrkrp/project_MDR_KRPgenomes_parsnp'
    #path_to_all_info = argv[1]
    dirpath=os.getcwd()
    os.chdir(path_to_all_info)
    genome_files = list_directories(path_to_all_info)
    os.chdir("/Users/gustavotamasco/mdrkrp/plasmids")
    plasmid_files = list_directories("/Users/gustavotamasco/mdrkrp/plasmids")


    '''Genomes'''
    #for genome in genome_files:
        #if "fna" in genome:
            #print(genome)
            #run_plasflow(genome)

    '''Eval Plasmids'''
    for organism in plasmid_files:
        if "plasflow_plasmids" in organism:
            run_plasclass(organism)
            run_blastn(organism)



#main
if __name__=='__main__':
    main()
