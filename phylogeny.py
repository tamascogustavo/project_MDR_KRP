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


#functionsandclasses
def list_directories(path):
    '''
    This function lists all file in a the path
    :param path:path of a dir
    :return: allfiles path in a list
    '''
    all_dirs = listdir(path)
    return all_dirs


def list_files(path,organism):
    '''
    This function lists all file in a the path
    :param path:path of adir
    :return:all files path in a list
    '''

    complete_path = "{}/{}".format(path,organism)
    directories = list_directories(complete_path)
    print(directories)
    #files=[fileforfileinlistdir(complete_path)ifisfile(join(complete_path,file))]
    #print(listdir(complete_path))


def get_file_m(path):
    """Thisfunctioncreatesthepathforallyourfiles

    input:isastringwiththepathtoyourfoldercontainingallfiles
    Out:isalistofstringscontainingthepathofallfiles

    Obs:nocasopegeuisoreadspoissãoosarquivosondeosreadsnãoalinhadosforamremovidos
    """
    files=[]
    for r,d,f in os.walk(path):
        for file in f:
            if file.startswith("reads_")and"montagem"in file:
                files.append(os.path.join(r,file))
    return files

def get_file_final(path):
    """Thisfunctioncreatesthepathforallyourfiles

    input:isastringwiththepathtoyourfoldercontainingallfiles
    Out:isalistofstringscontainingthepathofallfiles

    Obs:nocasopegeuisoreadspoissãoosarquivosondeosreadsnãoalinhadosforamremovidos
    """
    files = []
    for r,d,f in os.walk(path):
        for file in f:
            if file.startswith("only")and"montagem"in file:
                files.append(os.path.join(r,file))
    return files


def main():
    """Main code of the script"""

    #Getthefiles
    path_to_all_info = '/Users/gustavotamasco/Google Drive/Shared drives/Projeto MDR KRP/Dados_Sequenciamento/'
    dirpath=os.getcwd()
    os.chdir(path_to_all_info)
    directories=list_directories(path_to_all_info)

    '''Genomes'''
    genomes_path = "{}{}".format(path_to_all_info,directories[0])
    os.chdir(genomes_path)
    genome_dir = list_directories(genomes_path)
    for organism in genome_dir:
        list_files(genomes_path,organism)

    #Converttofastq
    #fq_files=fastq_converter(files)
    #Getthefqfiles
    #dirpath=os.getcwd()
    #os.chdir(dirpath)
    #fq=get_fq(dirpath)
    #run_canu(fq)


#main
if __name__=='__main__':
    main()
