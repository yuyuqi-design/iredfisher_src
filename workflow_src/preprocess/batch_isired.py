import os
import subprocess as sp
from enzyme_evolver.workflow_src.preprocess import isired
from typing import List


def batch_isired(codes: list, db: str, base_folder: str, data_folder: str) -> List[str]:
    """
    this function is search the best-matched template for the sequence and check whether the
    template is imine reductase in a batch mode.
    @param codes: a list of names of the protein sequence
    @param db: the pdb structure database for searching the template structures
    @param base_folder: the directory where you execute the function
    @param data_folder: the directory where the input data saved
    @return: a list of names of the ired protein sequence
    """
    ired_codes = []
    #limit the number of sequences to rank to 100
    codes = codes[:100]
    for code in codes:
        # make a folder to create alignment between each sequence and the pdb database
        code_folder = f'{data_folder}/{code}'
        if not os.path.isdir(code_folder):
            os.mkdir(code_folder)
        # copy the fasta file to the subfolder
        sp.run(f'cp {data_folder}/{code}.fasta {code_folder}', shell=True)
        # enter the code_folder
        os.chdir(code_folder)
        # check whether the sequence is imine reductase
        ired = isired.isired(seqcode=code, db=f'../..')
        # return the base_folder
        os.chdir(base_folder)
        # if it is imine reductase,  keep the sequence
        if ired:
            ired_codes.append(ired)
        # remove the subfoldder
        sp.run(f'rm -rf {code_folder}', shell=True)

    return ired_codes
