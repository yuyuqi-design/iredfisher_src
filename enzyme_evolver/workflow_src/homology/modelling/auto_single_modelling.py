import os
import shutil
from enzyme_evolver.workflow_src.homology.modelling import homology
import subprocess as sp
import pandas as pd


def auto_single_modelling(seqcode: str, base_folder: str, data_folder: str,
                          db_folder='enzyme_evolver/database') -> pd.DataFrame:
    """
    test the whole workflow to generate the model based on the best template from
    PDB tructure database and collect the information of template protein
    @param seqcode
    @param base_folder: the directory where you execute the function
    @param data_folder: the directory where the input data saved
    @param db_folder: the directory where pdb_95.bin saved
    @return df_best_template: a dataframe containing the template information for all of the sequences with column:
    # sequence, PDB,sequence_identity,Evalue,coverage, description, organism, seqcode,model, ramachandran_score
    """
    # seq_folder is the folder for running the modelling job for the seqcode.fasta
    seqcode = seqcode.strip()
    seq_folder = f'{data_folder}/{seqcode}'
    try:
        # make subdirectory based on code
        os.mkdir(f"{seq_folder}")
    except FileExistsError:
        shutil.rmtree(f'{seq_folder}')
        os.mkdir(f"{seq_folder}")
    # copy the sequence file from data_folder into the seq_folder
    sp.run(f"cp {data_folder}/{seqcode}.fasta {seq_folder}", shell=True)
    # enter the modelling directory
    os.chdir(seq_folder)
    # run modelling and collect the template information
    df_best_template = homology.auto_single_modelling(seqcode=seqcode, db=db_folder)
    os.chdir(base_folder)
    # copy the modelled structure to data_folder
    sp.run(f'cp {seq_folder}/{seqcode}.pdb {data_folder}', shell=True)
    # remove seq_folder
    sp.run(f'rm -r {seq_folder}', shell=True)

    return df_best_template
