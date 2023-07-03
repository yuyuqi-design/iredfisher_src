import os
import shutil
from enzyme_evolver.workflow_src.homology.modelling import homology
import subprocess as sp


def template_single_modelling(seqcode, template, base_folder, data_folder):
    """
    the function is to generate the model based on the given template protein
    @param seqcode: the sequence code. e.g. test
    @param template: the template structure. e.g. 5g6r.pdb
    @param base_folder: the home directory
    @param data_folder: the directory where the input data saved
    @return df_best_template: a dataframe containing the template information for all of the sequences with column:
    # sequence, sequence_identity, ramachandran_score
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
    # cpy the template file from data_folder into  the seq_folder
    sp.run(f"cp {data_folder}/{template} {seq_folder}", shell=True)
    # enter the modelling directory
    os.chdir(seq_folder)
    # run modelling and collect the template information
    df_best_template = homology.template_single_modelling(seqcode=seqcode, template=template)
    os.chdir(base_folder)
    # copy the modelled structure to data_folder
    sp.run(f'cp {seq_folder}/{seqcode}.pdb {data_folder}', shell=True)
    # remove seq_folder
    sp.run(f'rm -r {seq_folder}', shell=True)

    return df_best_template
