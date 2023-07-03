import os
import shutil
from enzyme_evolver.workflow_src.homology.modelling import homology
import subprocess as sp


def homodimer_modelling(seqcode, base_folder, data_folder, template_dimer='dimerT.pdb', template_monomer='singleT.pdb',
                        starting_residue='9', end_residue='294'):
    """
    this function is to generate the homodimeric model based on given dimeric and monomeric template structures
    @param seqcode
    @param base_folder: the directory where you execute the function
    @param data_folder: the directory where the input data saved
    @param template_dimer:  the name of the dimeric template structure file
    @param template_monomer: the name of the monomeric template structure file
    @param starting_residue: the index of the first residue in monomeric template structure file
    @param end_residue: the index of the last residue in monomeric template structure file
    @return df_best_template, a dataframe on the template information, sequence identity, ramachandran score etc.
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
    sp.run(f"cp {data_folder}/{template_dimer}  {seq_folder}", shell=True)
    sp.run(f"cp {data_folder}/{template_monomer}  {seq_folder}", shell=True)
    # enter the modelling directory
    os.chdir(seq_folder)
    # run modelling and collect the template information
    df_best_template = homology.homodimer_modelling(seqcode=seqcode, template_dimer=template_dimer,
                                                    template_monomer=template_monomer,
                                                    starting_residue=starting_residue, end_residue=end_residue)
    os.chdir(base_folder)
    # copy the modelled structure to data_folder
    sp.run(f'cp {seq_folder}/{seqcode}.pdb {data_folder}', shell=True)
    # remove seq_folder
    sp.run(f'rm -r {seq_folder}', shell=True)

    return df_best_template
