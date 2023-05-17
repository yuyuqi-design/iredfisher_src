import os
import pandas as pd
from enzyme_evolver.workflow_src.homology.modelling import template_single_modelling
from enzyme_evolver.workflow_src.homology.batch import batch


def batch_single_template_modelling(fasta_filename='', template='', base_folder='', data_folder=''):
    """
    this function is to is to generate the model based on the given template protein in a batch mode
    @param fasta_filename: a sequence file containing multiple sequences
    @param base_folder: the home directory
    @param data_folder: the directory where the input data is saved
    @param template: the name of the template structure
    @return: df_best_template: a dataframe containing the template information for all of the sequences with column:
    # sequence, sequence_identity, ramachandran_score
    """
    os.chdir(data_folder)
    codes = batch.split_sequences(fasta_filename=fasta_filename, data_folder=data_folder)
    df_template = pd.DataFrame()
    for seqcode in codes:
        # seqcode, template_PDB,sequence_identity,Evalue,coverage,description, organism, model, ramachandran_score
        df_template_code = template_single_modelling.template_single_modelling(seqcode=seqcode, template=template,
                                                                    base_folder=base_folder, data_folder=data_folder)
        df_template = pd.concat([df_template, df_template_code])

    batch.batch_template_info(df_best_template=df_template, datafolder=data_folder)

    return df_template, codes
