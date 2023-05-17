import os
import pandas as pd
from enzyme_evolver.workflow_src.homology.modelling import auto_multiple_modelling
from enzyme_evolver.workflow_src.homology.batch import batch


def batch_auto_multiple_modelling(fasta_filename, base_folder, data_folder, db_folder='enzyme_evolver/database'):
    """
    this function is to generate the model based on the top 3 template from
    PDB tructure database in a batch mode
    @param fasta_filename:a sequence file containing multiple sequences
    @param base_folder: the home directory
    @param data_folder:the directory where the input data is saved
    @param db_folder:the directory where the pdb_95.bin is saved
    @return: a dataframe containing the template information for all of the sequences with column:
    # sequence, PDB,sequence_identity,Evalue,coverage, description, organism, seqcode,model, ramachandran_score
    """
    os.chdir(data_folder)
    codes = batch.split_sequences(fasta_filename=fasta_filename, data_folder=data_folder)
    df_template = pd.DataFrame()
    for seqcode in codes:
        # seqcode, template_PDB,sequence_identity,Evalue,coverage,description, organism, model, ramachandran_score
        df_template_code = auto_multiple_modelling.auto_multiple_modelling(seqcode=seqcode, base_folder=base_folder,
                                                                  data_folder=data_folder, db_folder=db_folder)
        df_template = pd.concat([df_template, df_template_code])

    batch.batch_template_info(df_best_template=df_template, datafolder=data_folder)

    return df_template, codes
