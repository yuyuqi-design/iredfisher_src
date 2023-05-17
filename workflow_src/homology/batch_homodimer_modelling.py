import os
import pandas as pd
from enzyme_evolver.workflow_src.homology.modelling import homodimer_modelling
from enzyme_evolver.workflow_src.homology.batch import batch


def batch_homodimer_modelling(codes, base_folder, data_folder, template_dimer='dimerT.pdb',
                              template_monomer='singleT.pdb', starting_residue='9', end_residue='294'):
    """
    this function is to this function is to generate the homodimeric model based on given
     dimeric and monomeric template structures in a batch mode
    @param codes:a list of protein names for modelling
    @param base_folder: the home directory
    @param data_folder: the directory where the input data is saved
    @param template_dimer: the name of the dimeric template structure file
    @param template_monomer: the name of the monomeric template structure file
    @param starting_residue: the index of the first residue in monomeric template structure file
    @param end_residue: the index of the last residue in monomeric template structure file
    @return: a dataframe containing the template information for all of the sequences with column:
    # sequence, PDB,sequence_identity,Evalue,coverage, description, organism, seqcode,model, ramachandran_score
    """

    os.chdir(data_folder)
    df_template = pd.DataFrame()
    for seqcode in codes:
        # seqcode, template_PDB,sequence_identity,Evalue,coverage,description, organism, model, ramachandran_score
        df_template_code = homodimer_modelling.homodimer_modelling(seqcode=seqcode, base_folder=base_folder,
                        data_folder=data_folder, template_dimer=template_dimer, template_monomer=template_monomer,
                        starting_residue=starting_residue, end_residue=end_residue)

        df_template = pd.concat([df_template, df_template_code])

    batch.batch_template_info(df_best_template=df_template, datafolder=data_folder)

    return df_template, codes
