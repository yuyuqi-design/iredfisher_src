from enzyme_evolver.workflow_src.preprocess import clean_fasta_file
from enzyme_evolver.workflow_src.preprocess import batch_isired
from enzyme_evolver.workflow_src.homology import batch_homodimer_modelling
from enzyme_evolver.workflow_src.docking import IREDFisher_docking, IREDFisher_scoring
from enzyme_evolver.workflow_src.homology.batch import batch
from typing import List


def iredfisher_ranking(fasta_filename: str, lig_structures: list,  N_ID: str, base_folder: str,
                       data_folder: str, db='enzyme_evolver/database') -> List[str]:
    """
    this function is to run  dock a list of ligand structures to a list of protein structures in batch mode and
    calculate the IREDFisher score for each complex structure.
    @param fasta_filename: a sequence file containing multiple sequences
    @param lig_structures: ligand structures. e.g. ['lig1.pdb', 'lig2.pdb', ...]
    @param N_ID: the atom ID of the active N atom
    @param base_folder: the directory where you execute the function
    @param data_folder: the directory where the input data saved
    @return:
    """

    # reference protein as the first protein
    rec_names = ['ref_complex']

    #remove the special characters in the heading of sequences and save as allseq.fasta file
    clean_fasta_file.clean_fasta_file(fasta_filename=fasta_filename, data_folder=data_folder)
    # split the sequences in the fasta file
    codes = batch.split_sequences(fasta_filename='allseq.fasta', data_folder=data_folder)
    codes = [s.strip() for s in codes]
    # checking whether they are imine reductase (IRED) and retain the IRED sequences for the next step
    ired_codes = batch_isired.batch_isired(codes=codes, db=db, base_folder=base_folder, data_folder=data_folder)
    # write out the IRED codes
    with open(f'{data_folder}/ireds.txt', 'w') as ired:
        for i in ired_codes:
            ired.write(f'{i}\n')

    # run modelling of the imine reductase
    df_template, ired_codes = batch_homodimer_modelling.batch_homodimer_modelling(
        codes=ired_codes, data_folder=data_folder, base_folder=base_folder)

    # all proteins for docking
    rec_names = rec_names + ired_codes

    # run docking in batch mode
    IREDFisher_docking.iredfisher_docking(rec_names=rec_names, lig_structures=lig_structures, N_ID=N_ID,
                                          base_folder=base_folder, data_folder=data_folder)
    # run rescoring in batch mode
    IREDFisher_scoring.iredfisher_scoring(rec_names=rec_names, lig_structures=lig_structures, base_folder=base_folder,
                                          data_folder=data_folder)

    return rec_names
