from enzyme_evolver.workflow_src.docking.scoring.IREDFisher import rescore
from typing import List
import pandas as pd


def iredfisher_scoring(rec_names: List[str], lig_structures: List[str], base_folder: str, data_folder: str) -> None:
    """
    this function is to run rescore for docking poses in a batch mode
    @param rec_names: protein names. e.g. ['pro1', 'pro2', ...]
    @param lig_structures: ligand structures. e.g. ['lig1.pdb', 'lig2.pdb', ...]
    @param base_folder: the directory where you execute the function
    @param data_folder: the directory where the input data saved
    @return: data_folder: the directory where the input data saved
    """

    # rescore each complex structure,
    for lig_structure in lig_structures:
        # get the name of the ligand structure
        lig_name = lig_structure.split('.')[0]
        # create a dataframe of scores from ligand to proteins
        df_lig_name = pd.DataFrame()

        # run rescore
        for rec_name in rec_names:
            #get the sequence
            # Open sequence file for reading
            with open(f'{data_folder}/{rec_name}.fasta', 'r') as file:
                # Read the contents of the file into a string
                rec_sequence = file.read()
            #rescore
            binding_mode_file = f'{rec_name}_{lig_name}_binding_poses.pdb'
            iredfisher_score, vina_ranking = rescore.run_rescore(binding_mode_file=binding_mode_file,
                                                                 base_folder=base_folder, data_folder=data_folder)
            vina_ranking = int(vina_ranking)
            # append the values to the dataframe
            df_lig_name = df_lig_name.append({'Enzyme': rec_name, 'IREDFisher_score': iredfisher_score,
                                              'Binding_mode_vina_ranking': vina_ranking,  'Sequence': rec_sequence},
                                             ignore_index=True)

        # sort the dataframe by the IREDFisher_score
        df_lig_name = df_lig_name.sort_values(by='IREDFisher_score', ascending=True)
        df_lig_name = df_lig_name.reset_index(drop=True)

        # layout
        # first column is the Enzyme name
        df_lig_name.insert(0, 'Enzyme', df_lig_name.pop('Enzyme'))
        # second column is the iredfisher score
        df_lig_name.insert(1, 'IREDFisher_score', df_lig_name.pop('IREDFisher_score'))
        # third column is the vina ranking number
        df_lig_name.insert(2, 'Binding_mode_vina_ranking', df_lig_name.pop('Binding_mode_vina_ranking'))
        df_lig_name.insert(3, 'Sequence', df_lig_name.pop('Sequence'))
        # write out the dataframe
        out_name = f'{data_folder}/{lig_name}_IREDFisher_ranking.csv'
        df_lig_name.to_csv(out_name)


if __name__ == "__main__":
    data_folder = '.'
    base_folder = '.'
    with open(f'{data_folder}/receptors_rm.txt') as f:
        rec_names = [line.strip() for line in f.readlines()]
        lig_structures = ['7b']
        iredfisher_scoring(rec_names=rec_names, lig_structures=lig_structures,
                           base_folder=base_folder, data_folder=data_folder)
