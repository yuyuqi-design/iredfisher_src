from enzyme_evolver.workflow_src.docking import IREDFisher_scoring
import os


def test_iredfisher_rescore():
    """
    this function is to test the rescoring by IREDFisher of all complex structure from docking
    """

    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    rec_names = ['5ocm_pro', 'IRED1', 'IRED2']
    lig_structures = ['LIG1.pdb', 'LIG2.pdb']
    IREDFisher_scoring.iredfisher_scoring(rec_names=rec_names, lig_structures=lig_structures, base_folder=base_folder,
                                     data_folder=data_folder)

    out_file1 = f'{data_folder}/LIG1_IREDFisher_ranking.csv'
    out_file2 = f'{data_folder}/LIG2_IREDFisher_ranking.csv'

    assert os.path.isfile(out_file1) == True
    assert os.path.isfile(out_file2) == True

