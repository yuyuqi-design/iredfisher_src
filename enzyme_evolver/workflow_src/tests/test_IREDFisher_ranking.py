from enzyme_evolver.workflow_src import IREDFisher_ranking
import os


def test_IREDFisher_ranking():
    """
    test the function IREDFisher_ranking.py
    """
    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    # data_folder = '/home/g02808yy/data/webserver/NC/enzyme_evolver/database/ef2e9417-6253-40f3-9340-60ad67460023'
    fasta_filename = 'allseq.fasta'
    lig_structures = ['LIG2.pdb']
    N_ID = '19'

    IREDFisher_ranking.iredfisher_ranking(fasta_filename=fasta_filename, lig_structures=lig_structures,
                                          N_ID=N_ID, base_folder=base_folder, data_folder=data_folder,
                                          db='/home/g02808yy/data/webserver/NC/enzyme_evolver/database')

    out_file2 = f'{data_folder}/LIG2_IREDFisher_ranking.csv'
    assert os.path.isfile(out_file2) == True
