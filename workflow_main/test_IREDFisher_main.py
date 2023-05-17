import os
from enzyme_evolver.workflow_main.IREDFisher_main import IREDFisher
from enzyme_evolver.mongo.default_connection import make_default_connection


def test_IREDFisher_main():
    """
    test the IREDFisher workflow, the output file should be generated
     in ../database/5ed7fdd8-de19-4ec5-a97d-a355dac7cbd3
    """

    # connect to the mongo database
    make_default_connection()
    #input sequence
    fasta_filename = 'test_iredfisher_sequence.fasta'
    #input substrate of interest
    lig_structures = ['test_lig.pdb']
    # the ID of the N atom in imine bond
    N_ID = '7'
    # the working folder name under enzyme_evolver/database folder
    folder_id = '5ed7fdd8-de19-4ec5-a97d-a355dac7cbd3'

    # change to the parent folder to run iredfisher test
    os.chdir(f'../../')

    # run the workflow
    iredfisher = IREDFisher(fasta_filename=fasta_filename,lig_structures=lig_structures,N_ID=N_ID,
                            folder_id=folder_id, iredfisher_database='')
    iredfisher.run_IREDFisher()

    # return to this folder
    os.chdir(f'enzyme_evolver/workflow_main')
    #check if the iredfisher ranking file is genereated or not
    assert os.path.isfile(f'../database/5ed7fdd8-de19-4ec5-a97d-a355dac7cbd3/test_lig_IREDFisher_ranking.csv')  == True