import os
from enzyme_evolver.workflow_src.homology.modelling import homodimer_modelling


def test_homodimer_modelling():
    """test the generation of homology models for homodimer based on given template structure"""
    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    seqcode = 'test'
    # start the run
    os.chdir(data_folder)
    df_best_template = homodimer_modelling.homodimer_modelling(seqcode=seqcode, base_folder=base_folder,
                                                               data_folder=data_folder)
    # end the run
    print(df_best_template)
    # test if the output file exsit
    assert os.path.isfile(f'{data_folder}/{seqcode}.pdb') == True
