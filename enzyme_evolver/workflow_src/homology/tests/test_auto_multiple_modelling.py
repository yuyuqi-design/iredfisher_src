import os
from enzyme_evolver.workflow_src.homology.modelling import auto_multiple_modelling


def test_auto_single_modelling():
    """
    test the whole workflow to generate the model based on the best template from
    PDB tructure database
    """
    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    seqcode = 'test'
    # start the run
    os.chdir(data_folder)
    df_best_template = auto_multiple_modelling.auto_multiple_modelling(seqcode=seqcode, base_folder=base_folder,
                                                                       data_folder=data_folder, db_folder=base_folder)
    os.chdir(base_folder)
    assert os.path.isfile(f'{data_folder}/{seqcode}.pdb') == True
    assert df_best_template['template'].values == '5g6r'
