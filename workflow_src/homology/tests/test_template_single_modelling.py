import os
from enzyme_evolver.workflow_src.homology.modelling import template_single_modelling


def test_template_single_modelling():
    """
        test the whole workflow to generate the model based on a given template structure
    """
    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    seqcode = 'test'
    # start the run
    os.chdir(data_folder)
    df_best_template = template_single_modelling.template_single_modelling(seqcode=seqcode, template='5g6r.pdb',
                                                                           base_folder=base_folder,
                                                                           data_folder=data_folder)
    os.chdir(base_folder)
    assert os.path.isfile(f'{data_folder}/{seqcode}.pdb') == True
    assert df_best_template['sequence_identity'].to_string(index=False).strip() == '100'
