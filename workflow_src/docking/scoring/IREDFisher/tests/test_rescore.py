from enzyme_evolver.workflow_src.docking.scoring.IREDFisher import rescore
import os


def test_rescore():
    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    binding_mode_file = 'IRED1_LIG1_binding_poses.pdb'
    iredfisher_score, vina_ranking = \
        rescore.run_rescore(binding_mode_file=binding_mode_file, base_folder=base_folder, data_folder=data_folder)

    assert iredfisher_score == -33.4
    assert vina_ranking == 1
