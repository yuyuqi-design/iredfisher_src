from enzyme_evolver.workflow_src.docking.scoring.IREDFisher import residues_8a
import os

def test_residues_8a():
    """
    test the funtion residues_8a.py
    """

    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    os.chdir(data_folder)
    binding_mode_file = 'IRED1_LIG1_binding_poses.pdb'
    acid_number, basic_number, his_number, tyr_number, residues= \
        residues_8a.residues_8a(binding_mode_file=binding_mode_file, ligname='UNK', best_ranking=1)
    os.chdir(base_folder)

    assert acid_number == 2
    assert basic_number == 0
    assert his_number == 1
    assert tyr_number == 2
    assert len(residues) > 0