from enzyme_evolver.workflow_src.docking.scoring.IREDFisher import best_pose

def test_best_pose():
    binding_mode_file='test_data/ref_complex_marvinjs_untitled_file_binding_poses.pdb'
    Vina_score, best_mode_number, Distance = best_pose.best_pose(binding_mode_file=binding_mode_file)
    assert Vina_score == '10'
    assert Distance == '66'
    assert best_mode_number == 2