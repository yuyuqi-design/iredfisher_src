from enzyme_evolver.workflow_src.docking.scoring import vina_scores


def test_vina_scores():
    binding_mode_file = 'test_data/IRED1_LIG1_binding_poses.pdb'
    all_score_float = vina_scores.vina_scores(binding_mode_file=binding_mode_file)
    assert len(all_score_float) == 9