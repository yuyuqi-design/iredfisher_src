from enzyme_evolver.workflow_src.docking.scoring.IREDFisher import distance_c4_n

def test_distance_c4_n():
    binding_mode_file = 'test_data/IRED1_LIG1_binding_poses.pdb'
    all_dist_float = distance_c4_n.distance_c4_n(binding_mode_file=binding_mode_file)
    assert len(all_dist_float) == 9
