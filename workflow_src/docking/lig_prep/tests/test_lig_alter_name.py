from enzyme_evolver.workflow_src.docking.lig_prep.lig_alter_name import lig_alter_name


def test_lig_alter_name():
    """
    test the lig_alter_name.py
    @return:
    """
    lig_structure = 'LIG1.pdb'
    lig_alter_name(lig_structure=lig_structure, N_ID='6')
