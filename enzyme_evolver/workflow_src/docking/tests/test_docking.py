from enzyme_evolver.workflow_src.docking.docking import docking
import os


def test_docking():
    """
    test the whole workflow to run docking
    """
    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data/single_docking'
    rec_name = '5ocm_pro'
    lig_structure = 'LIG.pdb'
    ref_ligand = 'ref_lig.pdb'
    cof_structure = 'ref_cof.pdb'
    add_h = True
    lig_name = lig_structure.split('.')[0]
    # start the run
    os.chdir(data_folder)
    docking.run_docking(rec_name=rec_name, lig_structure=lig_structure, ref_ligand=ref_ligand, base_folder=base_folder,
                        data_folder=data_folder, cof_structure=cof_structure, add_h=add_h)
    os.chdir(base_folder)
    assert os.path.isfile(f'{data_folder}/{rec_name}_{lig_name}_binding_poses.pdbqt') == True
    assert os.path.isfile(f'{data_folder}/{rec_name}_{lig_name}_binding_poses.pdb') == True
