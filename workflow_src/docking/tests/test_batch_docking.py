from enzyme_evolver.workflow_src.docking import batch_docking
import os


def test_batch_docking():
    """
    test run docking in batch mode
    """
    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    rec_names = ['5ocm_pro', 'IRED1', 'IRED2']
    lig_structures = ['LIG1.pdb', 'LIG2.pdb']
    ref_ligand = 'ref_lig.pdb'
    ref_cof = 'ref_cof.pdb'
    add_h = True
    gen_3D = False
    Rg = 1.7
    batch_docking.batch_docking(rec_names=rec_names, lig_structures=lig_structures, ref_ligand=ref_ligand,
                                base_folder=base_folder, data_folder=data_folder, cof_structure=ref_cof,
                                add_h=add_h, gen_3D=gen_3D, Rg=Rg)
