from enzyme_evolver.workflow_src.docking import IREDFisher_docking
import os


def test_iredfisher_docking():
    base_folder = os.getcwd()
    data_folder = base_folder + '/test_data'
    rec_names = ['5ocm_pro', 'IRED1', 'IRED2']
    lig_structures = ['LIG2.pdb']
    ref_ligand = 'ref_lig.pdb'
    N_ID = '19'
    ref_cof = 'ref_cof.pdb'
    add_h = True
    gen_3D = False
    Rg = 1.7
    IREDFisher_docking.iredfisher_docking(rec_names=rec_names, lig_structures=lig_structures, ref_ligand=ref_ligand,
                                          N_ID=N_ID, base_folder=base_folder, data_folder=data_folder,
                                          cof_structure=ref_cof, add_h=add_h, gen_3D=gen_3D, Rg=Rg)
