from enzyme_evolver.workflow_src.docking.docking import docking
from enzyme_evolver.workflow_src.docking.batch import rec_align
from enzyme_evolver.workflow_src.docking.lig_prep import lig_alter_name
import os
import subprocess as sp

def iredfisher_docking(rec_names: list, lig_structures: list,  N_ID: str, base_folder: str, data_folder: str,
                       ref_ligand: str = 'ref_lig.pdb', cof_structure: str = 'ref_cof.pdb', add_h: bool = False,
                       gen_3D: bool = False, Rg: float = 1.7, iredfisher_database='') -> None:
    """
    this function is to dock a list of ligand structures to a list of protein structures in batch mode
    @param rec_names: protein names. e.g. ['pro1', 'pro2', ...]
    @param lig_structures: ligand structures. e.g. ['lig1.pdb', 'lig2.pdb', ...]
    @param ref_ligand: name of the reference ligand with extension. e.g. ref_lig.pdb
    @param N_ID: the atom ID of the N atom in C=N bond
    @param base_folder: the directory where you execute the function
    @param data_folder: the directory where the input data saved
    @param cof_structure: name of the cofactor structure file with an extension. e.g. ref_cof.pdb
    @param add_h: whether to add hydrogens to the ligand. Default is False
    @param gen_3D: generate 3D conformation for the input ligand. default is False
    @param Rg: the ratio to the gyration radius of the input ligand. default is 1.7
    @param iredfisher_database: iredfisher structure database.
    @return:
    """

    # align the protein structures to the first protein structure
    os.chdir(data_folder)
    if not iredfisher_database:
        rec_align.rec_align(rec_names)
    #change the residue name of the ligand and the N atom
    for lig_structure in lig_structures:
        lig_alter_name.lig_alter_name(lig_structure=lig_structure, N_ID=N_ID)
    os.chdir(base_folder)

    # docking ligand structures to protein structures
    for lig_structure in lig_structures:
        lig_name = lig_structure.split('.')[0]
        for rec_name in rec_names:
            docking.run_docking(rec_name=rec_name, lig_structure=lig_structure, ref_ligand=ref_ligand,
                                base_folder=base_folder, data_folder=data_folder, cof_structure=cof_structure,
                                add_h=add_h, gen_3D=gen_3D, Rg=Rg,iredfisher_database=iredfisher_database)

            # delete the docking directories
            docking_folder = f'{data_folder}/{rec_name}_{lig_name}'
            sp.run(f'rm -r {docking_folder}', shell=True)

