import os
import shutil
import subprocess as sp
from enzyme_evolver.workflow_src.docking.rec_prep import put_cof, rec_prep
from enzyme_evolver.workflow_src.docking.lig_prep import lig_convert_mol2, lig_prep
from enzyme_evolver.workflow_src.docking.docking_box import center_by_ligand, size
from enzyme_evolver.workflow_src.docking.configuration import docking_configuration
from enzyme_evolver.workflow_src.docking.docking import run_vina, binding_mode


def run_docking(rec_name: str, lig_structure: str, ref_ligand: str, base_folder: str = '',
                data_folder: str = '', cof_structure: str = '', add_h: bool = False,
                gen_3D: bool = False, Rg: float = 1.7, iredfisher_database='') -> None:
    """
    run docking workflow: protein preparation,  ligand preparation, creating docking box and run vina
    @rtype: object
    @param rec_name: the name of the protein structure without extension. e.g. pr1
    @param lig_structure: the name of the ligand file with extension. eg. lig1.pdb, lig1.sdf
    @param ref_ligand: name of the reference ligand with extension. e.g. ref_lig.pdb
    @param base_folder: the directory where you execute the function
    @param data_folder: the directory where the input data saved
    @param cof_structure: name of the cofactor structure file with an extension. e.g. ref_cof.pdb
    @param add_h: whether to add hydrogens to the ligand. Default is False
    @param gen_3D: generate 3D conformation for the input ligand. default is False
    @param Rg: the ratio to the gyration radius of the input ligand. default is 1.7
    @param iredfisher_database: iredfisher structure database.
    @return: None
    """

    # make a folder for the docking
    lig_name = lig_structure.split('.')[0]
    docking_folder = f'{data_folder}/{rec_name}_{lig_name}'
    try:
        # make subdirectory based on code
        os.mkdir(f"{docking_folder}")
    except FileExistsError:
        shutil.rmtree(f'{docking_folder}')
        os.mkdir(f"{docking_folder}")

    # copy the files from data_folder into the docking_folder
    sp.run(f"cp {data_folder}/{rec_name}.pdb* {docking_folder}", shell=True)
    sp.run(f"cp {data_folder}/{lig_structure} {docking_folder}", shell=True)
    sp.run(f"cp {data_folder}/{ref_ligand} {docking_folder}", shell=True)
    if cof_structure:
        sp.run(f"cp {data_folder}/{cof_structure} {docking_folder}", shell=True)

    # enter the docking folder
    os.chdir(docking_folder)
    # put the cofactor into the protein if there is one
    if cof_structure:
        put_cof.put_cof(rec_name=rec_name, cof_structure=cof_structure)
    # protein preparation if the protein is not from the iredfisher database
    if not iredfisher_database:
        rec_prep.rec_prep(rec_name=rec_name)

    # ligand preparation
    # convert the input ligand to mol2 format
    lig_name, mol2_structure = lig_convert_mol2.lig_convert_mol2(lig_structure=lig_structure,
                                                                 add_h=add_h, gen_3D=gen_3D)
    # prepare the ligand
    lig_prep.lig_prep(lig_structure=mol2_structure)

    # calculate the docking box
    # calculate the coordinate of the center of the box
    box_center = center_by_ligand.center_by_ligand(ref_lig_file=ref_ligand)
    # calculate the size of this box
    box_size = size.size(lig_structure=mol2_structure, Rg=Rg)

    # write the configuration file for docking
    docking_configuration.write_config(rec_name=rec_name, lig_name=lig_name, box_center=box_center, box_size=box_size)

    # run autodock vina
    run_vina.run_vina(configure_file='config.txt', logfile='log.txt')

    # write the binding modes
    binding_mode.binding_mode(rec_name=rec_name, lig_name=lig_name)

    # copy the binding mode to the datafolder
    sp.run(f'cp {rec_name}_{lig_name}_binding_poses.pdb {data_folder}', shell=True)
    sp.run(f'cp {rec_name}_{lig_name}_binding_poses.pdbqt {data_folder}', shell=True)

    # get back to the base folder
    os.chdir(base_folder)
