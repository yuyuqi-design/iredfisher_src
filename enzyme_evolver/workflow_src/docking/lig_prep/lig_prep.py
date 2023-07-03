import os
from pathlib import Path


def lig_prep(lig_structure) -> str:
    """
    prepare the ligand for docking using prepare_ligand4.py from autodock tools
    @param lig_structure: ligand structure, e.g. lig1.mol2
    @return: the name of the output structure e.g. lig1.pdbqt
    """

    # get the path of python scripts for docking
    file_abs_path = Path(__file__).resolve()
    script_folder = f'{file_abs_path.parents[1]}/scripts'

    # the name of the input structure
    lig_name = lig_structure.split('.')[0]

    structure_input = f'{lig_name}.mol2'
    # the name of the output structure
    structure_pdbqt = f'{lig_name}.pdbqt'

    # command for preparation of protein file
    comand_rec_prepare = f'python2.5 {script_folder}/prepare_ligand4.py -l {structure_input} -o {structure_pdbqt}'

    try:
        os.system(comand_rec_prepare)
    except:
        pass

    return structure_pdbqt
