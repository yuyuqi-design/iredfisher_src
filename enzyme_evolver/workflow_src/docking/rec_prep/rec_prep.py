import os
from pathlib import Path


def rec_prep(rec_name: str) -> str:
    """
    prepare the protein for docking using prepare_receptor4.py from autodock tools
    @param rec_name: protein structure name, e.g. pro1
    @return: the name of the output structure
    """

    # get the path of python scripts for docking
    file_abs_path = Path(__file__).resolve()
    script_folder = f'{file_abs_path.parents[1]}/scripts'

    # the name of the input structure
    structure_input = f'{rec_name}.pdb'
    # the name of the output structure
    structure_output = f'{rec_name}.pdbqt'

    # command for preparation of protein file
    comand_rec_prepare = f'python2.5 {script_folder}/prepare_receptor4.py -r {structure_input} -A checkhydrogens ' \
                         f'-U deleteAltB -o {structure_output}'

    try:
        os.system(comand_rec_prepare)
        # sp.run(f'cp {pdbqt} /home/g02808yy/data/webserver/IREDFisher/enzyme_evolver/database/Public', shell=True)
    except:
        pass

    return structure_output
