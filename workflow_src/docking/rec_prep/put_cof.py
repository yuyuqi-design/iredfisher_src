from pathlib import Path
import subprocess as sp


def put_cof(rec_name: str, cof_structure: str = 'ref_cof.pdb') -> None:
    """
    put cofactor into the protein before the protein preparation
    @param rec_name: the protein structure name without an extension, e.g. pro1
    @param cof_structure: the cofactor structure name with the extension, e.g. ref_cof.pdb
    @return:
    """

    # the name of the protein structure
    rec_structure = f'{rec_name}.pdb'
    # remove the END tag in the last line of protein file
    # get the path of shell script to remove the tag
    file_abs_path = Path(__file__).resolve()
    script_folder = f'{file_abs_path.parents[1]}/scripts'
    # remove the tag in the protein using the shell script
    rm_tag = f'sh {script_folder}/rm_END.sh {rec_structure}'
    sp.run(rm_tag, shell=True)
    # remove the tag in the cofactor file using the shell script
    rm_tag2 = f'sh {script_folder}/rm_END.sh {cof_structure}'
    sp.run(rm_tag2, shell=True)

    # get the lines of the cofactor structure file
    with open(cof_structure) as cof:
        cof_str = cof.read()
    # open the protein structure file and append the cofactor lines to the end of the protein lines
    f = open(rec_structure, 'a')
    f.write(cof_str)
    f.close()
