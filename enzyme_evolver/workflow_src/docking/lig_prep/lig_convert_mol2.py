import subprocess as sp
from typing import Tuple


def lig_convert_mol2(lig_structure, add_h=False, gen_3D=False) -> Tuple:
    """
    convert input ligand structure to mol2 format by babel
    @param lig_structure: the name of the lig. e.g. lig1.pdb or lig1.sdf etc
    @param add_h: add hydrogens to the input ligand. default is False
    @param gen_3D: generate 3D conformation for the input ligand. default is False
    @return: the name of the ligand without extension. e.g. lig1
    @return: the output structure in mol2 structure. e.g. lig1.mol2
    """

    # the name of the input structure
    input_structure = str(lig_structure)
    lig_name = input_structure.split('.')[0]
    # the name of the output structure
    mol2_structure = f'{lig_name}.mol2'

    # whether to add hydrogens or generate 3D conformations
    if (add_h) and (not gen_3D):
        # add hydrogen and not generate 3D conformations
        command_babel = f'babel {input_structure} {mol2_structure} -h -p 7'
        sp.run(command_babel, shell=True)
    if (not add_h) and (gen_3D):
        # not add hydrogen and  generate 3D conformations
        command_babel = f'babel {input_structure} {mol2_structure} --gen3D -p 7'
        sp.run(command_babel, shell=True)
    if (add_h) and (gen_3D):
        # add hydrogen and generate 3D conformations
        command_babel = f'babel {input_structure} {mol2_structure} --gen3D -h -p 7'
        sp.run(command_babel, shell=True)
    if (not add_h) and (not gen_3D):
        command_babel = f'babel {input_structure} {mol2_structure}'
        sp.run(command_babel, shell=True)

    return lig_name, mol2_structure
