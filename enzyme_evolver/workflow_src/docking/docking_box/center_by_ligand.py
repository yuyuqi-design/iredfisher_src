from pymol import cmd


def center_by_ligand(ref_lig_file) -> list:
    """
    use a reference ligand to set the center of the docking box
    @param: ref_lig_file: name of the reference ligand. e.g. ref_lig.pdb
    @return: the coordinate of the center of mass of the reference ligand. e.g. ['x', 'y' ,'z']
    """

    # name of the reference ligand
    ref_lig_name = ref_lig_file.split('.')[0]
    # load the reference ligand into pymol
    cmd.load(ref_lig_file)
    # get the coordinate of the box center
    box_center = cmd.centerofmass(ref_lig_name)
    # unload the file from pymol
    cmd.delete('all')

    return box_center
