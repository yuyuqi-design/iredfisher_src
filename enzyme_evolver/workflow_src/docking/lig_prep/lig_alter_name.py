from pymol import cmd


def lig_alter_name(lig_structure: str, N_ID: str) -> None:
    """
    this function is to change the residue name in the ligand structure file
    and change the atom N to N01  in C=N bond of imine structure.
    @param lig_structure: ligand structure file. e.g. LIG1.pdb
    @param N_ID: the atom ID of the active N atom
    @return: None
    """

    # load the structure into pymol
    cmd.load(filename=lig_structure)
    # select the N atom
    cmd.select('N_atom', f'id {N_ID}')
    # alter the atom name to N01
    cmd.alter('N_atom', "name='N01'")
    # alter the residue name to 'UNK'
    cmd.alter('all', "resn='UNK'")
    # save the structure
    cmd.save(lig_structure)
