from pymol import cmd
from pymol import stored
from typing import Tuple


def residues_8a(binding_mode_file: str = 'IR-29_1a_binding_poses.pdb',
                ligname: str = 'UNK', best_ranking: int = 2) -> Tuple[int, int, int, int, list]:
    """
    This function is to get residues from protein file around 8 angstrom of the bound ligand
    @param binding_mode_file: the protein file with ligand bound. e.g. IR-29_1a_binding_poses.pdb
    @param ligname: the name of the ligand. e.g UNK
    @param best_ranking: the ranking of the best pose
    @return: the number of acidic residues, basic residues, histidine residues, tyrosine residues and a list
            of residues from the protein within 8 angstrom of the bound ligand
    """

    # create a list for the storage of residues within 8 angstrom of the ligand
    stored.list = []
    # load the structure file
    cmd.load(f"{binding_mode_file}")
    # object name of the structure file in pymol
    if len(cmd.get_object_list()) == 1:
        obj = binding_mode_file.split('.')[0]
    elif int(best_ranking) < 10:
        obj = binding_mode_file.split('.')[0] + '_000'+str(best_ranking)
    else:
        obj = binding_mode_file.split('.')[0] + '_00'+str(best_ranking)
    # select ligand from the pdb file by residue name ligname
    cmd.select('lig', obj + ' and resn ' + ligname)
    # select r for residues within 8 angstrom of ligand
    cmd.select('r', 'lig around 8')
    # select only named alpha atoms from r and save them into stored.list
    cmd.iterate("(name ca and r)", "stored.list.append((resi,resn))")
    # count basic and acidic and histidine residue numbers
    basic_number = 0
    acid_number = 0
    his_number = 0
    tyr_number = 0
    for i in stored.list:
        # print(i[1])
        # if i[1] == 'ASP' or i[1] == 'GLU' or i[1] == 'HIS' or i[1] == 'TYR':
        # if there is aspartate or glutamate,acidic number is added by 1
        if i[1] == 'ASP' or i[1] == 'GLU':
            acid_number = acid_number + 1
        # if there is lysine or arginine, basic number is added by 1
        if i[1] == 'LYS' or i[1] == 'ARG':
            basic_number = basic_number + 1
        # if there is a histidine, histidine number is added by 1
        if i[1] == 'HIS':
            his_number = his_number + 1
        if i[1] == 'TYR':
            tyr_number = tyr_number + 1
    # delete all objects from pymol
    cmd.delete('all')
    # return all acidic, basic, histidine residue number and the list of all residues
    return acid_number, basic_number, his_number, tyr_number, stored.list
