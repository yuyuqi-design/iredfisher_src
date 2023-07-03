import numpy as np


def box_center_residues(rec, residue_ID1, residue_ID2, residue_ID3, residue_ID4) -> list:
    """
    use 4 reference residues in protein to set the center of the docking box
    @param rec: the name of the protein file. e.g. pro1
    @param residue_ID1: the residue ID of first reference residue in protein
    @param residue_ID2: the residue ID of second reference residue in protein
    @param residue_ID3: the residue ID of third reference residue in protein
    @param residue_ID4: the residue ID of fourth reference residue in protein
    @return:the coordinate of the center of mass of the reference ligand. e.g. ['x', 'y' ,'z']
    """
    rec_file = f'{rec}.pdb'
    with open(rec_file) as receptor:
        # read the lines of the protein file
        content = receptor.readlines()
        # get the coordinate of the CA atoms of the four reference residues.
        for line in content:
            if (line.strip()[22:27].strip() == str(residue_ID1).strip()) & (line.strip()[13:15] == 'CA') & \
                    (line.strip().startswith('ATOM')):
                cord1 = line.strip()[32:55].split()
            #                     print(cord1)
            elif (line.strip()[22:27].strip() == str(residue_ID2).strip()) & (line.strip()[13:15] == 'CA') & \
                    (line.strip().startswith('ATOM')):
                cord2 = line.strip()[32:55].split()
            #                     print(cord2)
            elif (line.strip()[22:27].strip() == str(residue_ID3).strip()) & (line.strip()[13:15] == 'CA') & \
                    (line.strip().startswith('ATOM')):
                cord3 = line.strip()[32:55].split()
            #                     print(cord3)
            elif (line.strip()[22:27].strip() == str(residue_ID4).strip()) & (line.strip()[13:15] == 'CA') & \
                    (line.strip().startswith('ATOM')):
                cord4 = line.strip()[32:55].split()
        # get the coordinate of the center of four CA atoms
        cord = list(map(float, cord1)) + list(map(float, cord2)) + list(map(float, cord3)) + list(map(float, cord4))
        arry_cord = np.array(cord).reshape(4, 3)
        cord_aver = np.average(arry_cord, axis=0)
        box_center = list(cord_aver)

        return box_center
