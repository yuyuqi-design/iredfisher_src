def write_config(rec_name, lig_name, box_center, box_size):
    """
    create configuration file for docking
    @rtype: object
    @param rec_name: the name of the protein structure without extension. e.g. pr1
    @param lig_name: the name of the ligand file without extension. e.g. lig1
    @param box_center: the coordinate of center of the docking box. e.g. ['x', 'y', 'z']
    @param box_size: the length of the cubic docking box in angstrom, e.g. '15'
    @return: none
    """
    # the name of the protein and ligand
    rec_pdbqt = f'{rec_name}.pdbqt'
    lig_pdbqt = f'{lig_name}.pdbqt'
    # coordinate of the center of the docking box
    center_x = str(box_center[0])
    center_y = str(box_center[1])
    center_z = str(box_center[2])
    # write the configuration file for docking
    with open('config.txt', 'w') as config:
        # the receptor file declare
        config.write(f'receptor = {rec_pdbqt} \n')
        # the ligand file declare
        config.write(f'ligand = {lig_pdbqt} \n')
        config.write('\n')
        # the docking box center declare
        config.write(f'center_x = {center_x}\n')
        config.write(f'center_y = {center_y}\n')
        config.write(f'center_z = {center_z}\n')
        config.write('\n')
        config.write(f'size_x = {box_size}\n')
        config.write(f'size_y = {box_size}\n')
        config.write(f'size_z = {box_size}\n')
        config.write('\n')
        config.write('exhaustiveness = 10\n')
        # config.write('num_modes = 20')
