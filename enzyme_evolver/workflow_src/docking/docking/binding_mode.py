import subprocess as sp
from pymol import cmd
from pathlib import Path


def binding_mode(rec_name: str, lig_name: str):
    """
    split the binding mode into individual object by pymol
    @param rec_name:
    @param lig_name:
    @return:
    """
    file_abs_path = Path(__file__).resolve()
    script_folder = f'{file_abs_path.parents[1]}/scripts'
    try:
        sp.run(f'sh {script_folder}/binding_mode.sh {rec_name} {lig_name}', shell=True)
        mode_name = f"{rec_name}_{lig_name}_binding_poses.pdbqt"
        obj_name = f"{rec_name}_{lig_name}_binding_poses"
        cmd.load(f"{mode_name}")
        cmd.split_states(f"{obj_name}")
        cmd.delete(f"{obj_name}")
        cmd.multisave(f"{obj_name}.pdb")
        cmd.delete('all')
        # sp.run(f"rm {folder}{mode_name}",shell=True)
        sp.run(f"sed -i '/HEADER/d' {obj_name}.pdb", shell=True)
    except:
        pass
