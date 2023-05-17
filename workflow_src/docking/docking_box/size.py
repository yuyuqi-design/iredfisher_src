from pathlib import Path
import subprocess as sp


def size(lig_structure: str, Rg=1.7) -> str:
    """
    calculate the size of the box by using the gyration radius of the input ligand
    @return: the length of the cubic docking box in angstrom, string type
    """
    # get the path of perl scripts for calculation of box size
    file_abs_path = Path(__file__).resolve()
    script_folder = f'{file_abs_path.parents[1]}/scripts'

    # apply the script eBoxSize.pl to calculate the gyration radius of the input ligand
    box_cmd = f'perl {script_folder}/eBoxSize.pl {lig_structure}'
    gyration_radius = sp.run(box_cmd, shell=True, stdout=sp.PIPE).stdout.strip().decode("utf-8")
    # set the length of docking box
    size = Rg * float(gyration_radius)
    size = str(size)

    return size
