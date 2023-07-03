from math import sqrt
import subprocess as sp
from typing import List


def distance_two_atoms(c1: list, c2: list) -> float:
    """
    this function is to calculate the distance between coordinate 1 and 2
    @param c1: coordinate of the first atom
    @param c2: coordinate of the second atom
    @return: distance: the distance between the two atoms
    """
    distance = sqrt(
        (float(c1[0]) - float(c2[0])) ** 2 + (float(c1[1]) - float(c2[1])) ** 2 + (float(c1[2]) - float(c2[2])) ** 2)
    return distance


def distance_c4_n(binding_mode_file: str, string1: str = 'C4N NAP', string2: str = 'N01 UNK') -> List[float]:
    """
    this function is to calculate the distance between C4 atom in NADPH and the N atom in the ligand
    @param binding_mode_file: structure file containing protein and ligand poses
    @param string1: The string matching the C4 atom in NADPH
    @param string2: the string matching the N atom in the ligand
    @return: a list of distances for all docking poses
    """

    # get the coordinate of C4 atom of NAP from the structure file
    get_p1_line = f"grep --no-filename {string1} {binding_mode_file}"
    p1_line = sp.run(get_p1_line, capture_output=True, encoding="utf-8", shell=True)
    p1 = p1_line.stdout[32:54].split()

    # get the coordinate of N atom of the cofactor
    get_p2_line = f"grep --no-filename {string2} {binding_mode_file}"
    p2_line = sp.run(get_p2_line, capture_output=True, encoding="utf-8", shell=True)
    p2s = p2_line.stdout

    # make a list to save C4-N distances in each docking pose
    all_dist = []
    for line in p2s.split('\n'):
        if line:
            p2 = line[32:54].split()
            all_dist.append(distance_two_atoms(p1, p2))

    # turn every element in the list into float type
    all_dist_float = [float(i) for i in all_dist]

    if all_dist:
        # get the distance of the first docking pose
        first = all_dist[0]
        # get the minimum C4-N distance.
        minimum = min(all_dist_float)
        # get the rank of the minimum distance
        index_min = all_dist_float.index(minimum)
        # return the first distance, minimum distance, the rank of the minimum disntance and distances
        # for all docking poses
    return all_dist_float

