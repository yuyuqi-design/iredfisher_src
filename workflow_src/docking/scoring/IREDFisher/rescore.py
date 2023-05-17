from enzyme_evolver.workflow_src.docking.scoring.IREDFisher import residues_8a
from enzyme_evolver.workflow_src.docking.scoring.IREDFisher import best_pose
from typing import Tuple
import os


def run_rescore(binding_mode_file: str, base_folder: str, data_folder: str) -> Tuple[float, int]:
    """
    rescore the docking pose
    @param binding_mode_file: the file containing protein and docking poses
    @param base_folder: the directory where you execute the function
    @param data_folder: the directory where the input data saved
    @return: iredfisher_score and vina_ranking
    """

    # enter the data folder
    os.chdir(data_folder)

    # get the vina score and ranking of the pose with the smallest C4-N distance and distance
    if os.path.isfile(binding_mode_file) == True:
        vina_score, vina_ranking, distance = best_pose.best_pose(binding_mode_file=binding_mode_file)

        # get the number of acidic, basic residues, histidine residues, tyrosine residues and a list of all residues within
        # 8 angstrom of the binding pose
        acid_number, basic_number, His_number, tyr_number, all_residues = residues_8a.residues_8a(
            binding_mode_file=binding_mode_file, ligname='UNK', best_ranking=vina_ranking)

        His_apprearence = 0
        # if histidine appear in the residues, his_number=1
        if His_number > 0:
            His_apprearence = 1
        # otherwise give a penalty of -5
        if His_number == 0:
            His_apprearence = 0

        # the IREDFisher scoring function
        iredfisher_score = 4.0 * float(vina_score) + 1.0 * float(acid_number) - 9.0 * float(
                His_apprearence) + 9.0 * float(basic_number)
        iredfisher_score = "%.2f" % iredfisher_score
        iredfisher_score = float(iredfisher_score)
    else:
        iredfisher_score = 999.0
        vina_ranking = 99

    # return to the base folder
    os.chdir(base_folder)

    return iredfisher_score, vina_ranking
