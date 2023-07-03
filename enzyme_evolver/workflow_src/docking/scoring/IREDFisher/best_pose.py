from enzyme_evolver.workflow_src.docking.scoring.IREDFisher import distance_c4_n
from enzyme_evolver.workflow_src.docking.scoring import vina_scores
import pandas as pd
import os
from typing import Tuple


def best_pose(binding_mode_file: str = 'IRED2_LIG1_binding_poses.pdb') -> Tuple[str, int, str]:
    """
    this function is to select the smallest C4-N distance after removing the poses with C4-N between 3.5 and 6.0 A.
    @param binding_mode_file: the file containing protein and docking poses, pdb format
    @return: Vina_score (the vina score), best_mode_number (ranking of the pose with the smallest C4-N distance),
            Distance (the smallest C4-N distance)
    """

    # get all C4-N distances in docking poses
    all_dist_float = distance_c4_n.distance_c4_n(binding_mode_file=binding_mode_file, string1='C4N NAP',
                                                 string2='N01 UNK')
    all_score_float = []
    data_tuples = ()
    if len(all_dist_float) > 0:
        # get the vina scores of poses
        all_score_float = vina_scores.vina_scores(f'{binding_mode_file}')
        # combine the C4-N distance and corresponding scores as tuples. e.g [(4.6, -6.0), (4.0, -6.5),...]
        data_tuples = list(zip(all_dist_float, all_score_float))
        # turn the tuples into dataframe with column names 'Distance', 'Score'
        df = pd.DataFrame(data_tuples, columns=['Distance', 'Vina_score'])

        # extract the poses with C4-N distance between 3.5 and 6.0
        df_dist = df[(df['Distance'] >= 3.5) & (df['Distance'] <= 6.0)]
        # rank the dataframe by C4-N distance in the remaining poses
        df_dist_sort = df_dist.sort_values(by='Distance')
        # get the smallest C4-N distance in the remaining poses
        Distance = (df_dist_sort.iloc[:1]['Distance']).to_string()

        # if there is no poses with C4-N between 3.5 and 6.0, the Distance is given a big penalty as 66
        if 'Series' in Distance:
            Distance = '0 66'

        # get the vina score of the pose
        Vina_score = (df_dist_sort.iloc[:1]['Vina_score']).to_string()
        # if there is no poses with C4-N between 3.5 and 6.0, the Vina_score is given a big penalty as 10
        if 'Series' in Vina_score:
            Vina_score = '0 10'

        # get the corresponding the ranking number of the selected pose
        best_mode_number = int(Distance.split()[0]) + 1

        # remove the index
        Vina_score = Vina_score.split()[1]
        Distance = Distance.split()[1]

        return Vina_score, best_mode_number, Distance
