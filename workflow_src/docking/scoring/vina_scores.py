import subprocess as sp
from typing import List


def vina_scores(binding_mode_file: str = 'IR-92_1a_best_mode.pdb') -> List[float]:
    """
    to get all vina scores of all docking poses
    @param binding_mode_file: the file containing protein and all docking poses
    @return: a list of vina scores
    """
    # get the lines containing the score value
    file = f'{binding_mode_file}qt'
    get_score_line = f"grep --no-filename 'REMARK VINA RESULT' {file}"
    score_line = sp.run(get_score_line, capture_output=True, encoding="utf-8", shell=True)
    score_lines = score_line.stdout

    # make a list to save scores
    all_score = []
    for line in score_lines.split('\n'):
        if line:
            score = line.split()[3]
            # print(score)
            all_score.append(score)
    # turn element in score list into float type
    all_score_float = [float(i) for i in all_score]

    return all_score_float
