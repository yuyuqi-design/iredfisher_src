import subprocess as sp
from pathlib import Path


def split_sequences(fasta_filename='allseq.fasta', data_folder=''):
    """
    this function is to split a fasta file containing multipe sequences into each individual ones
    @param fasta_filename: the input sequence file e.g. allseq.fasta
    @return: a list of sequence code. e.g. ['test1', 'test2', 'test3']
    """
    file_abs_path = Path(__file__).resolve()
    script_folder = f'{file_abs_path.parents[1]}/shell_scripts'
    spec_char_script = f'{script_folder}/rm_specialChara.sh'
    sp.run(f'sh {spec_char_script} {data_folder}/{fasta_filename}', shell=True)

    # split the fasta sequences into individual files and save all fasta titles into codes.txt file
    splitSeq_script = f'{script_folder}/split_sequences2.sh'
    # print(command_rm)
    command_split = f'sh {splitSeq_script} {data_folder}/{fasta_filename} {data_folder} > {data_folder}/codes0.txt'

    sp.run(command_split, shell=True)
    sp.run(f"sort {data_folder}/codes0.txt |uniq > {data_folder}/codes.txt",shell=True)
    with open(f'{data_folder}/codes.txt') as file_codes:
        codes = file_codes.readlines()

    return codes

def batch_template_info(df_best_template, datafolder='.'):
    """
    this function reads a dataframe containing the information of the template protein and output them into a csv file
    @param df_best_template: a dataframe containing the information of the template protein
    @return: none
    """
    df_best_template.to_csv(f'{datafolder}/template_info.csv',index=False)
