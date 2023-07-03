import os
import pandas as pd
from urllib import request
import subprocess as sp
from typing import Tuple


def template_des(seqcode: str) -> Tuple[str, str]:
    """
    this function is to extract the information of the template protein by reading the fasta file
    . for example, the description/function and the origin organism
    """
    code = seqcode.strip()
    with open(code) as f:
        # read the fasta file
        all_codes = f.readlines()
        # extract the first line of fasta file
        firs_line = all_codes[0]
        # extract the description from first line
        organism = firs_line.strip().split('|')[3]
        description = firs_line.strip().split('|')[2]
    return description, organism


def find_template(seqcode):
    """
    This function is to select best template to build homology model
    """
    # prf including the alignment, sequence identiy, E-vaule etc information
    seqcode = seqcode.strip()
    result_file = seqcode + '.prf'
    # read potential template file from prf file
    df = pd.read_csv(result_file, skiprows=6, delim_whitespace=True,
                     names=["index", "PDB", "SX", "start_template", 'end_template', 'seq_startmatch', 'seq_endmatch',
                            'template_startmatch', 'template_endmatch', 'coverage', 'sequence_identity', 'Evalue',
                            'alignment'])
    # sort the template by the sequence identity
    df_sort = df.sort_values(by='sequence_identity', ascending=False)
    # the biggest sequence identity
    biggest_SI = float(df_sort.iloc[0, 10])
    # if template exist in PDB database
    if biggest_SI > 0.0:
        # create a spreadsheet for the detail for modelling
        df_code_seqID = df_sort[
                            ["PDB", "sequence_identity", 'seq_startmatch', 'seq_endmatch', 'template_startmatch',
                             'template_endmatch',
                             'coverage']][:10]
        df_code_seqID['template_id'] = df_code_seqID['PDB'].str[:4]
        df_code_seqID = df_code_seqID.drop_duplicates(subset=['template_id'])
        # get the top 3 template protein
        templates = list(df_code_seqID.iloc[:3, 0])
        # the best template, for example 5g6rA
        template1 = str(df_sort.iloc[0, 1])
        # pdb code of the best template, for example 5g6r
        template1_PdbCode = template1[0:4]
        # chain of the best template, for example, A
        template1_chain = template1[-1]
        # the second template
        if len(templates) >= 2:
            template2 = templates[1]
            template2_PdbCode = template2[0:4]
            template2_chain = template2[-1]
            download_fasta2 = f'wget https://www.rcsb.org/fasta/entry/' + template2_PdbCode
            sp.run(download_fasta2, shell=True)

        # if second template does not exist
        else:
            template2 = ''
            template2_PdbCode = ''
            template2_chain = ''
        # the third template
        if len(templates) >= 3:
            template3 = templates[2]
            template3_PdbCode = template3[0:4]
            template3_chain = template3[-1]
            download_fasta3 = f'wget https://www.rcsb.org/fasta/entry/' + template3_PdbCode
            sp.run(download_fasta3, shell=True)
        # if third template does not exist
        else:
            template3 = ''
            template3_PdbCode = ''
            template3_chain = ''
        # download  the sequence file of the first template protein
        download_fasta = f'wget https://www.rcsb.org/fasta/entry/' + template1_PdbCode
        sp.run(download_fasta, shell=True)
        try:
            request.urlretrieve('http://files.rcsb.org/download/' + template1_PdbCode + '.pdb',
                                template1_PdbCode + '.pdb')
            request.urlretrieve('http://files.rcsb.org/download/' + template2_PdbCode + '.pdb',
                                template2_PdbCode + '.pdb')
            request.urlretrieve('http://files.rcsb.org/download/' + template3_PdbCode + '.pdb',
                                template3_PdbCode + '.pdb')
        except:
            pass

    # if no available template
    else:
        templates = ''
        template1 = ''
        template1_PdbCode = ''
        template1_chain = ''
        template2 = ''
        template2_PdbCode = ''
        template2_chain = ''
        template3 = ''
        template3_PdbCode = ''
        template3_chain = ''
        f = open('info.csv', 'w')
        f.write(f'{seqcode}:no reliable template was found\n')
        f.close()
        print('Sorry, there is no reliable template  in the database\n')
        pass

    return templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, \
        template2_PdbCode, template3_chain, template3_PdbCode


def template_info(seqcode):
    seqcode = seqcode.strip()
    result_file = seqcode + '.prf'
    # read potential template file from prf file
    df = pd.read_csv(result_file, skiprows=6, delim_whitespace=True,
                     names=['index', 'PDB', 'SX', "start_template", 'end_template', 'seq_startmatch', 'seq_endmatch',
                            'template_startmatch', 'template_endmatch', 'coverage', 'sequence_identity', 'Evalue',
                            'alignment'])
    # sort the template by the sequence identity
    df_sort = df.sort_values(by='sequence_identity', ascending=False)
    # create a spreadsheet for template with collumns: PDB,sequence_identity,Evalue,coverage,
    # description, organism, seqcode,model
    # the biggest sequence identity
    df_best_template = df_sort.head(1)[['PDB', 'sequence_identity', 'Evalue', 'coverage']]
    template1_PdbCode = df_best_template['PDB'].to_string().split()[1][:4].strip()
    description, organism = template_des(template1_PdbCode)
    df_best_template['template'] = template1_PdbCode
    df_best_template['description'] = description
    df_best_template['organism'] = organism
    df_best_template['sequence'] = seqcode
    df_best_template['model'] = seqcode+'.pdb'
    first_column = df_best_template.pop('sequence')
    second_column_template_pdb = df_best_template.pop('PDB')
    df_best_template.insert(0, 'sequence', first_column)
    df_best_template.insert(1, 'template_PDB', second_column_template_pdb)
    return df_best_template


if __name__ == '__main__':
    # from pathlib import Path
    base_folder = os.getcwd()
    data_folder = '../tests/test_data'
    seqcode = f'test'
    os.chdir(data_folder)
    find_template(seqcode)
    os.chdir(base_folder)
