import re


def clean_fasta_file(fasta_filename, data_folder):
    """
    Reads in a FASTA file, removes special characters from the heading,
    and limits the heading length to 12 characters and returns the cleaned
    FASTA sequences as a string.
    @param fasta_filename: a sequence file containing multiple sequences
    @param data_folder: the directory where the input data saved
    """
    clean_lines = []
    with open(f'{data_folder}/{fasta_filename}', 'r') as f:
        f_in = f.readlines()
        for line in f_in:
            if line.startswith('>'):  # heading line
                heading = line[1:].strip()  # remove the ">" and any whitespace
                heading = re.sub(r'[^_a-zA-Z0-9]', '', heading)  # remove special characters
                heading = heading[:12]  # limit the heading length to 12 characters
                heading = f'>{heading}\n'
                clean_lines.append(heading)
                  # write the cleaned heading to the output file
            else:  # sequence line
                 clean_lines.append(line)  # write the sequence to the output file
    with open(f'{data_folder}/allseq.fasta','w') as fout:
        fout.writelines(clean_lines)
