import modeller


def fasta2pir(seqcode):
    """
    This function is to convert the sequence from FASTA to PIR format which is
    recognized by Modeller
    """
    seqcode = str(seqcode).strip()
    e = modeller.environ()
    a = modeller.alignment(e, file=seqcode + '.fasta', alignment_format='FASTA')
    a.write(file=str(seqcode).strip() + '_0.pir', alignment_format='PIR')

    tempofile = str(seqcode).strip() + '_0.pir'
    filename = str(seqcode).strip() + '.pir'
    with open(tempofile) as f:
        with open(filename, 'w') as f2:
            seq = f.readlines()[3:]
            seq = ''.join(line.strip() for line in seq)
            # print(seq)
            seq_pir = f'>P1;{seqcode}\n' \
                      f'sequence:{seqcode}:::::::0.00: 0.00\n' \
                      f'{seq}'
            f2.write(seq_pir)


if __name__ == '__main__':
    # from pathlib import Path
    seqcode = f'../tests/test_data/test'
    fasta2pir(seqcode)
