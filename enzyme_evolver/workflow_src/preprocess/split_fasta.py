import sys

def split_fasta(infile):
    seq_id = None
    seq = []

    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                if seq_id is not None:
                    write_seq(seq_id, seq)
                    seq.clear()

                seq_id = line[1:13]
            else:
                seq.append(line)

        if seq_id is not None:
            write_seq(seq_id, seq)

def write_seq(seq_id, seq):
    outfile = f'{seq_id}.fasta'
    with open(outfile, 'w') as f:
        f.write(f'>{seq_id}\n')
        f.write(''.join(seq) + '\n')
