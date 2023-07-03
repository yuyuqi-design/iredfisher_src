import modeller
import os


def create_alignment(seqcode: str, db='enzyme_evolver/database'):
    """
        This function is to make sequence alignment between user's sequence and
        proteins in PDB database (pdb_95.pir which is updated in the end of November, 2019.
    """
    seqcode = str(seqcode).strip()
    env = modeller.environ()

    # -- Read in the sequence database
    sdb = modeller.sequence_db(env)

    # -- Now, read in the binary database
    sdb.read(seq_database_file=db + '/pdb_95.bin', seq_database_format='BINARY',
             chains_list='ALL')

    # -- Read in the target sequence/alignment
    aln = modeller.alignment(env)
    aln.append(file=seqcode + '.pir', alignment_format='PIR', align_codes='ALL')

    # -- Convert the input sequence/alignment into
    #   profile format
    prf = aln.to_profile()

    # -- Scan sequence database to pick up homologous sequences
    prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
              gap_penalties_1d=(-500, -50), n_prof_iterations=1,
              check_profile=False, max_aln_evalue=0.01)

    # -- Write out the profile in text format
    prf.write(file=seqcode + '.prf', profile_format='TEXT')

    # -- Convert the profile back to alignment format
    aln = prf.to_alignment()

    # -- Write out the alignment file
    aln.write(file=seqcode + '.ali', alignment_format='PIR')


if __name__ == '__main__':
    # from pathlib import Path
    masterpath = os.getcwd()
    os.chdir(f'{masterpath}/enzyme_evolver/workflow_src_2/homology/tests/test_data')
    seqcode = 'test'
    create_alignment(seqcode, db=f'{masterpath}/enzyme_evolver/database/pdb_95.pir')
    os.chdir(masterpath)
