from enzyme_evolver.workflow_src.homology import batch_homodimer_modelling


def test_batch_homodimer_modelling():
    """
    test the whole workflow to generate the dimeric model based on the given dimeric template 'dimerT.pdb'
    and monomeric template 'singleT.pdb'
    """
    base_folder = '/home/g02808yy/data/webserver/NC/enzyme_evolver/workflow_src/homology/tests'
    data_folder = base_folder + '/test_data'
    fasta_filename = 'allseq.fasta'

    batch_homodimer_modelling.batch_homodimer_modelling(fasta_filename=fasta_filename, base_folder=base_folder,
                                                        data_folder=data_folder)
