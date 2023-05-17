from enzyme_evolver.workflow_src.homology import batch_multiple_modelling


def test_batch_best_multiple_modelling():
    """
    test the whole workflow to generate the model based on the top 3 template from
    PDB tructure database
    """
    base_folder = '/home/g02808yy/data/webserver/NC/enzyme_evolver/workflow_src/homology/tests'
    data_folder = base_folder + '/test_data'
    fasta_filename = 'allseq.fasta'
    batch_multiple_modelling.batch_auto_multiple_modelling(fasta_filename=fasta_filename, base_folder=base_folder,
                                                           data_folder=data_folder, db_folder=base_folder)
