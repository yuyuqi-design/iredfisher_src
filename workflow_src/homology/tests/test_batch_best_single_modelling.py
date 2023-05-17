from enzyme_evolver.workflow_src.homology import batch_best_single_modelling


def test_batch_best_single_modelling():
    """
    test the whole workflow to generate the model based on the best template from
    PDB tructure database
    """
    base_folder = '/home/g02808yy/data/webserver/NC/enzyme_evolver/workflow_src/homology/tests'
    data_folder = base_folder + '/test_data'
    fasta_filename = 'allseq.fasta'
    batch_best_single_modelling.batch_auto_single_modelling(fasta_filename=fasta_filename, base_folder=base_folder,
                                                            data_folder=data_folder, db_folder=base_folder)
