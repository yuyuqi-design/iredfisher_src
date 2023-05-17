from enzyme_evolver.workflow_src.homology import batch_single_template_modelling


def test_batch_single_template_modelling():
    """
    test the whole workflow to generate the model based on the given template
    """
    base_folder = '/home/g02808yy/data/webserver/NC/enzyme_evolver/workflow_src/homology/tests'
    data_folder = base_folder + '/test_data'
    fasta_filename = 'allseq.fasta'
    template = '5g6r.pdb'
    batch_single_template_modelling.batch_single_template_modelling(fasta_filename=fasta_filename,
                                                                    base_folder=base_folder, data_folder=data_folder,
                                                                    template=template)
