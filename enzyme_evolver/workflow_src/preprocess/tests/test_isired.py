from enzyme_evolver.workflow_src.preprocess import isired


def test_isired():
    ired_code = isired.isired(seqcode='test',
                             db='/home/g02808yy/data/webserver/NC/enzyme_evolver/workflow_src/homology/tests')

    assert ired_code == 'test'
