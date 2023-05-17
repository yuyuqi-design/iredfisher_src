from enzyme_evolver.workflow_src.preprocess import batch_isired


def test_isired():
    ired_code = batch_isired('')

    assert ired_code == 'test'