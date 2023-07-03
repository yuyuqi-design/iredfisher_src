from enzyme_evolver.workflow_src.homology.seqformat import fasta2pir
from enzyme_evolver.workflow_src.homology.find_template import create_alignment, find_template
from enzyme_evolver.workflow_src.homology.alignment import align2template
from enzyme_evolver.workflow_src.homology.modelling import build_models, homology
from enzyme_evolver.workflow_src.homology.evaluation import evaluation
import os


def test_sequence_convert():
    """ test converting of fasta to pir format """
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the test
    os.chdir(data_folder)
    # convert sequence format
    fasta2pir.fasta2pir(seqcode)
    os.chdir(function_folder)
    # whether the pir file is generated
    assert os.path.isfile(f'{data_folder}/{seqcode}.pir') == True


def test_create_alignmnent():
    """test the sequence alignment to the sequences in PDB database"""
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the test
    os.chdir(data_folder)
    create_alignment.create_alignment(seqcode, db=f'{function_folder}')
    os.chdir(function_folder)
    assert os.path.isfile(f'{data_folder}/{seqcode}.ali') == True
    assert os.path.isfile(f'{data_folder}/{seqcode}.prf') == True


def test_find_template():
    """select the top 3 best template according to the sequence identity
    and test top 3 best template structure of the input sequence
    """
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the test
    os.chdir(data_folder)
    templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain,\
        template2_PdbCode, template3_chain, template3_PdbCode = find_template.find_template(seqcode)
    df_best_template = find_template.template_info(seqcode)
    os.chdir(function_folder)
    assert template1 == '5g6rA'
    assert template2 == '6eodF'
    assert template3 == '5ojlA'
    assert len(df_best_template) == 1


def test_template_info():
    """test the generation of homology models for homodimer based on given template structure"""
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the run
    os.chdir(data_folder)
    df_best_template = find_template.template_info(seqcode)
    os.chdir(function_folder)
    assert df_best_template['template'].values == '5g6r'


def test_align2d():
    """test the alignment of sequence to one template sequence"""
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the test
    os.chdir(data_folder)
    align2template.align2d(seqcode, template1='5g6rA', template1_PdbCode='5g6r', template1_chain='A')
    os.chdir(function_folder)
    assert os.path.isfile(f'{data_folder}/{seqcode}-5g6rA.ali') == True
    assert os.path.isfile(f'{data_folder}/{seqcode}-5g6rA.pap') == True


def test_align_multi():
    """test the alignment of sequence to top 3 template sequences"""
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the test
    os.chdir(data_folder)
    align2template.salign(seqcode, template1_PdbCode='5g6r', template1_chain='A', template2_PdbCode='6eod',
                          template2_chain='F', template3_PdbCode='5ojl', template3_chain='A')
    align2template.align2d_mult(seqcode)
    os.chdir(function_folder)
    assert os.path.isfile(f'{data_folder}/{seqcode}-multi.ali') == True
    assert os.path.isfile(f'{data_folder}/{seqcode}-multi.pap') == True


def test_align_dimer():
    """test the alignment of sequence to dimeric template protein"""
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the test
    os.chdir(data_folder)
    sequence_identity = align2template.align_dimer(seqcode, start_residue='9',
                                                   end_residue='294', N_chains='2')
    os.chdir(function_folder)
    assert os.path.isfile(f'{data_folder}/{seqcode}dimer.ali') == True
    assert sequence_identity == 60


def test_build_model():
    """test the generation of homology models based on one single template"""
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # tart the test
    os.chdir(data_folder)
    build_models.build_model(seqcode, template1='5g6rA')
    os.chdir(function_folder)
    # end the test
    assert os.path.isfile(f'{data_folder}/{seqcode}.pdb') == True


def test_build_multi_models():
    """test the generation of homology models based on top 3 template"""
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the run
    os.chdir(data_folder)
    build_models.build_multi_models(seqcode, template1='5g6rA', template2='6eodF', template3='5ojlA')
    os.chdir(function_folder)
    # end the run
    assert os.path.isfile(f'{data_folder}/{seqcode}.pdb') == True


def test_evaluation_models():
    """test the evaluation of homology models generated my modeller"""
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # start the run
    os.chdir(data_folder)
    evaluation.procheck(seqcode)
    ramachandran_score = evaluation.rama_core(seqcode)
    os.chdir(function_folder)
    # test the ramachandran score
    assert ramachandran_score > 90.0


def test_homodimer_modelling():
    """
        test the whole workflow to generate the homodimeric model based
        on a given template structure
    """
    function_folder = os.getcwd()
    data_folder = function_folder + '/test_data'
    seqcode = 'test'
    # tart the run
    os.chdir(data_folder)
    df_best_template = homology.homodimer_modelling(seqcode)
    os.chdir(function_folder)
    assert os.path.isfile(f'{data_folder}/{seqcode}.pdb') == True
    assert df_best_template['sequence_identity'].to_string(index=False).strip() == '60'
