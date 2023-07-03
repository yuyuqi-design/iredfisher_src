from enzyme_evolver.workflow_src.homology.seqformat import fasta2pir
from enzyme_evolver.workflow_src.homology.find_template import create_alignment, find_template
import pandas as pd


def isired(seqcode='test', db='enzyme_evolver/database'):
    """
    this function is search the best-matched template for the sequence and check whether the
    template is imine reductase
    @param seqcode: name of the protein sequence
    @param db: the pdb structure database for searching the template structures
    @return: df_best_template: a dataframe on the template information, sequence identity, ramachandran score etc.
    """
    seqcode = seqcode.strip()
    # convert the sequence test.fasta to test.pir
    fasta2pir.fasta2pir(seqcode)
    # create alignment between the target sequence and the sequences in PDB database
    create_alignment.create_alignment(seqcode, db=db)
    # search the template
    templates, template1, template2, template3, template1_chain, template1_PdbCode, \
        template2_chain, template2_PdbCode, template3_chain, template3_PdbCode = find_template.find_template(seqcode)
    # if the best template is one of the imine reductase crystal sturctures
    if template1_PdbCode in ['3zgy', '3zhb', '4d3f', '4d3d', '4d3s', '4oqy', '4oqz', '4gmg', '4gmf', '4oqz', '5a9s',
                             '5a9r', '5g6r', '5g6s', '7xr5', '5fwn', '5a9t', '5kvq', '5kvs', '5ocm', '5ojl', '6eod',
                             '6eoh', '6toe','6eoi', '6to4', '6jit', '6rqa', '6jiz', '6grl', '6smt', '6sle', '6h7p',
                             '6jiz','7a3w', '7wnn', '7wnw', '7wxe', '7xe8', '7og3', '7osn']:
        return seqcode
    else:
        return ''