from enzyme_evolver.workflow_src.homology.seqformat import fasta2pir
from enzyme_evolver.workflow_src.homology.find_template import create_alignment, find_template
from enzyme_evolver.workflow_src.homology.find_template import analyze_template
from enzyme_evolver.workflow_src.homology.alignment import align2template
from enzyme_evolver.workflow_src.homology.modelling import build_models
from enzyme_evolver.workflow_src.homology.evaluation import evaluation
import pandas as pd


def auto_single_modelling(seqcode='test', db='enzyme_evolver/database'):
    """
    this function is to create homology model for the target sequence
    using the best matched template structure in PDB database
    @param seqcode: test
    @param db: the pdb structure database for searching the template structures
    @return: df_best_template: a dataframe on the template information, sequence identity, ramachandran score etc.
    """
    # convert the sequence test.fasta to test.pir
    fasta2pir.fasta2pir(seqcode)
    # create alignment between the target sequence and the sequences in PDB database
    create_alignment.create_alignment(seqcode, db=db)
    # get the top 3 template structures of the target seqeunce
    templates, template1, template2, template3, template1_chain, template1_PdbCode, \
        template2_chain, template2_PdbCode, template3_chain, template3_PdbCode = find_template.find_template(seqcode)
    df_best_template = find_template.template_info(seqcode)
    if templates:
        # align the sequence to the best template structure
        align2template.align2d(seqcode, template1=template1, template1_PdbCode=template1_PdbCode,
                               template1_chain=template1_chain)
        # build models
        build_models.build_model(seqcode, template1=template1)
    else:
        pass
    # use procheck to evaluate the modelled structure
    evaluation.procheck(seqcode)
    # get the ramachandran core score of the structure
    ramachandran_score = evaluation.rama_core(seqcode)
    # add the ramachandran score to the spreadsheet containing template and model information
    df_best_template['ramachandran_score'] = ramachandran_score
    return df_best_template


def auto_multiple_modelling(seqcode, db='enzyme_evolver/database'):
    """
    this function is to create homology model for the target sequence
    using the top 3 template structure in PDB database
    @param seqcode: the sequence code e.g. test
    @param db: the pdb structure database for searching the template structures
    @return: df_best_template: a dataframe on the template information, sequence identity, ramachandran score etc.
    """
    # convert the sequence test.fasta to test.pir
    fasta2pir.fasta2pir(seqcode)
    # create alignment between the target sequence and the sequences in PDB database
    create_alignment.create_alignment(seqcode, db=db)

    try:
        # find the top 3 template
        templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, \
            template2_PdbCode, template3_chain, template3_PdbCode = find_template.find_template(seqcode)
        if templates and template1 and template2 and template3:
            align2template.salign(seqcode, template1_PdbCode=template1_PdbCode, template1_chain=template1_chain,
                                  template2_PdbCode=template2_PdbCode, template2_chain=template2_chain,
                                  template3_PdbCode=template3_PdbCode, template3_chain=template3_chain)
            align2template.align2d_mult(seqcode)
            # auto_modeller.trim_TerGap(str(seqcode).strip() + '-' + 'multi' + '.ali')
            build_models.build_multi_models(seqcode, template1=template1, template2=template2, template3=template3)
        else:
            pass
    except:
        pass
    # Get the information of the best template
    df_best_template = find_template.template_info(seqcode)
    # use procheck to evaluate the modelled structure
    evaluation.procheck(seqcode)
    # get the ramachandran core score of the structure
    ramachandran_score = evaluation.rama_core(seqcode)
    df_best_template['ramachandran_score'] = ramachandran_score

    return df_best_template


def template_single_modelling(seqcode, template='5g6r.pdb'):
    """
    this function is to create homology model for the target sequence
    using a given template struture
    @param seqcode: test
    @param template: the template structure file
    @return df_best_template: a dataframe on the template information, sequence identity, ramachandran score etc.
    """
    # convert the sequence test.fasta to test.pir
    fasta2pir.fasta2pir(seqcode)
    # get the chain name of the template
    chain_name = analyze_template.chain_name_template(file=template)
    # template
    template1 = template.split('.')[0] + chain_name
    # template code
    template1_PdbCode = template.split('.')[0]
    # template chain
    template1_chain = chain_name
    # align the target sequence to the template protein
    sequence_identity = align2template.align2d(seqcode, template1, template1_PdbCode, template1_chain)
    # build model
    build_models.build_model(seqcode, template1)

    # use procheck to evaluate the modelled structure
    evaluation.procheck(seqcode)
    # get the ramachandran core score of the structure
    ramachandran_score = evaluation.rama_core(seqcode)
    # collect the information about the sequence identity and model quality
    l_seqcode = []
    l_seqcode.append(seqcode)
    l_sequence_identity = []
    l_sequence_identity.append(sequence_identity)
    l_ramachandran_score = []
    l_ramachandran_score.append(ramachandran_score)
    data = {'sequence': l_seqcode, 'sequence_identity': l_sequence_identity,
            'ramachandran_score': l_ramachandran_score}

    df_best_template = pd.DataFrame(data=data)
    return df_best_template


def homodimer_modelling(seqcode, template_monomer='singleT.pdb', template_dimer='dimerT.pdb',
                        starting_residue='9', end_residue='294'):
    """
    this funtion is to create dimeric model for the target sequence
    using a given template struture
    @param seqcode:
    @param template_dimer: the dimer template file: e.g. dimerT.pdb
    @param template_monomer: the monomer template file: e.g. singleT.pdb
    @param starting_residue: the first residue number in monomer template file: e.g. 9
    @param end_residue:  the last residue number in monomer template file: e.g. 294
    @return: a dataframe on the template information, sequence identity, ramachandran score etc.
    """
    # convert the sequence test.fasta to test.pir
    fasta2pir.fasta2pir(seqcode)
    # create the dimeric alignment
    sequence_identity = align2template.align_dimer(seqcode, start_residue=starting_residue,
                                                   end_residue=end_residue, template_monomer=template_monomer)
    # build the model
    build_models.build_dimer(seqcode, template_dimer=template_dimer)
    # use procheck to evaluate the modelled structure
    evaluation.procheck(seqcode)
    # get the ramachandran core score of the structure
    ramachandran_score = evaluation.rama_core(seqcode)
    # collect the information about the sequence identity and model quality
    l_seqcode = []
    l_seqcode.append(seqcode)
    l_sequence_identity = []
    l_sequence_identity.append(sequence_identity)
    l_ramachandran_score = []
    l_ramachandran_score.append(ramachandran_score)
    data = {'sequence': l_seqcode, 'sequence_identity': l_sequence_identity,
            'ramachandran_score': l_ramachandran_score}

    df_best_template = pd.DataFrame(data=data)

    return df_best_template
