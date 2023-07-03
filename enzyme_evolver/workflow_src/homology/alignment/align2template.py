import modeller


def align2d(code, template1, template1_PdbCode, template1_chain):
    """
    for homology modelling using one template
    this function is from Modeller 9.23 to align the user's sequence to template sequence
    """
    code = code.strip()
    env = modeller.environ()
    aln = modeller.alignment(env)
    mdl = modeller.model(env, file=template1_PdbCode + '.pdb',
                         model_segment=('FIRST:' + template1_chain, 'LAST:' + template1_chain))
    # add the template sequence
    aln.append_model(mdl, align_codes=template1, atom_files=template1_PdbCode + '.pdb')
    # add the sequence to model
    aln.append(file=str(code).strip() + '.pir', align_codes=str(code).strip())
    # align the sequence to the template
    aln.align2d()
    # write out the alignment file
    aln.write(file=code + '-' + template1 + '.ali', alignment_format='PIR')
    aln.write(file=code + '-' + template1 + '.pap', alignment_format='PAP')

    # calculate the sequence identity  between the target sequence and the template protein
    env = modeller.environ()
    env.io.atom_files_directory = ['../atom_files']
    # Read all entries in this alignment:
    aln = modeller.alignment(env, file=code + '-' + template1 + '.ali')
    sequence_identity = int(aln[0].get_sequence_identity(aln[1]))
    return sequence_identity


def salign(seqcode, template1_PdbCode='', template1_chain='', template2_PdbCode='', template2_chain='',
           template3_PdbCode='', template3_chain=''):
    """
    for homology modelling using multiple templates
    this function is from Modeller 9.23 to create an intial alignment of the top 3 template proteins
    """

    env = modeller.environ()
    env.io.atom_files_directory = './:../atom_files/'

    aln = modeller.alignment(env)
    # read the sequences of template proteins in PDB format
    for (Code, chain) in \
            ((template1_PdbCode, template1_chain), (template2_PdbCode, template2_chain),
             (template3_PdbCode, template3_chain)):
        mdl = modeller.model(env, file=Code + '.pdb', model_segment=('FIRST:' + chain, 'LAST:' + chain))
        aln.append_model(mdl, atom_files=Code+'.pdb', align_codes=Code + chain)
    # align the sequences
    for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                        ((1., 0.5, 1., 1., 1., 0.), False, True),
                                        ((1., 1., 1., 1., 1., 0.), True, False)):
        aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50),
                   gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                   dendrogram_file='dedro.tree',
                   alignment_type='tree',  # If 'progresive', the tree is not
                   # computed and all structues will be
                   # aligned sequentially to the first
                   feature_weights=weights,  # For a multiple sequence alignment only
                   # the first feature needs to be non-zero
                   improve_alignment=True, fit=True, write_fit=write_fit,
                   write_whole_pdb=whole, output='ALIGNMENT QUALITY')
    # write out the alignment file
    aln.write(file=str(seqcode).strip() + 'salign.pap', alignment_format='PAP')
    aln.write(file=str(seqcode).strip() + 'salign.ali', alignment_format='PIR')
    # calculate the alignment quality score
    aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
               gap_gap_score=0, gap_residue_score=0, dendrogram_file='dedro.tree',
               alignment_type='progressive', feature_weights=[0] * 6,
               improve_alignment=False, fit=False, write_fit=True,
               write_whole_pdb=False, output='QUALITY')


def align2d_mult(code):
    """
        for homology modelling using multiple templates
        this function is from Modeller 9.23 to align input sequence to the aligned template sequences
    """
    code = str(code).strip()
    env = modeller.environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read aligned template proteins:
    aln = modeller.alignment(env)
    aln.append(file=code + 'salign' + '.ali', align_codes='all')
    aln_block = len(aln)

    # Read aligned sequence(s):
    aln.append(file=code + '.pir', align_codes=code)

    # Structure sensitive variable gap penalty sequence-sequence alignment:
    aln.salign(output='', max_gap_length=20,
               gap_function=True,  # to use structure-dependent gap penalty
               alignment_type='PAIRWISE', align_block=aln_block,
               feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
               gap_penalties_1d=(-450, 0),
               gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
               similarity_flag=True)
    # write out the alignment file containing input seqeunce and the sequences of template proteins
    aln.write(file=code + '-' + 'multi' + '.ali', alignment_format='PIR')
    aln.write(file=code + '-' + 'multi' + '.pap', alignment_format='PAP')


def align_dimer(code, start_residue, end_residue, N_chains='2',
                template_monomer='singleT.pdb'):
    """
    for homology modelling of sequences in homodimer form based on a give template structure file
    code: the input sequence code, for example, test.fasta would be test as "code"
    singleT.pdb: monomer form of the template structure with chain labelled as A
    dimerT.pdb: homodimeric form of the template structure with chain labelled as A  and B
    start_residue: the residue index of the first amino acid in singleT.pdb, for example, 9
    end_residue: the residue index of the last amino acid in singeT.pdb, for example, 294
    N_chains: number of chains. the default is 2 chains.
    """
    code = code.strip()
    # align the target sequence to the monomeric template and output the alignment file: code-singleTA.ali
    template1 = template_monomer.split('.')[0] + 'A'
    template1_PdbCode = template_monomer.split('.')[0]
    template1_chain = 'A'
    align2d(code, template1, template1_PdbCode, template1_chain)

    # Read the monomeric alignmnet file
    with open(code + '-singleTA.ali') as single:
        # write out a new dimeric alignment file by hacking the monomeric alignment file
        with open(code+'dimer.ali', 'w') as multi:
            ali = single.readlines()
            for line in ali:
                # get the line starting with structure to locate the template sequence
                if line.startswith('structure'):
                    # get the line number of the template protein start
                    line_struc_BEGIN = ali.index(line) + 1
                # get the line starting with sequence to locate the target sequence
                if line.startswith('sequence'):
                    # get the line number of the target sequence
                    line_sequence_BEGIN = ali.index(line) + 1
                    # get the line number of the template protein end
                    line_structure_END = line_sequence_BEGIN - 2
                    # get the line number where the file end
                    line_file_END = len(ali)
            # the whole aligned sequence of the structure file
            structure_sequence = ''.join(ali[line_struc_BEGIN:line_structure_END]).replace("*", "").strip()
            # the whole aligned sequence of the target protein
            target_sequence = ''.join(ali[line_sequence_BEGIN: line_file_END]).replace("*", "").strip()
            # write the headline of the dimeric alignment file
            multi.writelines(ali[0])
            # change of name of the template structure from singleT to dimerT
            multi.write(ali[1].replace('singleTA', 'dimerT'))
            # write thhe details of the template structure file
            multi.write('structureX:dimerT:   ' + start_residue + ':A:' + end_residue + ':B:::-1.00:-1.00\n')
            # write two times of the whole aligned sequence of the structure file separated by / and ended by *
            multi.write((structure_sequence + '/' + '\n') * (int(N_chains) - 1))
            multi.write((structure_sequence + '*' + '\n'))
            # write the details of the target sequence
            multi.writelines(ali[line_sequence_BEGIN - 2:line_sequence_BEGIN])
            # write two times of the whole aligned sequence separated by / and ended by *
            multi.write((target_sequence + '/' + '\n') * (int(N_chains) - 1))
            multi.write((target_sequence + '*' + '\n'))
    # calculate the sequence identity  between the target sequence and the template protein
    env = modeller.environ()
    env.io.atom_files_directory = ['../atom_files']
    # Read all entries in this alignment:
    aln = modeller.alignment(env, file=f'{code}-singleTA.ali')
    sequence_identity = int(aln[0].get_sequence_identity(aln[1]))
    return sequence_identity
