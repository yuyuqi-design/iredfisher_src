from pymol import cmd


def rec_align(recs):
    """
    align all the proteins to the first protein in the list
    @param recs: a list of names for protein structures, e.g. ['pro1', 'pro2', 'pro3']
    @return: None
    """
    #get the name of the first protein
    first_protein = recs[0].strip()
    #load proteins into pymol
    for protein in recs:
        protein = protein.strip() + '.pdb'
        cmd.load(protein)
    #align other proteins to the first protein
    cmd.alignto(first_protein)
    #save aligned proteins as name aln_protein.pdb
    for protein in recs:
        name = protein.strip().split('.')[0]
        cmd.save(name + '.pdb', name)
    #unload the proteins from pymol
    cmd.delete('all')
