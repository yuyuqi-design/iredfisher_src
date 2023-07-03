from modeller import environ
from modeller.scripts import complete_pdb


def chain_name_template(file):
    """
    this function is to read a pdb file and return the name of the first chain
    @param file: pdb file
    @return: the name of the first chain
    """
    env = environ()
    env.io.atom_files_directory = ['../atom_files']
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    mdl = complete_pdb(env, file)
    chain_name = [c.name for c in mdl.chains][0]

    return chain_name
