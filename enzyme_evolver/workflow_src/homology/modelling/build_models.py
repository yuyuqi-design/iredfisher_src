from modeller import *
from modeller.automodel import *
import os


def build_model(code, template1):
    """
    This function is to build a model based on a single template
    @param code: sequence code e.g. test
    @param template1: the name of the template protein, e.g. 5g6rA
    @return: none
    """

    env = environ()
    code = code.strip()
    a = automodel(env, alnfile=code + '-' + template1 + '.ali',
                  knowns=template1, sequence=code, assess_methods=(assess.DOPE, assess.GA341))
    a.starting_model = 1
    a.ending_model = 1
    try:
        a.make()
    except:
        pass
    # remove files that does not need
    try:
        os.remove(f'{code}.ini')
        os.remove(f'{code}.rsr')
        os.remove(f'{code}.sch')
        os.remove(f'{code}.V99990001')
        os.remove(f'{code}.D00000001')
    except FileNotFoundError:
        pass
    # rename the output model
    os.rename(f'{code}.B99990001.pdb', f"{code}.pdb")


def build_multi_models(code, template1, template2, template3):
    """
    This function is to build a model based on the top 3 templates
    @param code: sequence code e.g. test
    @param template1: the name of the first template protein, e.g. 5g6rA
    @param template2: the name of the first template protein, e.g. 6eodF
    @param template3: the name of the first template protein, e.g. 5ojlA
    @return: none
    """
    env = environ()
    a = automodel(env, alnfile=str(code).strip() + '-' + 'multi' + '.ali',
                  knowns=(template1, template2, template3), sequence=str(code).strip(),
                  assess_methods=(assess.DOPE, assess.GA341))
    a.starting_model = 1
    a.ending_model = 1
    a.make()
    # remove files that does not need
    try:
        os.remove(f'{code}.rsr')
        os.remove(f'{code}.sch')
        os.remove(f'{code}.V99990001')
        os.remove(f'{code}.D00000001')
    except FileNotFoundError:
        pass
    # rename the output model
    os.rename(f'{code}.B99990001.pdb', f"{code}.pdb")


class MyModel(automodel):
    """
    this class is to create constraint for modelling of dimeric protein
    """
    def special_restraints(self, aln):
        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
        s1 = selection(self.chains['A']).only_atom_types('CA')
        s2 = selection(self.chains['B']).only_atom_types('CA')
        self.restraints.symmetry.append(symmetry(s1, s2, 1.0))

    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)


def build_dimer(code, template_dimer='dimerT.pdb'):
    """
    this funtion is to create homodimeric protein based on the alignment file named code+ 'dimer.ali'
    and a dimer template strucutre named dimerT.pdb
    @param code: sequence code e.g. test
    @param template_dimer: the structure file of the dimer template
    @return: none
    """
    template_dimer = template_dimer.split('.')[0]
    code = code.strip()
    env = environ()
    # directories for input atom files
    env.io.atom_files_directory = ['.', '../atom_files']

    # Be sure to use 'MyModel' rather than 'automodel' here!
    a = MyModel(env,
                alnfile=code+'dimer.ali',  # alignment filename
                knowns=template_dimer,  # codes of the templates
                sequence=str(code).strip())  # code of the target

    a.starting_model = 1  # index of the first model
    a.ending_model = 1  # index of the last model
    # (determines how many models to calculate)
    a.make()  # do comparative modeling
    # remove files that does not need
    try:
        os.remove(f'{code}.ini')
        os.remove(f'{code}.rsr')
        os.remove(f'{code}.sch')
        os.remove(f'{code}.V99990001')
        os.remove(f'{code}.D00000001')
    except FileNotFoundError:
        pass
    # rename the output model
    os.rename(f'{code}.B99990001.pdb', f"{code}.pdb")
