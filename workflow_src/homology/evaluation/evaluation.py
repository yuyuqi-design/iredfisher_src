import subprocess as sp


def procheck(code):
    """
    this function is to use software procheck to evaluate the modelled structure
    @param code: sequence code e.g. test
    @return: none
    """
    code = code.strip()
    evaluate_cmd = f"procheck.scr {code}.pdb 2.0"
    sp.run(evaluate_cmd, shell=True)
    # remove files
    rm_cmd = f"rm procheck.prm *.ps *lot.log *out *sdh *.lan *.pln *sco nb.log *.nb *.rin secstr.log *new anglen.log " \
             f"clean.log"
    try:
        sp.run(rm_cmd, shell=True)
    except FileNotFoundError:
        pass


def rama_core(code):
    """
        this function is to get the ramachandran core score of the structure. A value over 90% is considered good.
        @param code: sequence code e.g. test
        @return: core_score
    """
    with open(f'{code}.sum') as fin:
        result = fin.readlines()
        for line in result:
            if '| Ramachandran plot:' in line:
                core = float(line.split()[3].strip('%'))
                return core
