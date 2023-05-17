import subprocess as sp


def run_vina(configure_file: str = 'config.txt', logfile: str = 'log.txt') -> None:
    """
    run docking by autodock vina
    @param configure_file: the configuration file for running vina
    @param logfile: the log file of running vina
    @return:
    """
    command = f'vina --config {configure_file} --out out.pdbqt --log {logfile}'
    try:
        sp.run(command, shell=True)
    except:
        pass
