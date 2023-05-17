from enzyme_evolver.workflow_src import IREDFisher_ranking
from enzyme_evolver.workflow_src.docking import IREDFisher_docking, IREDFisher_scoring
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
from enzyme_evolver import job_functions
from pathlib import Path
from typing import List


class IREDFisher():
    """
    the class is to run IREDFisher and add output files to mongo database
    """
    def __init__(self, fasta_filename: str, lig_structures: List[str], N_ID: str, folder_id: str,
                 iredfisher_database='', print_log=True):
        """
        the input parameters to run IREDFisher
        @param fasta_filename: a sequence file containing multiple sequences
        @param lig_structures: ligand structures. e.g. ['lig1.pdb', 'lig2.pdb', ...]
        @param N_ID: the atom ID of the N atom in C=N bond in ligand structure
        @param folder_id: the folder name where the output files will be saved
        @param iredfisher_database: default is None, which enable the user upload and rank their own sequence files.
                other two options are 'Public' and 'Screened'. Option 'Public' enable the user to rank IRED sequences
                collected from public database and option 'Screened' enable the user to rank IRED sequences collected
                from published paper.
        @param print_log: default is true to print the log information
        """
        self.fasta_filename = fasta_filename
        self.lig_structures = lig_structures
        self.N_ID = N_ID
        self.folder_id = folder_id
        self.iredfisher_database = iredfisher_database
        self.print_log = print_log
        # base_folder: the directory where you execute the function
        self.base_folder = Path.cwd()
        # data_folder: the directory where the input data saved
        self.data_folder = f'{self.base_folder}/enzyme_evolver/database/{folder_id}'
        # IREDFisher database folder for saving all job files
        self.database_folder = f'{self.base_folder}/enzyme_evolver/database'
        # launch the mongo database service
        self.db_file = Job.objects(folder_id=folder_id)[0]
        # update computation time
        self.db_file.update_notes(f'please bookmark this page. Calculation of each sequence will'
                                  f' take about 1.5 minutes')


    def run_IREDFisher(self):
        """
        main function to run IREDFisher
        @return:
        """

        job_functions.update_a_job_status(job=self.db_file, status='Job running')
        # copy the template structures and reference ligand and cofactor files into the datafolder folder
        database_functions.copy_references(data_folder=self.data_folder)
        # run IREDFisher ranking
        rec_names = []
        # if the user upload a sequence file
        if self.fasta_filename:
            job_functions.update_a_job_status(job=self.db_file, status='Job running for ranking provided sequences')
            codes = IREDFisher_ranking.iredfisher_ranking(fasta_filename=self.fasta_filename,
                                                          lig_structures=self.lig_structures, N_ID=self.N_ID,
                                                          base_folder=self.base_folder, data_folder=self.data_folder)
            rec_names = rec_names + codes


        # if the user choose to rank sequences from iredfisher database where all protein strucutres are already built
        # and prepared for molecular docking
        else:
            job_functions.update_a_job_status(job=self.db_file, status=f'Job running for ranking '
                                                                       f'{self.iredfisher_database} sequences')
            # get the names of the proteins
            rec_names = database_functions.prepared_recs(iredfisher_database=self.iredfisher_database,
                                                     database_folder=self.database_folder, data_folder=self.data_folder)

            # docking the ligand structures to the proteins
            IREDFisher_docking.iredfisher_docking(rec_names=rec_names, lig_structures=self.lig_structures,
                                                  N_ID=self.N_ID, base_folder=self.base_folder,
                                                  data_folder=self.data_folder,
                                                  iredfisher_database=self.iredfisher_database)

            # run iredfisher scoring
            IREDFisher_scoring.iredfisher_scoring(rec_names=rec_names, lig_structures=self.lig_structures,
                                                  base_folder=self.base_folder, data_folder=self.data_folder)


        #add files to the mongo database
        job_functions.add_3dmodel_files(job=self.db_file, rec_names=rec_names)
        job_functions.add_binding_mode_structures(job=self.db_file, rec_names=rec_names,
                                                  lig_structures=self.lig_structures)
        job_functions.add_ranking_file(job=self.db_file, lig_structures=self.lig_structures)
        job_functions.add_zipfile(job=self.db_file, base_folder=self.base_folder, data_folder=self.data_folder)
        job_functions.update_a_job_status(job=self.db_file, status='Job finished')
