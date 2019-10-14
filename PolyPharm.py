#!/usr/bin/python

'''
author:     Abhijeet Kapoor
contact:    kapoor.abhijeet@gmail.com
year: 2015
'''

import sys
import logging
import getpass
from optparse import OptionParser
import numpy as np
import datetime
import subprocess
import time 
from FilterSequence import Filter_Sequence
from SimilarityBLAST import Similarity_Search
from TargetLigandBindingDB import Target_Ligand_BindingDB
from CombineRedundant import Combine_Redundant
from FindCompounds import Find_Compounds
from FindTargets import Find_Targets

def Display_Menu():
    print "1. %s\n2. %s\n3. %s\n4. %s\n5. %s\n6. %s\n7. %s\n8. %s\n" % (option1, option2, option3, option4, option5, option6, option7, option8)
    answer= raw_input("What would you like to do (select one):")
    return answer

def Polypharm():
    #Get the Start time to time code
    start_time = time.clock()
    #Report the date, time and username#
    now=datetime.datetime.now()
    date=now.strftime("%B %d, %Y %H:%M:%S")
    logging.info("Run on: %s", date)
    logging.info("By user: %s", getpass.getuser())

    #######Display Welcome Screen and read user input#######
    display=True
    print("\n\n########## Welcome to PolyPharmacology Pipeline ##########\n\nList of things that you can do:\n ")
    while display:
	user_selection=Display_Menu()
    	if user_selection == "1":
		display=False
		logging.info("User Selected option 1: %s", option1)
		print info_option1, "\n"

		reference_sequence= raw_input("Path to reference sequence file:")
		target_identifier_csv = raw_input("Path to target identifier csv file:")

		logging.info("\tUser input reference sequence file: %s", reference_sequence)
		logging.info("\tUser input target identifier file: %s", target_identifier_csv)
		Filter_Sequence(reference_sequence,target_identifier_csv)

	elif user_selection == "2":
		display=False
		logging.info("User Selected option 2: %s", option2)
                print info_option2, "\n"

                database_sequence= raw_input("Path to FASTA formatted sequence file:")
                database_type = raw_input("Database Type (nucl or prot):")
                database_name = raw_input("Database Name:")

                logging.info("\tUser input FASTA sequence file: %s", database_sequence)
                logging.info("\tUser input database type: %s", database_type)
                logging.info("\tUser input database name: %s", database_name)

                proc=subprocess.Popen(['makeblastdb', '-in', database_sequence, '-out', database_name, '-dbtype', database_type, '-parse_seqids'],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (stdout_data, stderr_data)= proc.communicate()
                logging.info("BLAST+ command makeblastdb output: %s\n", stdout_data)
                if stderr_data != '':
                        logging.error("***BLAST+ command makeblastdb ERROR***\n%s", stderr_data)  
                	sys.exit("***ERROR creating database***: Check log_runinfo")
                else:
                	print "Success: Database Created!!"
        elif user_selection == "3":
                display=False
                logging.info("User Selected option 3: %s", option3)
                print info_option3, "\n"

                blast_type= raw_input("Which BLAST to use (blastp or blastn..etc.):")
                query_sequence= raw_input("Query sequence file in FASTA format:")
                database_path = raw_input("Path to local BLAST Database:")
                database_type = raw_input("Database Type (nucl or prot):")
                eval_cutoff = raw_input("Evalue threshold:")
                output_filename = raw_input("Output filename:")

                logging.info("\tUser input BLAST type: %s", blast_type)
                logging.info("\tUser input query sequence file: %s", query_sequence)
                logging.info("\tUser input database path: %s", database_path)
                logging.info("\tUser input database type: %s", database_type)
                logging.info("\tUser input Evalue cutoff: %s", eval_cutoff)
                logging.info("\tUser input output filename: %s", output_filename)

                Similarity_Search(blast_type,query_sequence,database_path,eval_cutoff,output_filename,database_type)

        elif user_selection == "4":
                display=False
                logging.info("User Selected option 4: %s", option4)
                print info_option4, "\n"

                bindingDB_file_path= raw_input("Path to BindingDB ligand-target affinity dataset file:")
                ic50_cutoff= raw_input("IC50 threshold (nM) to filter ligands:")
                ki_cutoff = raw_input("Ki threshold (nM) to filter ligands:")
                kd_cutoff = raw_input("Kd threshold (nM) to filter ligands:")
                ec50_cutoff = raw_input("EC50 threshold (nM) to filter ligands:")
                output_filename = raw_input("Output filename:")
		target_filename = raw_input("Path to tsv file with required targets name (optional):")
		database_path = raw_input("Path to local BLAST Database (optional):")
                database_type = raw_input("Database Type [nucl or prot](required if BLAST database passed):")

		try:
			if database_path.strip() != '':
				assert target_filename.strip() != ''
		except AssertionError:
			sys.exit("***ERROR: Path to tsv file with required targets name is required when passing local BLAST Database***")
			

                logging.info("\tUser input path to BindingDB ligand-target affinity dataset file: %s", bindingDB_file_path)
                logging.info("\tUser input IC50 threshold to filter ligands: %s", ic50_cutoff+" nM")
                logging.info("\tUser input Ki threshold to filter ligands: %s", ki_cutoff+" nM")
                logging.info("\tUser input Kd threshold to filter ligands: %s", kd_cutoff+" nM")
                logging.info("\tUser input EC50 threshold to filter ligands: %s", ec50_cutoff+" nM")
                logging.info("\tUser input output filename: %s", output_filename)
		logging.info("\tUser input Path to tsv file with required targets name (optional): %s", target_filename)
		logging.info("\tUser input database path: %s", database_path)
                logging.info("\tUser input database type: %s", database_type)

                Target_Ligand_BindingDB(bindingDB_file_path,ic50_cutoff,ki_cutoff,kd_cutoff,ec50_cutoff,output_filename,target_filename,database_path,database_type)

	elif user_selection == "5":
                display=False
                logging.info("User Selected option 5: %s", option5)
                print info_option5, "\n"

                blast_type= raw_input("Which BLAST to use (blastp or blastn..etc.):")
                query_sequence= raw_input("Query sequence file in FASTA format:")
                database_path = raw_input("Path to local BLAST Database:")
		target_filename = raw_input("Path to tSV file with target id and names:")
                t_l_map = raw_input("Path to target-ligand map tsv file:")
                eval_cutoff = raw_input("E-value threshold (leave blank for using identity and coverage):")
		p_identity_cutoff = raw_input("Percent identity cutoff (only if Evalue was not specified):")
		p_coverage_cutoff = raw_input("Percent hit (or query) coverage cutoff (only if Evalue was not specified):")
		ligand_cutoff = raw_input("Skip target ligand map if number of ligands is less than:")
                output_filename = raw_input("Output filename:")

		try:
                        if eval_cutoff.strip() == '':
                                assert p_identity_cutoff.strip() != '' and p_coverage_cutoff.strip() != ''
                except AssertionError:
                        sys.exit("***ERROR: Percent identity and hit (or query) coverage required when not specifying E-value***")

		try:
                        if eval_cutoff.strip() != '':
                                assert p_identity_cutoff.strip() == '' and p_coverage_cutoff.strip() == ''
                except AssertionError:
                        print "***WARNING: Percent identity and hit (or query) coverage is not used when E-value is provided***"

                logging.info("\tUser input BLAST type: %s", blast_type)
                logging.info("\tUser input query sequence file: %s", query_sequence)
                logging.info("\tUser input database path: %s", database_path)
		logging.info("\tUser input path to tSV file with target id and names: %s", target_filename)
                logging.info("\tUser input path to target-ligand map tsv file: %s", t_l_map)
                logging.info("\tUser input Evalue cutoff: %s", eval_cutoff)
		logging.info("\tUser input percent identity cutoff: %s", p_identity_cutoff)
		logging.info("\tUser input percent coverage cutoff: %s", p_coverage_cutoff)
		logging.info("\tUser input skip target ligand map if number of ligands is less than: %s", ligand_cutoff)
                logging.info("\tUser input output filename: %s", output_filename)

		Combine_Redundant(blast_type,query_sequence,database_path,target_filename,t_l_map,eval_cutoff.strip(),p_identity_cutoff.strip(),p_coverage_cutoff.strip(),int(ligand_cutoff),output_filename)
       
	elif user_selection == "6":
                display=False
                logging.info("User Selected option 6: %s", option6)
                print info_option6, "\n"

                target_file_path= raw_input("Path to target info tsv file:")
                smiles_file_path= raw_input("Path to directory containing ligand SMILES (smi) files:")
                input_targets = raw_input("Comma separated list of target identifier numbers:")
                fp_type = raw_input("Fingerprint type (options: fp2, fp3, fp4, maccs):")
                tc_cutoff = raw_input("Tanimoto cutoff:")
		not_in_targets = raw_input("Similar ligands not in targets (comma separated list; Optional):")

                logging.info("\tUser input path to target info tsv file: %s", target_file_path)
                logging.info("\tUser input path to directory containing ligand SMILES (smi) files: %s", smiles_file_path)
                logging.info("\tUser input comma separated list of target identifier numbers: %s", input_targets)
                logging.info("\tUser input fingerprint type (options: fp2, fp3, fp4, maccs): %s", fp_type)
                logging.info("\tUser input tanimoto cutoff: %s", tc_cutoff)
		logging.info("\tUser input similar ligands not in targets (comma separated list; Optional): %s", not_in_targets)

		Find_Compounds(target_file_path,smiles_file_path,input_targets,fp_type.strip().lower(),float(tc_cutoff),not_in_targets)

	elif user_selection == "7":
                display=False
                logging.info("User Selected option 7: %s", option7)
                print info_option7, "\n"

                target_file_path= raw_input("Path to target info tsv file:")
                smiles_file_path= raw_input("Path to directory containing ligand SMILES (smi) files:")
                query_smi_file = raw_input("Query SMILE file (One SMILE per line):")
                fp_type = raw_input("Fingerprint type (options: linear or circular):")
		if fp_type.strip().lower() == 'linear':
			package_type = raw_input("Use linear fingerprint from (options: openbabel or rdkit):")
		elif fp_type.strip().lower() == 'circular':
			package_type = 'rdkit'
		else:
			sys.exit("***ERROR: Invalid fingerprint type***")
		tc_cutoff = raw_input("Tanimoto cutoff:")
		K_KNN = raw_input("Input K (integer) for KNN scoring:")
                
                logging.info("\tUser input path to target info tsv file: %s", target_file_path)
                logging.info("\tUser input path to directory containing ligand SMILES (smi) files: %s", smiles_file_path)
                logging.info("\tUser input query SMILE file (One SMILE per line): %s", query_smi_file)
                logging.info("\tUser input fingerprint type (options: linear or circular): %s", fp_type)
		if fp_type.strip().lower() == 'linear':
			logging.info("\tUser input use linear fingerprint from (options: openbabel or rdkit): %s", package_type)
		else:
			logging.info("\tPackage to use was set to: %s", package_type)
                logging.info("\tUser input tanimoto cutoff: %s", tc_cutoff)
		logging.info("\tUser input K (integer) for KNN scoring: %s", K_KNN)

		if tc_cutoff.strip() == '':
			print 'No Tanimoto cutoff provided. A default of 0.6 will be used.\n'
			tc_cutoff = 0.6
			logging.info("Tanimoto cutoff was set to a default value: 0.6\n")

		if K_KNN.strip() == '':
                        K_KNN = 3
                        logging.info("K for KNN scoring was set to a default value: 3\n")
                
                Find_Targets(target_file_path,smiles_file_path,query_smi_file,fp_type.strip().lower(),package_type.strip().lower(),float(tc_cutoff), int(K_KNN))
		logging.info("Targets Fishing Complete!!\n")
	
        elif user_selection == "8":
		display=False
		print("Goodbye!")
	else:
		print("Invalid option! Try again.\n")

    end_time = time.clock()
    logging.info("\tTotal run time: %0.3f seconds process time" % (end_time-start_time,))

def main(argv):
    """Parse arguments."""
    description = "Python script to perform individual steps of PolyPharmacology pipeline."
    usage = "%prog [options]"
    version = "%prog: version 1.0 - created by Abhijeet Kapoor October 17, 2015"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    #parser.set_defaults(pdbfile=None, Lmethod='single', cophen=0.1, FpType='fp2', tanscore=2)
    options, args = parser.parse_args(args=argv[1:])

    #######Configuration for logging#######
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s: %(message)s', filename='log_runinfo', filemode='w')
    
    ######Call to Polypharm function######
    flag = Polypharm()
    return 0


######Global Variables######
option1 = "Filter Drugbank Target Sequence"
option2 = "Create a local BLAST database"
option3 = "Similarity search against local BLAST database"
option4 = "Create Target-Ligand-Map from BindingDB ligand-target affinity dataset"
option5 = "Combine reduntant target-ligand map entries"
option6 = "Find compounds for my target(s)"
option7 = "Find targets for my compound(s)"
option8 = "Exit"

info_option1 = "\nInput: \n1) Drugbank target sequence file in FASTA format (reference sequence)\n2) Drug target identifier CSV file from drugbank.\n\
\nDescription: This option creates a new FASTA sequence file containing reference sequences corresponding to drug target identifier in input csv file."

info_option2 = "\nInput: \n1) Sequence file in FASTA format for which database is created.\n2) Database type (nucl or prot).\n3) Output name of database\
\nDescription: This option creates a local BLAST database from the input FASTA sequence file."

info_option3 = "\nInput: \n1) Blast type (blastp, blastn, blastx etc.).\n2) Query sequence file in FASTA format.\n3) Path to local BLAST database.\n4) Database type (nucl or prot).\n\
5) E-value threshold.\n6) Output filename.\
\nDescription: This option performs BLAST similarity search using user specified blast program and query sequence file against a local BLAST database. The hit sequences (from the local database)\
 are saved as separate file (hits_output_filename.fasta) and a mapping is created between query sequences and corresponding hits (hit_query_map_output_filename.tsv)."

info_option4 = "\nInput: \n1) Path to BindingDB ligand-target affinity dataset tsv file.\n2) IC50 threshold (nM) to filter ligands.\n3) Ki threshold (nM) to filter ligands.\n4) Kd threshold (nM) to filter ligands.\n5) EC50 threshold (nM) to filter ligands.\n6) Output filename.\n7) Path to tsv file with required targets name (optional).\n8) Path to local BLAST database (optional).\n9) Database type (nucl or prot).\
\nDescription: This option screens the input BindingDB ligand-target affinity dataset file (in tsv format) and for each target in file creates a list of all ligands that binds to the target (output ligand smiles) satisfying either one of IC50, EC50, Ki, or Kd threshold in nanomolar units (i.e. activity < threshold). Additionally, target-ligand map is also saved for all targets provided as a separate input file (optional). The file must contain tab-separated values with target ids in 1st column and name in 2nd column, and one target entry per line. If a local database is also passed along with required target file then the sequences of required targets mapped to ligands are also saved. When passing local database, required target file is must."

info_option5 = "\nInput: \n1) Blast type (blastp, blastn, blastx etc.).\n2) Query sequence file in FASTA format.\n3) Path to local BLAST database.\n4) Path to tSV file with target id and names.\n\
5) Path to target-ligand map tsv file.\n6) E-value threshold (leave blank for using identity and coverage).\n7) Percent identity cutoff (only if Evalue was not specified).\n8) Percent hit (or query) coverage cutoff (only if Evalue was not specified).\n9) Ligand number cutoff to skip target-ligand map entry from output.\n10) Output filename.\
\nDescription: The target-ligand map created in Option 4 may have identical targets (but with different names; as in mutants) or highly homologous targets as different target-ligand entries. This option combines such redundant or homologous entries and creats a unique target-ligand map. BLAST similarity search is performed to identify reduandant/homologous targets using user specified blast program and query sequence file against a local BLAST database. If Evalue specified, than only Evalue threshold is used to identify homologous targets. If Evalue left blank then percent identity and hit (or query) coverage is used. The second option is more restrictive as setting identity and coverage to close to 100% will result in combining those targets only that are same but have different names."

info_option6 = "\nInput: \n1) Path to target info tsv file.\n2) Path to directory containing ligand SMILES (smi) files.\n3) Comma separated list of target identifier numbers.\n\
4) Fingerprint type (options: fp2, fp3, fp4, maccs).\n5) Tanimoto cutoff.\n6) Similar ligands not in targets (comma separated list; Optional).\
\nDescription: This option calculates chemical similarity (using input fingerprint type) between all the ligands of input targets and finds compounds that could potentially be active against all input targets based on the similarlity principle. Two molecules are defined similar if Tanimoto coefficient >= Tanimoto cutoff. Optionally, the identified similar ligands are further screened out if an exact match is found with ligands of targets in not_in_targets list."

info_option7 = "\nInput: \n1) Path to target info tsv file.\n2) Path to directory containing ligand SMILES (smi) files.\n3) Input query SMILE file.\n\
4) Fingerprint type (options: linear or circular). \n5) Use linear fingerprint from (options: openbabel or rdkit). \n6) Tanimoto cutoff.\n6) Input K (integer) for KNN scoring.\n\
\nDescription: This option calculates chemical similarity (using input fingerprint type) for each query SMILE (one SMILE per line) with all the ligands of targets in target info tsv file. Potential targets for a query compound are reported based on similarlity principle. Two molecules are defined similar if Tanimoto coefficient >= Tanimoto cutoff. Potential targets can be ranked according to one of the three scores assigned based on KNN scheme. The three scores include 1NN (Max Score), KNN (K is the input parameter), and Centroid Score (average similarity with all ligands of a target) to the query."

if __name__ == '__main__':
    sys.exit(main(sys.argv))

