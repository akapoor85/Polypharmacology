#!/usr/bin/python

import sys
import os
import logging
import numpy as np
import pandas as pds
import subprocess

def Target_Ligand_BindingDB(bindingDB_file_path,ic50_cutoff,ki_cutoff,kd_cutoff,ec50_cutoff,output_filename,target_filename,database_path,database_type):
    #Read following required data from bindingDB_file
    # Column 2: Ligand smiles
    # Column 5: Ligand (Monomer) ID
    # Column 7: Target name
    # Column 9: Ki (nM)
    # Column 10: IC50 (nM)
    # Column 11: Kd (nM)
    # Column 12: EC50 (nM)

    #bindingDB_data = np.genfromtxt(bindingDB_file_path, dtype=None, usecols=(1,6,8,9,10,11), skip_header=1, delimiter='\t') #dtype=None will guess the datatype of each column
    #Determine maximum row length (required otherwise panda will give error with variable row length, as it determines maximum row length by reading the name of columns and this can vary in bindingDB)
    ####Now Use pandas to read data as genfromtxt behaviour was not consistent; missing values are replaced by string nan
    #bindingDB_data_DF = pds.read_csv(bindingDB_file_path, dtype=None, usecols=(1,6,8,9,10,11), delimiter='\t', names=range(max_row_len), header=None, skiprows=1, engine='python') 
    #Even panda started failing with variable number of rows than column..passing names with after determining max_row_len also didn't fix the issue. Need to figure out how to use pandas in thi case.
    ###Convert panda dataframe to numpy array###
    #bindingDB_data= bindingDB_data_DF.as_matrix()

    file_input= open(bindingDB_file_path, 'r')
    firstline= file_input.readline()
    
    #max_row_len=0
    #Compile target-ligand map
    target_ligand_map={}
    count_line=1
    no_skipped=0  
    for line in file_input:
	line_split=line.split('\t')
	try:
		reactant_id=int(line_split[0])  #To escape inappropriate entries
	except ValueError:
		logging.warning("Inappropriate data entry. Skipped line number: %d. Reactant id: %s " % (count_line, line_split[0]))
		no_skipped+=1
		continue
        ligand_smiles= line_split[1]
        ligand_id= line_split[4].strip()
	#Skip entry if ligand SMILE is empty
	if ligand_smiles.strip() == '':
		logging.warning("Missing ligand SMILE. Skipped line number: %d. Reactant id: %s. Ligand id: %s " % (count_line, line_split[0], ligand_id))
		continue
        #Quit if ligand_id is empty
        if ligand_id == '':
		sys.exit('Missing ligand id at line number: %d. Reactant id: %s' % (count_line, ligand_id))
        target_name= line_split[6].strip().lower() #Cases where target names are same but differ in letter case will be combined together
        if line_split[8].strip() == '':   #Ki value was missing
                ki_val=100000000.0 #Set to very large value,100 mM; these ligands are not of interest
        else:
		#Strip < or > symbol in the field
		line_split[8]=line_split[8].strip('>')
		line_split[8]=line_split[8].strip('<')
                ki_val=float(line_split[8])

        if line_split[9].strip() == '':   #IC50 value was missing
                ic50_val=100000000.0 #Set to very large value,100 mM; these ligands are not of interest
        else:
		#Strip < or > symbol in the field
                line_split[9]=line_split[9].strip('>')
                line_split[9]=line_split[9].strip('<')
		ic50_val=float(line_split[9])
		
        if line_split[10].strip() == '':   #Kd value was missing
                kd_val=100000000.0 #Set to very large value,100 mM; these ligands are not of interest
        else:
		#Strip < or > symbol in the field
                line_split[10]=line_split[10].strip('>')
                line_split[10]=line_split[10].strip('<')
                kd_val=float(line_split[10])

	if line_split[11].strip() == '':   #EC50 value was missing
                ec50_val=100000000.0 #Set to very large value,100 mM; these ligands are not of interest
        else:
		#Strip < or > symbol in the field
                line_split[11]=line_split[11].strip('>')
                line_split[11]=line_split[11].strip('<')
                ec50_val=float(line_split[11])

	#if ligand activity is within threshold then compile data for target 
        if ki_val < float(ki_cutoff) or ic50_val < float(ic50_cutoff) or kd_val < float(kd_cutoff) or ec50_val < float(ec50_cutoff):
                if not target_ligand_map.has_key(target_name):
                        target_ligand_map[target_name]= [[ligand_smiles, ligand_id, ki_val, ic50_val, kd_val, ec50_val]]
                else:
                        target_ligand_map[target_name].append([ligand_smiles, ligand_id, ki_val, ic50_val, kd_val, ec50_val])
    	count_line+=1
    
    #OUTPUT: Write the target-ligand mapping to a file
    file_map = open('target_ligand_map_'+output_filename+'.tsv', 'w')
    file_map.write('Target Name\tLigands Smile\tLigand ID\tKi (nM)\tIC50 (nM)\tKd (nM)\tEC50 (nM)\n')
    for key in target_ligand_map:
	file_map.write(key)
	for entries in target_ligand_map[key]:
		for columns in entries:
			file_map.write('\t'+str(columns))
		file_map.write('\n')
    if no_skipped>0:
	print "WARNING: %d entries skipped due to bad data format. Check log for more details" % (no_skipped,)	   
    	logging.warning("\nNumbers of line skipped: %d\n", no_skipped)
    logging.info("Target-ligand mapping written to file: %s\n", 'target_ligand_map_'+output_filename+'.tsv')
    logging.info("Number of Targets with ligand mapped: %d\n", len(target_ligand_map))

    #Compile results for required targets
    if target_filename.strip() != '':
	try:
		required_target_name= np.genfromtxt(target_filename, usecols=(0,1),delimiter='\t',dtype='string')
    	except IOError:
    		sys.exit("Cannot locate file %s" % target_filename + ".\n***Invalid filename or path***")

    	logging.info("\nNumber of target entries in user provided target file %s: %d\n" % (target_filename, len(required_target_name)))

    	#OUTPUT: Write the required target-ligand mapping to a file
	#Also write the target ids for which mapping is done to a separate file. This file can then be used to extract sequences from local blast database
    	file_required_target = open('Required_target_ligand_map_'+output_filename+'.tsv', 'w')
	file_mapped_target_ids = open('Required_target_mapped_ids_'+output_filename, 'w')
    	file_required_target.write('Target Name\tLigands Smile\tLigand ID\tKi (nM)\tIC50 (nM)\tKd (nM)\tEC50 (nM)\n')
	t_not_found=[]
	t_seen=[] #To avoid copying already written data if same target name repeats for different target_id 
	
        for t_id_name in required_target_name:
		try:
			if target_ligand_map[t_id_name[1].strip().lower()]:
				if t_id_name[1].strip().lower() in t_seen:
					#Write only the target id and continue to avoid repeating data
					file_mapped_target_ids.write(t_id_name[0]+'\n')
					continue
				else:
					t_seen.append(t_id_name[1].strip().lower())
				file_mapped_target_ids.write(t_id_name[0]+'\n')
				file_required_target.write(t_id_name[1])     #Change this to lower-case later.
				for entries in target_ligand_map[t_id_name[1].strip().lower()]:
					for columns in entries:
						file_required_target.write('\t'+str(columns))
					file_required_target.write('\n')
		except KeyError:
			t_not_found.append(t_id_name)

	logging.info("Required target-ligand mapping written to file: %s\n", 'Required_target_ligand_map_'+output_filename+'.tsv')
	logging.info("Number of required target-ligand maps written: %d\n" % (len(t_seen),))
	if len(t_not_found)>0:
		print "WARNING: Not all required targets were found. Check log for more details."
		logging.warning("\nNumber of targets in user provided target file for which no map is created: %d\n" % (len(t_not_found,)))
		logging.warning("\nList of target id and names not found:\n\tID\tName\n")
		for id_name in t_not_found:
			logging.warning(id_name[0]+"\t"+id_name[1])

	file_mapped_target_ids.close()
	#If local BLAST database passed then also extract sequences of required mapped targets
	if database_path.strip() != '':
		#Extract hit sequences to file hits_output_filename
		proc1=subprocess.Popen(['blastdbcmd', '-db', database_path, '-dbtype', database_type, '-entry_batch', 'Required_target_mapped_ids_'+output_filename ,'-outfmt','%f', '-out', 'Seq_Required_target_'+output_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		(stdout_data_p1, stderr_data_p1)= proc1.communicate()
		if stderr_data_p1 != '':
			logging.error("***BLAST hit Sequence retrieval ERROR***\n%s", stderr_data_p1)
			sys.exit("***ERROR running BLAST command blastdbcmd***: Check log_runinfo")
		else:
			logging.info("\nRequired target sequences for targets mapped to ligands written to file: %s\n", 'Seq_Required_target_'+output_filename)
		



