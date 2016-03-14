#!/usr/bin/python

import sys
import os, glob
import logging
import numpy as np
import subprocess
from SimilarityBLAST import Run_BLAST
import matplotlib.pyplot as plt


def Id_Name(target_filename):
    id_names={} 
    file_id_name=np.genfromtxt(target_filename, dtype=None, delimiter = '\t')
    for line in file_id_name:
	t_id = line[0].strip()
	t_id_name = line[1].strip().lower()
	if not id_names.has_key(t_id):
		id_names[t_id]=t_id_name
	else:
		print "***ERROR: Multiple entry for target id: %s***" % (t_id)
		sys.exit("First Entry: id: %s, target name: %s\n" % (t_id, id_names[t_id])+ "Second Entry: id: %s, target name: %s\n" % (t_id, t_id_name))
    logging.info("\nTarget_ids->Target_names map created!!")
    logging.info("\nNumber of target ids mapped: %d\n" % (len(id_names),))
    return id_names

def Name_Ligand(t_l_map):
    names_ligand={}
    #file_name_ligand=np.genfromtxt(t_l_map, dtype=None, delimiter = '\t', usecols=(0,1,2), skip_header=1)
    file_name_ligand = open(t_l_map, 'r')
    firstline= file_name_ligand.readline()

    for line in file_name_ligand:
        line=line.split('\t')
        t_name = line[0].strip().lower()
        t_ligand_smile = line[1]
        ligand_id = int(line[2])
	if t_name.strip() != '':
		current_t_name= t_name
	elif t_name.strip() == '':
		t_name = current_t_name

        if not names_ligand.has_key(t_name):
                names_ligand[t_name]=[[t_ligand_smile], [ligand_id]]
        else:
                names_ligand[t_name][0].append(t_ligand_smile)
		names_ligand[t_name][1].append(ligand_id)
    logging.info("\nTarget_names->ligand map created!!")
    logging.info("\nNumber of target names mapped: %d\n" % (len(names_ligand),))
    return names_ligand

def Target_Target_Map(hit_map, T_id_to_name, eval_cutoff, p_identity_cutoff, p_coverage_cutoff):
    '''This function builds a prospective target-target combination list.'''
    t_t_map={}
    #Determine the combining criteria to use
    if eval_cutoff != '':
        use_iden = False
    elif p_identity_cutoff != '' and p_coverage_cutoff != '':
	use_iden = True

    #For each hit in hit map, find name of all targets to combine
    for h_id in hit_map:
	h_id_name = T_id_to_name[h_id]
	#logging.info("\nCompiling prospective target list for: %s\n", h_id_name)
        #Go through all matching targets for h_id
	for t_matched in hit_map[h_id]:
		t_m_id = t_matched[0] #Query id
		#Skip this entry if h_id and t_m_id are same
		if h_id == t_m_id:
			continue
		t_m_name = T_id_to_name[t_m_id]
		#Skip this entry if h_id and t_m_id both point to same target name
		if h_id_name == t_m_name:
                        continue
		#If names are different then create an entry for this target provided it satisfies the identity and coverage cutoff (if set)
		if use_iden:
			p_iden = float(t_matched[1]) #percent identity
			q_len = int(t_matched[3])  #query length
			h_len = int(t_matched[4])  #hit length
			#Determine whether q_coverage or h_coverage to use
			if q_len < h_len:
				u_cov = float(t_matched[8])  #query coverage
			else:
				u_cov = float(t_matched[7])  #hit coverage
			#logging.info("\nMatched target name: %s, p_iden: %s, u_cov: %s\n" % (t_m_name,str(p_iden),str(u_cov)))
			#Skip this entry if identity or coverage criteria are not satisfied
			if p_iden < float(p_identity_cutoff) or u_cov < float(p_coverage_cutoff):
				continue
			
		if not t_t_map.has_key(h_id_name):
			t_t_map[h_id_name]=[t_m_name]
		else: #If the entry for h_id already exists then add new t_m_name only if not previously added
			if t_m_name not in t_t_map[h_id_name]:
				#Before appending make sure that a previous h_id (for target with same name) did not result in 'Matches itself entry', if so then pop
				if t_t_map[h_id_name][0] == 'Matches Itself Only':
					pop_item = t_t_map[h_id_name].pop(0)   #This will pop entry at index 0
					logging.info("\n%s entry popped for target: %s\n" % (pop_item, h_id_name))
				t_t_map[h_id_name].append(t_m_name)
	#Before moving on to next hit make sure that h_id is not missing from t_t_map. This would happen if all t_matched entries of h_id correspond to same target names causing all entries to be skipped
	if not t_t_map.has_key(h_id_name):
		t_t_map[h_id_name]=['Matches Itself Only'] 
    return t_t_map

def Unique_Target_Map(targets_to_combine):
    u_t_l_map={}
    skip_target=[] #Keeps track of which targets are already combined
    t_unamb=0
    t_itself=0
    for t_name in targets_to_combine:
	logging.info("\nCompiling unique target list for: %s\n", t_name)
	#Skip if this target has already been combined with a previous target
	if t_name in skip_target:
		logging.info("\nSKIPPING %s. This target has already been combined.\n", t_name)
		continue
	elif targets_to_combine[t_name][0] == 'Matches Itself Only':
		logging.info("\nDONE %s. It matches itself only.\n", t_name)
		t_itself+=1
		if not u_t_l_map.has_key(t_name):
			u_t_l_map[t_name]=['Matches Itself Only', 'unambiguous', 0, 0]
		else:
			sys.exit('Multiple entry for target: %s' % (t_name))
		continue
	#Make a list of all targets for t_name (including itself) that are to be combined
	t_list= targets_to_combine[t_name] + [t_name]
	#Update skip_target list with t_name
	skip_target.append(t_name)
	#Create a list that will hold the final list of unambigous targets that will be combined with t_name
	unamb_t_list=[]
	#For t_name go through all prospective targets to combine
	logging.info("\nGoing through all prospective targets of: %s\n", t_name)
	for combine_t_name in targets_to_combine[t_name]:
		#Skip if this target has already been combined with a previous target
		if combine_t_name in skip_target:
			logging.warning("\nSKIPPING %s. Either it was already combined with some other unambiguous target or was not combined and exists as individual entry.\n", combine_t_name)
			continue
		#Extract target list for combine_t_name (including itself)
		combine_t_list= targets_to_combine[combine_t_name] + [combine_t_name]
		#Now check if t_name and combine_t_name can be unambigously combined. This is defined as two targets for which prospective target lists are same
		if set(t_list) == set(combine_t_list):
			#Update skip_target list with combine_t_name
			skip_target.append(combine_t_name)
			#Update unamb_t_list with combine_t_name
			unamb_t_list.append(combine_t_name)
		else: 
			logging.warning("\nSKIPPING %s. This prospective target cannot be unambigously combined.\n", combine_t_name)
                        
	#Check if unamb_t_list is same as targets_to_combine[t_name]
	if set(targets_to_combine[t_name]) == set(unamb_t_list):
		logging.info("\nDONE %s. Target list was unambiguous.\n", t_name)
		t_unamb+=1
		if not u_t_l_map.has_key(t_name):
			u_t_l_map[t_name]=[unamb_t_list, 'unambiguous', len(targets_to_combine[t_name]), len(set(targets_to_combine[t_name])-set(unamb_t_list))]  #Unambiguous tag here means all prospective target entries mapped each other
        	else:
			sys.exit('Multiple entry for target: %s' % (t_name))
	else:
		logging.warning("\nDONE %s. One or more targets in target list cannot be unambiguously combined.\n", t_name)
		logging.warning("\nDIFFERENCE: %s\n", set(targets_to_combine[t_name])-set(unamb_t_list))

		if not u_t_l_map.has_key(t_name):
			u_t_l_map[t_name]=[unamb_t_list, 'ambiguous', len(targets_to_combine[t_name]), len(set(targets_to_combine[t_name])-set(unamb_t_list))]
		else:
			sys.exit('Multiple entry for target: %s' % (t_name))
    logging.info("\nUnambiguously combined entries: %d\n", t_unamb)
    logging.info("\nMatching self only entries: %d\n", t_itself)
    logging.info("\nTotal entries for unique target ligand maps (excluding targets not found in BLAST output): %d\n", len( u_t_l_map))

    return u_t_l_map		

def Combine_ligands(unique_t_l_map, append_to_uni_t_l_map, T_name_to_ligand, lig_cutoff):
    count_targets=0
    skipped_target=0 #Keeps track of number of target-ligand map not satisfying lig_cutoff
    skipped_match_itself=0 #How many that are skipped match itself only
    skipped_missing_blast=0 #How many that are skipped are missing in BLAST output
    file_targets_info = open('Targets_for_ligands.tsv', 'w') #For each .smi file write their corresponding targets
    file_targets_info.write('Target_ID'+ '\t' + 'No_of_Ligands' + '\t' + 'Target_Name(s)')
    x_t_num=[] 
    y_l_num =[] #x, y are for plotting purpose
    #Create an Output Folder. 
    if not os.path.isdir("Target_Ligands"):
	os.mkdir("Target_Ligands")
    else:
	if os.listdir("Target_Ligands") != []:
		smifiles= glob.glob("Target_Ligands/*.smi")
		for filename in smifiles:
			os.remove(filename)

    #Go through all targets in unique_t_l_map and combine ligand entries removing all redundant entries of a ligand
    for u_target_name in unique_t_l_map:
	ligand_seen=[] #holds ligand id, For each target keeps track of unique ligand
	t_names=[] #List of target names to write corresponding to key count_targets
	ligand_smiles={} #Holds non_redundant ligand smiles for each target; Key: Ligand id 

	#append to t_names the current target name
	t_names.append(u_target_name)
	#Go through the ligand list for the target u_target_name and populate ligand_smiles with non-redundant smiles
	for l_s, l_id in zip(T_name_to_ligand[u_target_name][0], T_name_to_ligand[u_target_name][1]):
		if l_id in ligand_seen: #This entry already exits so skip
			continue
		else:
			ligand_seen.append(l_id)
			if not ligand_smiles.has_key(l_id):
				ligand_smiles[l_id]=l_s
			else:
				sys.exit('***ERROR: Multiple entry for ligand id: %d. Target: %s' % (l_id,u_target_name))
	#Now go through all the redundant targets of target u_target_name and combine those as well
	if unique_t_l_map[u_target_name][0] == 'Matches Itself Only':
		pass
	elif len(unique_t_l_map[u_target_name][0]) == 0:
		pass
	else:
		for t_redun_name in unique_t_l_map[u_target_name][0]:
			#append to t_names the current target name
        		t_names.append(t_redun_name)
			#Go through the ligand list for the target t_redun_name and populate ligand_smiles with non-redundant smiles
			for l_s, l_id in zip(T_name_to_ligand[t_redun_name][0], T_name_to_ligand[t_redun_name][1]):
				if l_id in ligand_seen: #This entry already exits so skip
					continue
				else:
					ligand_seen.append(l_id)
					if not ligand_smiles.has_key(l_id):
						ligand_smiles[l_id]=l_s
					else:
						sys.exit('***ERROR: Multiple entry for ligand id: %d. Target: %s' % (l_id,t_redun_name))
					
	#Now that all the target ligands are combined, write this list if number of ligand >= lig_cutoff or skip
	if len(ligand_smiles) >= lig_cutoff:
		count_targets+=1
		#Write the target info
		file_targets_info.write(str(count_targets)+ '\t' + str(len(ligand_smiles)))
		for names_t in t_names:
			file_targets_info.write('\t'+names_t)
		file_targets_info.write('\n')
	
		x_t_num.append(count_targets)
		y_l_num.append(len(ligand_smiles))	

		#Write Ligand smiles in Target_Ligands subdirectory
		file_ligands_smi = open('Target_Ligands/'+ str(count_targets)+'.smi', 'w')
		for key in ligand_smiles:
			file_ligands_smi.write(ligand_smiles[key] + " "+ str(key)+'\n')
		file_ligands_smi.close()
	else:
		skipped_target+=1
		if unique_t_l_map[u_target_name][0] == 'Matches Itself Only':
			skipped_match_itself+=1
			logging.warning("\nNumber of ligands (%d) is less than lig_cutoff (%d). Skipping target: %s (Matches Itself Only)\n" % (len(ligand_smiles), lig_cutoff, u_target_name))
		else:
			logging.warning("\nNumber of ligands (%d) is less than lig_cutoff (%d). Skipping target: %s\n" % (len(ligand_smiles), lig_cutoff, u_target_name))
		

    #Repeat the same for append_to_uni_t_l_map
    if len(append_to_uni_t_l_map) != 0:
	for u_target_name in append_to_uni_t_l_map:
        	ligand_seen=[] #holds ligand id, For each target keeps track of unique ligand
        	ligand_smiles={} #Holds non_redundant ligand smiles for each target; Key: Ligand id 

        	#Go through the ligand list for the target u_target_name and populate ligand_smiles with non-redundant smiles
        	for l_s, l_id in zip(T_name_to_ligand[u_target_name][0], T_name_to_ligand[u_target_name][1]):
                	if l_id in ligand_seen: #This entry already exits so skip
                        	continue
                	else:
                        	ligand_seen.append(l_id)
                        	if not ligand_smiles.has_key(l_id):
                                	ligand_smiles[l_id]=l_s
                        	else:
                                	sys.exit('***ERROR: Multiple entry for ligand id: %d. Target: %s' % (l_id,u_target_name))

		#Now that all the target ligands are combined, write this list if number of ligand >= lig_cutoff or skip
		if len(ligand_smiles) >= lig_cutoff:
			count_targets+=1
                	#Write the target info
			file_targets_info.write(str(count_targets)+'\t'+ str(len(ligand_smiles))+ '\t' +u_target_name+'\n')

			x_t_num.append(count_targets)
                	y_l_num.append(len(ligand_smiles))

                	#Write Ligand smiles in Target_Ligands subdirectory
                	file_ligands_smi = open('Target_Ligands/'+ str(count_targets)+'.smi', 'w')
                	for key in ligand_smiles:
                        	file_ligands_smi.write(ligand_smiles[key] + " "+ str(key)+'\n')
                	file_ligands_smi.close()
        	else:
			skipped_target+=1
			skipped_missing_blast+=1
                	logging.warning("\nNumber of ligands (%d) is less than lig_cutoff (%d). Skipping target: %s (Not in BLAST output)\n" % (len(ligand_smiles), lig_cutoff, u_target_name))	

    logging.info("\nNumber of target-ligands map skipped: %d\n" % (skipped_target,))
    logging.info("\n%d skipped targets are those that match itself only\n" % (skipped_match_itself,))    
    logging.info("\n%d skipped targets are those that were missing in BLAST output\n" % (skipped_missing_blast,))
    logging.info("\nTarget-Ligand Mapping Finished.\n\tLigand Smiles for each target are in subdirectory 'Target_Ligands'.\n\tCorresponding target info is in file 'Targets_for_ligands.tsv'") 
    logging.info("\nNumber of target-ligands map written in subdirectory 'Target_Ligands': %d\n" % (count_targets,))
    
    #Plot Target-Ligand map statistic
    plt.bar(x_t_num, y_l_num, alpha=0.4)
    plt.xlabel('Target Number')
    plt.ylabel('Number of Ligands')
    plt.savefig('Target-Ligand-Stat.png')
    plt.close('all')

    binW=100 #Bin-width
    plt.hist(y_l_num, bins=int((max(y_l_num)-min(y_l_num))/binW) ,color='b', alpha=0.5)
    plt.xlabel('Number of Ligands')
    plt.ylabel('Frequency')
    plt.title('Bin Width: '+ str(binW))
    plt.savefig('Target-Ligand-Stat-Hist.png')
    plt.close('all')
    


def Combine_Redundant(blast_type,query_sequence,database_path,target_filename,t_l_map,eval_cutoff,p_identity_cutoff,p_coverage_cutoff,lig_cutoff,output_filename):
    #Create a dictionary that maps target id to protein target names
    logging.info("\nCreating a map between target_ids and target_names using file: %s\n", target_filename)
    T_id_to_name = Id_Name(target_filename)

    #Create a dictionary that maps target name to associated ligand smiles
    logging.info("\nCreating a map between target_names and ligands using file: %s\n", t_l_map)
    T_name_to_ligand = Name_Ligand(t_l_map)

    #Call to Run_BLAST
    if eval_cutoff != '':
	Run_BLAST(blast_type,query_sequence,database_path,eval_cutoff,output_filename)
    elif p_identity_cutoff != '' and p_coverage_cutoff != '':
	default_eval_cutoff= '1e-120'  #Set it high enough to reduce unnecessary blast output
	logging.info("\nRunning BLAST similarity search with e-value cutoff set to: %s\n", default_eval_cutoff)
	Run_BLAST(blast_type,query_sequence,database_path,default_eval_cutoff,output_filename)
    else:
        sys.exit("***ERROR: Call to RunBLAST failed. One of eval_cutoff, p_identity_cutoff, or p_coverage_cutoff is invalid***")
    
    #Compile hit-query mapping
    hit_map={}
    output_data = np.genfromtxt(output_filename, dtype=None) #dtype=None will guess the datatype of each column
    logging.info("\nCompiling hit-query map from BLAST similarity search output\n")
    for line in output_data:
        query_id = line[0][line[0].find('|')+1:].strip() 
        hit_id = line[1]
        p_identity = line[2]
        a_length = line[3]
        query_length = line[4]
        hit_length = line[5]
        s_a_query = line[6]
        e_a_query = line[7]
        s_a_hit = line[8]
        e_a_hit = line[9]
        evalue= line[10]
        hit_coverage = '%.2f' % (100*(float((e_a_hit-s_a_hit+1))/hit_length))  #Not using a_length as it counts gaps also
        query_coverage = '%.2f' % (100*(float((e_a_query-s_a_query+1))/query_length))

        if not hit_map.has_key(hit_id):
                hit_map[hit_id]= [[query_id,p_identity,a_length,query_length,hit_length,s_a_hit,e_a_hit,hit_coverage,query_coverage,evalue]]
        else:
                hit_map[hit_id].append([query_id,p_identity,a_length,query_length,hit_length,s_a_hit,e_a_hit,hit_coverage,query_coverage,evalue])

    #OUTPUT:Write the hit-query mapping to a file
    file_map = open('hit_query_map_'+output_filename+'.tsv', 'w')
    file_map_90 = open('hit_query_map_'+output_filename+'_HitCoverage90.tsv', 'w')
    file_map.write('Hit_ID\tQuery_ID\tPercent_Identity\tAlignment_length\tQuery_length\tHit_length\tHit_alignment_start\tHit_alignment_end\tHit_coverage\tQuery_Coverage\tEvalue\n')
    file_map_90.write('Hit_ID\tQuery_ID\tPercent_Identity\tAlignment_length\tQuery_length\tHit_length\tHit_alignment_start\tHit_alignment_end\tHit_coverage\tQuery_Coverage\tEvalue\n')
    hits_90 =[]  #List of hits with coverage >=90
    for key in hit_map:
        write_key=True    #Keeps track of when to write key in file_map_90
        file_map.write(key)
        for entries in hit_map[key]:
                if float(entries[7]) >= 90.0:  #If hit coverage is >=90; entries[7] is entered as string
                        write_map_90=True
                else:
                        write_map_90=False
                for columns in entries:
                        file_map.write('\t'+str(columns))
                        if write_map_90:
                                if write_key:
                                        file_map_90.write(key)
                                        hits_90.append(key)
                                        write_key=False       #key is written for only first entry
                                file_map_90.write('\t'+str(columns))
                file_map.write('\n')
		if write_map_90:
                        file_map_90.write('\n')
    
    logging.info("\nHit query map written to file: %s\n" % ('hit_query_map_'+output_filename+'.tsv'))

    #Create a dictionary which keeps track of prospective targets that are combined together (for each hit in the hit-query map)
    logging.info("\nCreating a list of prospective targets to combine for each hit in the hit-query map\n")
    targets_to_combine = Target_Target_Map(hit_map, T_id_to_name, eval_cutoff, p_identity_cutoff, p_coverage_cutoff)
    logging.info("\nFinished creating prospective target list to combine\n")

    #Check if all targets from the user input are in the BLAST output (not all sequences will be in output if e-value is very high).
    #Make sure targets not in BLAST output are also included in final result (unique target ligand map)
    logging.info("\nNumber of hits (target_ids) in hit-query map: %d\n" % (len(hit_map),))
    logging.info("\nNumber of target names for which prospective list was created: %d\n" % (len(targets_to_combine),))
    append_to_uni_t_l_map = set()
    if len(targets_to_combine) < len(T_name_to_ligand):
	logging.info("\nNumber of targets in user provided target-ligand map file that are not present in BLAST output: %d\n" % (len(T_name_to_ligand)-len(targets_to_combine),))
        logging.info("\nThe following targets will not be combined to any other targets:\n")
	append_to_uni_t_l_map = set(T_name_to_ligand)-set(targets_to_combine) #These entries will later be part of unique target-ligand map
        for t_e_name in append_to_uni_t_l_map:
		logging.info(t_e_name)

    #OUTPUT: Write the prospective targets that are to be combined to a file
    file_t_combine = open('prospective_targets_combined_'+output_filename+'.tsv', 'w')
    for key in targets_to_combine:
	file_t_combine.write(key)
	for entries in targets_to_combine[key]:
		file_t_combine.write('\t'+entries)
	file_t_combine.write('\n')

    logging.info("\nProspective target list written to file: %s\n" % ('prospective_targets_combined_'+output_filename+'.tsv'))

    #Find unambiguous target combination 
    logging.info("\nCreating a unique target ligand map from prospective targets list.\n")
    unique_t_l_map = Unique_Target_Map(targets_to_combine)
    logging.info("\nFinished creating unique target ligand map.\n")  

    #OUTPUT: Write the final targets that are to be combined to a file
    file_ut_combine = open('Final_targets_combined_'+output_filename+'.tsv', 'w')
    file_ut_combine.write('Target_Name\tCombined_Unambiguously\tNo_of_targets_to_combine\tNumber_of_targets_not_combined\tCombined_With_Targets\n')
    for key in unique_t_l_map:
        file_ut_combine.write(key)
	if unique_t_l_map[key][0] == 'Matches Itself Only':
		file_ut_combine.write('\t1\tMatches Itself Only\t0\t0\n')
	elif len(unique_t_l_map[key][0]) == 0:
		file_ut_combine.write('\t'+str(0)+'\tAll targets ambiguous or already assigned unambiguously\t'+str(unique_t_l_map[key][2])+'\t'+str(unique_t_l_map[key][3])+'\n')
	else:
		if unique_t_l_map[key][1] == 'unambiguous':
			file_ut_combine.write('\t'+str(1)+'\t'+str(unique_t_l_map[key][2])+'\t'+str(unique_t_l_map[key][3]))
		elif unique_t_l_map[key][1] == 'ambiguous':
			file_ut_combine.write('\t'+str(0)+'\t'+str(unique_t_l_map[key][2])+'\t'+str(unique_t_l_map[key][3]))
		else:
			sys.exit('Expected ambiguity entry not found')
		for name in unique_t_l_map[key][0]:
			file_ut_combine.write('\t'+name)
		file_ut_combine.write('\n') 

    #Write to final targets all targets which were not in BLAST output
    if len(append_to_uni_t_l_map) != 0:
	logging.info("\nTotal entries for unique target ligand maps (including targets not found in BLAST output): %d\n" % (len(unique_t_l_map) + len(append_to_uni_t_l_map),)) 
	for key in append_to_uni_t_l_map:
		file_ut_combine.write(key)
		file_ut_combine.write('\t1\t0 (not in BLAST output)\t0\t0\n')

    logging.info("\nFinal unique target list written to file: %s\n" % ('Final_targets_combined_'+output_filename+'.tsv'))

    #Now Combine ligands based on unique target-ligand map
    logging.info("\nCombining ligand smiles based on unique target ligand map.\n")
    Combine_ligands(unique_t_l_map, append_to_uni_t_l_map, T_name_to_ligand, lig_cutoff)






    	
