#!/usr/bin/python

from __future__ import division
import sys
import os, glob
import logging
import numpy as np
import pybel as pb
import matplotlib.pyplot as plt
from matplotlib import cm, colors


def Target_Ligands(t_iden_list, smiles_file_path):
    t_l_s={} #key is the target id
    for t_id in t_iden_list:
	try:
                t_smi_file = open(os.path.join(smiles_file_path, t_id+'.smi'), 'r')
        except:
                sys.exit("Cannot locate file %s" % os.path.join(smiles_file_path, t_id+'.smi') + ".\n***Invalid target identifier or path***")

        #Compile the list of ligand smiles for target t_id
        for line in t_smi_file:
                line=line.split()
                l_smiles= line[0]
                l_id = int(line[1])
                if not t_l_s.has_key(t_id):
                        t_l_s[t_id] = [[l_smiles], [l_id]]
                else:
                        t_l_s[t_id][0].append(l_smiles)
                        t_l_s[t_id][1].append(l_id)
        t_smi_file.close()
    return t_l_s


def CalcFpt(t_smi_list, FPType):
    fpts_list=[]
    l_in_fmt = 'smi' #Input format is SMILES
    for l_smi in t_smi_list:
	#First convert SMILES string to pybel mol
	py_mol = pb.readstring(l_in_fmt, l_smi)
	fpts_list.append(py_mol.calcfp(fptype=FPType))
    return fpts_list

def CalcTan(T_fpts):
    tanmat= np.array([l1 | l2 for l1 in T_fpts for l2 in T_fpts]) #this will be a N*N length array, where N = Total No of ligands in all input targets
    tanmat.shape= (len(T_fpts), len(T_fpts)) # Change the shape to N X N

    #norm_tc_val = (tc_val - np.min(tc_val))/ (np.max(tc_val)-np.min(tc_val))
    #Lets plot all similar ligands
    #print r_ind+1, c_ind+1
    #sl= plt.scatter(r_ind+1, c_ind+1, s=norm_tc_val*15, alpha=0.5) #+1 since we are plotting ligand number not indices 
    #plt.xticks(range(max(max(r_ind), max(c_ind))+2))
    #plt.yticks(range(max(max(r_ind), max(c_ind))+2))
    #plt.savefig('Similar_ligand.png')
	
    return tanmat

def Plot_TM(tc_mat, outfile, no_lig_tar=[], tick_labels=[]):
    f1=plt.figure(figsize=(8.0,8.0))
    #f1.subplots_adjust(left=0.18,bottom=0.15)
    plt.hold(True)
    ax = plt.gca()
    #cmap= colors.ListedColormap(['#FFFFFF', '#FF99FF','#FF00FF', '#B2FFFF','#006400'])
    cmap= colors.ListedColormap(['#FFFFFF', '#FF99FF','#FF00FF','#00FFFF','#FFD700','#FF0000','#0000CD','#006400'])
    #bounds=[0.0,0.2,0.4,0.6,0.8,1.0]
    bounds=[0.0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    img= plt.imshow(tc_mat,interpolation='nearest',origin='lower', vmin=0,vmax=1, cmap=cmap, norm=norm)
    if len(tick_labels) == 0:
	plt.xlabel('Ligand Number')
	plt.ylabel('Ligand Number')
    else:
	plt.xlabel('Ligand Id')
	plt.ylabel('Ligand Id')
	ax.set_xticks(range(len(tc_mat)))
	ax.set_yticks(range(len(tc_mat)))
	ax.set_xticklabels(tick_labels, rotation='vertical', fontsize=7)
	ax.set_yticklabels(tick_labels, fontsize=7)

    cbar=plt.colorbar(img, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds,shrink=.8, pad=.02)
    cbar.set_label('Tanimoto Similarity')
    #for t in cbar.ax.get_yticklabels():
	#t.set_fontsize(5)
    #plot lines defining different target boundaries 
    if len(no_lig_tar) != 0:
	start_ind = 0
	for l_no in no_lig_tar[:-1]: #Leave the last entry
		#t_indices = range(start_ind, start_ind+l_no)
		t_indices = range(start_ind+l_no)
		#Plot the horizontal line
		plt.plot(range(len(tc_mat)), [(start_ind+l_no-0.5)]*len(tc_mat), 'k--')  #0.5 because we want to position lines on the edge not center
		#Plot the vertical line
		#plt.plot([(start_ind+l_no-0.5)]*l_no, t_indices, 'k--')
		plt.plot([(start_ind+l_no-0.5)]*len(t_indices), t_indices, 'k--')
		plt.xlim(xmin=0)
		plt.ylim(ymin=0)  
		start_ind+= l_no
	
    plt.savefig(outfile, dpi=300)
    plt.close('all')    

def make_autopct(values):
    def my_autopct(pct):
	total = sum(values)
	val = int(round(pct*total/100.0))
	return '{p:.1f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def Find_Compounds(target_file_path,smiles_file_path,input_targets,fp_type,tc_cutoff, not_in_targets):
    #craf_file = open('Target_Ligands/595.smi', 'r')
    #craf_ids=[]
    #for line in craf_file:
	#line=line.split()
	#if line[1] not in craf_ids:
	#	craf_ids.append(int(line[1]))
    #Open target info file
    try:
	t_i_file= open(target_file_path, 'r')
    except IOError:
       sys.exit("Cannot locate file %s" % target_file_path + ".\n***Invalid filename or path***")

    #Make a list of input target identifiers
    t_identifier= input_targets.split(",")   
    if t_identifier == ['']: # No targets were passed
	sys.exit("***Exiting. No input target identifiers.***")
    elif len(set(t_identifier)) == 1:
	sys.exit("***Exiting. Not enough non-redundant targets to work with. Input at least 2 target identifiers.***")

    if len(set(t_identifier)) != len(t_identifier):
	logging.warning("\nFound redundant target identifiers in input target list.\n")
	logging.warning("\nWorking with non-redundant target identifiers set: %s\n", set(t_identifier))

    t_identifier = set(t_identifier)

    #Make a list of not_in_targets identifiers
    not_in_t= not_in_targets.split(",")
    if not_in_t != ['']: 
	if len(set(not_in_t)) != len(not_in_t):
		logging.warning("\nFound redundant target identifiers in NOT_IN_TARGETS list.\n")
        	logging.warning("\nWorking with non-redundant target identifiers set in NOT_IN_TARGETS: %s\n", set(not_in_t))

    not_in_t = set(not_in_t)	

    #Compile the list of ligands for all input targets
    t_l_smi= Target_Ligands(t_identifier, smiles_file_path) #key is the target id
    #for t_id in t_identifier:
        #Compile the list of ligand smiles for target t_id
	#for line in t_smi_file:
	#	line=line.split()
	#	l_smiles= line[0]
	#	l_id = int(line[1])
	#	if not t_l_smi[t_id].has_key(l_id):
	#		t_l_smi[t_id][l_id] = l_smiles
	#	else:
	#		sys.exit('Multiple entry for ligand id: %d for target id: %s' % (l_id, t_id)) #There must only be one entry per ligand id
	#t_smi_file.close()

    #Get the t_l_smi dict keys in ascending order of number of entries per target
    #sorted_tkeys_smi = sorted(t_l_smi, key= lambda k: len(t_l_smi[k]))

    #Get the t_l_smi dict keys in ascending order of target ids
    sorted_tkeys_smi = sorted(t_l_smi, key= lambda k: int(k))
    
    #Calculate the fingerprints for ligands of input targets
    logging.info("Using Open Babel to compute %s fingerprint of compounds for all input targets\n" % (fp_type))
    t_l_fpts={} #key is target id
    all_smiles_fpts=[] #List of ligand fpts in the order of target ids 
    t_l_no =[] #number of ligands per target (sorted by key)
    for t_id in sorted_tkeys_smi:
	t_l_fpts[t_id]= CalcFpt(t_l_smi[t_id][0], fp_type) #Pass smiles list
	all_smiles_fpts.extend(t_l_fpts[t_id])
	t_l_no.append(len(t_l_fpts[t_id]))

    #Calculate Tanimoto Coefficient between ligands of each target
    logging.info("Using Open Babel to compute Tanimoto coefficient between ligands of same target\n")
    #Go through each target
    for t_id in sorted_tkeys_smi:
	TC_matrix = CalcTan(t_l_fpts[t_id]) 
	#Plot the resulting matrix
	Plot_TM(np.triu(TC_matrix), 'TC_Matrix_'+t_id) #Take only upper half as matrix is symmteric
	logging.info("Number of ligands for target id %s: %d\n" % (t_id,len(t_l_fpts[t_id])))
	logging.info("Tanimoto matrix plot for target id %s: %s\n" % (t_id, 'TC_Matrix_'+t_id+'.png'))
	
    #Calculate tanimoto coefficient between ligands of all target
    logging.info("Using Open Babel to compute Tanimoto coefficient between ligands of all input targets\n")
    TC_matrix_all= CalcTan(all_smiles_fpts)
    #Plot the resulting matrix
    Plot_TM(np.triu(TC_matrix_all), 'TC_Matrix_all', t_l_no)
    logging.info("Tanimoto matrix plot for ligands of all input targets: %s\n" % ('TC_Matrix_all.png'))

    #Find similar ligands between two input targets. Use resulting ligands to find similar ligands with 3rd target and so on
    #TC_matrix_two_tar = CalcTan(t_l_fpts[sorted_tkeys_smi[0]]+t_l_fpts[sorted_tkeys_smi[1]])
    #Keep only upper half with diagonal set to 0
    #TC_two_tar_upper = np.triu(TC_matrix_two_tar, k=1)
    #Get the row and column index where tanimoto >= cutoff
    #r_ind_two_tar, c_ind_two_tar = np.where(TC_matrix_two_tar >= tc_cutoff)
    

    
    #Keep only upper half with diagonal set to 0
    TC_mat_all_upper = np.triu(TC_matrix_all, k=1) 
    #Get the row and column index where tanimoto >= cutoff
    r_ind, c_ind = np.where(TC_mat_all_upper >= tc_cutoff)
    #Get Tanimoto value for r_ind and c_ind
    tc_val_rc = np.array([TC_mat_all_upper[row,col] for row,col in zip(r_ind, c_ind)])    

    #Map r_ind, c_ind to target
    t_indices_rc = {}  #TC_mat_all indices corresponding to a particular target. Key: target id
    start_ind = 0
    for t_id, l_no in zip(sorted_tkeys_smi, t_l_no):
        t_indices_rc[t_id] = range(start_ind, start_ind+l_no)
        start_ind+= l_no

    #Create two lists, r_t_id (holds target id corresponding to ligand indices of r_ind) and c_t_id
    r_t_id=[]
    c_t_id=[]
    for row,col in zip(r_ind, c_ind):
	#Find target id of row
	for key in t_indices_rc:
		if row in t_indices_rc[key]:
			r_t_id.append(key)
			break

	#Find target id of col
        for key in t_indices_rc:
                if col in t_indices_rc[key]:
                        c_t_id.append(key)
                        break

    #Check if all entries were made
    try:
	assert len(r_t_id) == len(r_ind) and len(c_t_id) == len(c_ind)
    except AssertionError:
	sys.exit("***ERROR: Number of r_ind (or c_ind) does not match r_t_id (or c_t_id)***")

    #OUTPUT: Write similar ligands information
    logging.info("Writing information about similar ligands (pairwise) for all input targets to file: Similar_ligands_pair.tsv\n")
    file_sim_lig = open('Similar_ligands_pair.tsv', 'w')
    same_l_id={} #keeps track of same ligands occuring in multiple targets
    simi_l_id_t={} #Will be used to identify ligands that are similar and occur in multiple input targets. Key: ligand id
    simi_l_fpts=[] #Saves the fpts correspoding to similar ligands (will be used in case NOT list is provided)
    simi_l_targets=[] #Saves the targets info correspoding to similar ligands (to use in output and case NOT list is provided)
    simi_l_ids=[] #Saves the ligand ids correspoding to similar ligands (to use in output and case NOT list is provided)
    file_sim_lig.write('Ligand_ID\tLigand_ID\tTarget_ID\tTarget_ID\tTanimoto Similarity\n')

    for rind,cind,rtar,ctar,tc in zip(r_ind, c_ind, r_t_id, c_t_id, tc_val_rc):
	#Get the index of rind and cind in t_indices_rc. This will be the index for that particular ligand in t_l_smi
	index_rind = t_indices_rc[rtar].index(rind)
	index_cind = t_indices_rc[ctar].index(cind)
	#Now use this to figure out ligand id 
	l_id_r = t_l_smi[rtar][1][index_rind]
	l_id_c = t_l_smi[ctar][1][index_cind]

	#Skip the ligand entry if they belong to same target
	if rtar == ctar:
		#logging.info("Skipping similar ligands: %d and %d (ligands belong to same target ids: %s and %s)\n" % (l_id_r, l_id_c, rtar, ctar))
		continue
	elif l_id_r == l_id_c:
		#Make entries to same_l_id and simi_l_id_t 
		if not same_l_id.has_key(l_id_r):
			same_l_id[l_id_r]=[rtar, ctar]
		else:
			if rtar not in same_l_id[l_id_r]:
				same_l_id[l_id_r].append(rtar)
			if ctar not in same_l_id[l_id_r]:
                                same_l_id[l_id_r].append(ctar)

		if not simi_l_id_t.has_key(l_id_r):
                        simi_l_id_t[l_id_r]=[rtar, ctar]
                else:
                        if rtar not in simi_l_id_t[l_id_r]:
                                simi_l_id_t[l_id_r].append(rtar)
                        if ctar not in simi_l_id_t[l_id_r]:
                                simi_l_id_t[l_id_r].append(ctar)
		#Append fpts, target id, and ligand id for only row entry to avoid redundancy
		if l_id_r not in simi_l_ids:
			simi_l_ids.append(l_id_r)
			simi_l_targets.append(rtar)
			simi_l_fpts.append(t_l_fpts[rtar][index_rind])
	else:
		#Make entries to simi_l_id_t for both row and col
		if not simi_l_id_t.has_key(l_id_r):
                        simi_l_id_t[l_id_r]=[rtar, ctar]
                else:
                        if rtar not in simi_l_id_t[l_id_r]:
                                simi_l_id_t[l_id_r].append(rtar)
                        if ctar not in simi_l_id_t[l_id_r]:
                                simi_l_id_t[l_id_r].append(ctar)

		if not simi_l_id_t.has_key(l_id_c):
                        simi_l_id_t[l_id_c]=[rtar, ctar]
                else:
                        if rtar not in simi_l_id_t[l_id_c]:
                                simi_l_id_t[l_id_c].append(rtar)
                        if ctar not in simi_l_id_t[l_id_c]:
                                simi_l_id_t[l_id_c].append(ctar)

		#Append fpts, target id, and ligand id for both row and col entry
		if l_id_r not in simi_l_ids:
			simi_l_ids.append(l_id_r)
			simi_l_targets.append(rtar)
			simi_l_fpts.append(t_l_fpts[rtar][index_rind])

		if l_id_c not in simi_l_ids:
			simi_l_ids.append(l_id_c)
			simi_l_targets.append(ctar)
			simi_l_fpts.append(t_l_fpts[ctar][index_cind])

	file_sim_lig.write(str(l_id_r) + '\t' + str(l_id_c) + '\t' + rtar + '\t' + ctar + '\t' + str(tc) +'\n')
    file_sim_lig.close()
    logging.info("Finished writing information about similar ligands (pairwise) for all input targets to file: Similar_ligands_pair.tsv\n")

    #OUTPUT: Write information about same ligand occuring in multiple targets 
    logging.info("Number of ligands with same ligand id occuring in two or more input targets: %d\n" % (len(same_l_id),))
    logging.info("Information about these ligands written to file: Ligand_In_Multiple_Targets.tsv\n")
    file_same_lig = open('Ligand_In_Multiple_Targets.tsv', 'w')
    file_same_lig.write('Ligand_ID\tTarget_Ids\n')
    for key in same_l_id:
	file_same_lig.write(str(key))
	for t_id in same_l_id[key]:
		file_same_lig.write('\t'+t_id)
	file_same_lig.write('\n')
    file_same_lig.close()

    logging.info("Number of similar ligands that can potentially target two or more input targets: %d\n" % (len(simi_l_id_t),))

    #Filter similar ligands if similarity found with ligands of targets in not_in_targets list
    skip_simi_ligands=[] #Will hold ids of ligands to skip
    if not_in_t != set(['']):
	logging.info("Identifying if any similar ligands should be skipped due to similarity with ligands of targets in not_in_targets list.\n")
	#Compile the list of ligands for all targets in not_in_t
	t_l_smi_not= Target_Ligands(not_in_t, smiles_file_path) #key is the target id

	#Get the t_l_smi_not dict keys in ascending order of target ids
	sorted_tkeys_smi_not = sorted(t_l_smi_not, key= lambda k: int(k))

	#Calculate the fingerprints for ligands of targets in not_in_targets list
	logging.info("Using Open Babel to compute %s fingerprint of compounds for all targets in not_in_targets list\n" % (fp_type))
	t_l_fpts_not={} #key is target id
	t_l_no_not = [len(simi_l_ids)] #number of ligands per target (sorted by key) in not_in_targets list with first entry for no of similar ligands identified

	for t_id in sorted_tkeys_smi_not:
        	t_l_fpts_not[t_id]= CalcFpt(t_l_smi_not[t_id][0], fp_type) #Pass smiles list
        	simi_l_fpts.extend(t_l_fpts_not[t_id])
        	t_l_no_not.append(len(t_l_fpts_not[t_id]))

	#Calculate Tanimoto Coefficient between ligands of each target in not_in_targets list
	logging.info("Using Open Babel to compute Tanimoto coefficient between ligands of same target in not_in_targets list\n")
	#Go through each target
	for t_id in sorted_tkeys_smi_not:
		TC_matrix = CalcTan(t_l_fpts_not[t_id])
        	#Plot the resulting matrix
        	Plot_TM(np.triu(TC_matrix), 'TC_Matrix_NOT_'+t_id) #Take only upper half as matrix is symmteric
		logging.info("Number of ligands for target id %s in not_in_targets list: %d\n" % (t_id,len(t_l_fpts_not[t_id])))
		logging.info("Tanimoto matrix plot for target id %s in not_in_targets list: %s\n" % (t_id, 'TC_Matrix_NOT_'+t_id+'.png'))

	
	#Calculate tanimoto coefficient between all ligands in simi_l_fpts
	logging.info("Using Open Babel to compute Tanimoto coefficient between all identified smiliar ligands and ligands of targets in not_in_targets list\n")
	TC_matrix_simi_not= CalcTan(simi_l_fpts)
	#Plot the resulting matrix
	Plot_TM(np.triu(TC_matrix_simi_not), 'TC_matrix_simi_not', t_l_no_not)
	logging.info("Tanimoto matrix plot for all identified similar ligands and ligands of targets in not_in_targets list: TC_matrix_simi_not.png\n")

	#Keep only upper half with diagonal set to 0
	TC_simi_not_upper = np.triu(TC_matrix_simi_not, k=1)
	#Get the row and column index where tanimoto >= cutoff
	#r_ind_not, c_ind_not = np.where(TC_simi_not_upper >= tc_cutoff)
	r_ind_not, c_ind_not = np.where(TC_simi_not_upper >= 0.99)
    	#Get Tanimoto value for r_ind_not and c_ind_not
    	tc_val_rc_not = np.array([TC_simi_not_upper[row,col] for row,col in zip(r_ind_not, c_ind_not)])

	#Map c_ind_not to target_ids
	t_indices_c_not = {}  #TC_matrix_simi_not indices corresponding to a particular target. Key: target id
	start_ind = 0 
	for t_id_not, l_no_not in zip(['0']+sorted_tkeys_smi_not, t_l_no_not): #0 is padded as t_id for similar ligands (first entry of t_l_no_not)  
		t_indices_c_not[t_id_not] = range(start_ind, start_ind+l_no_not)
		start_ind+= l_no_not

	#Create a list, c_t_id_not that holds target id corresponding to ligand indices of c_ind_not
	c_t_id_not=[]
	for col in c_ind_not:
		#Find target id of col
		for key in t_indices_c_not:
			if col in t_indices_c_not[key]:
				c_t_id_not.append(key)
				break
	#Check if all entries were made
	try:
		assert len(c_t_id_not) == len(c_ind_not) 
	except AssertionError:
		sys.exit("***ERROR: Number of c_ind_not does not match c_t_id_not***")
	
	#Here we are interested in the similarity of first t_l_no_not[0] r_ind_not (r_ind_not in range 0 to t_l_no_not[0]-1) with c_ind_not (excluding first t_l_no_not[0] entries).
	#This tells us which similar ligands are also similar to ligands of targets in not_in_targets list
	rq_rind_list =  range(0, t_l_no_not[0])
	rq_cind_list =  range(t_l_no_not[0], len(TC_simi_not_upper))

	#Create an Output Folder for smiliar ligands that are to be skipped. 
	dir_skip = "Skipped_Ligands"
	if not os.path.isdir(dir_skip):
		os.mkdir(dir_skip)
	else:
		if os.listdir(dir_skip) != []:
			smifiles= glob.glob(os.path.join(dir_skip, '*.smi'))
			for filename in smifiles:
				os.remove(filename)

	logging.info("SMILE file for skipped ligand is in subdirectory Skipped_Ligands. Each SMILE file lists the SMILE of ligands in not_in_targets list that are similar to skipped ligand with the\
first SMILE entry for the skipped ligand. Format\nFilename: SkippedLigandId_TargetIDs.smi\nSMILE Identifier: SkippedLigandID_TargetID/LigandID_TargetID_TanimotoSimilarity.\n")
	logging.info("Following identified similar ligands are skipped due to similarity with ligands of targets in not_in_targets list\nLigand_ID_Sim\tLigand_ID_Not\tTC\n")

	files_skip_lig ={} #File handles for smiliar ligands to be skipped
	for rind,cind,ctar,tc in zip(r_ind_not, c_ind_not, c_t_id_not, tc_val_rc_not):
		if rind in rq_rind_list and cind in rq_cind_list:
			#Get the index of cind in t_indices_c_not. This will be the index for that particular ligand in t_l_smi_not
			index_cind = t_indices_c_not[ctar].index(cind)
			#Now use this to figure out ligand id and smiles
			l_id_c_not = t_l_smi_not[ctar][1][index_cind]
			l_smi_c_not = t_l_smi_not[ctar][0][index_cind]

			#rind is also the index of ligand id in simi_l_ids
			row_lid = simi_l_ids[rind]
			if row_lid not in skip_simi_ligands:
				skip_simi_ligands.append(row_lid)

			#Get the target id for rind
			row_tid = simi_l_targets[rind]
			#Get the index of row_lid in t_l_smi ligands list
                        index_lid = t_l_smi[row_tid][1].index(row_lid)
                        #Use this index to get the ligand SMILES
                        row_l_smi = t_l_smi[row_tid][0][index_lid]

			#Check if row ligand occurs in multipe receptor
			if row_lid in same_l_id:
				use_r_lid = str(row_lid)+'_m_'+'_'.join(same_l_id[row_lid])
			else:
				use_r_lid = str(row_lid)+'_'+row_tid
						
			logging.info("%s\t%s\t%0.3f\n" % (use_r_lid, str(l_id_c_not)+'_'+ctar, tc))
			#OUTPUT: Write to file 
			if not files_skip_lig.has_key(use_r_lid): #When opening the file make two SMILE entries: row_l_smi and l_smi_c_not 
				files_skip_lig[use_r_lid] = open(os.path.join(dir_skip, use_r_lid+'.smi'), 'w')
				files_skip_lig[use_r_lid].write(row_l_smi+ ' ' + use_r_lid+'\n') 
				files_skip_lig[use_r_lid].write(l_smi_c_not+ ' ' + str(l_id_c_not)+'_'+ctar+'_'+str('%0.3f' % tc)+'\n')
			else:
				files_skip_lig[use_r_lid].write(l_smi_c_not+ ' ' + str(l_id_c_not)+'_'+ctar+'_'+str('%0.3f' % tc)+'\n')	

	logging.info("Number of identified similar ligands to be skipped: %d\n" % (len(skip_simi_ligands),))

    #OUTPUT: Write similar ligands that are active against all input targets skipping any if required
    logging.info("Ligands smiles file for ligands that could potentially act on all input targets: Similar_Ligands_All.smi\n")
    logging.info("Information about similar ligands and their potential targets: Similar_Ligands_Targets.tsv\n")
    file_similar_smiles = open('Similar_Ligands_All.smi', 'w')
    file_similar_targets = open('Similar_Ligands_Targets.tsv', 'w')
    file_similar_targets.write('Ligand_ID\tPotential_Target_IDs\n')
    count_lig = 0 #Number of ligands that could be active against all targets
    t_count={} #number of similar ligands chosen per target
    l_id_bar = [] #for bar plot
    l_id_bar_t_count = [] #for bar plot
    active_l_all_t_fpts =[] #for final TM plot of selected ligands that are potentially active against all input targets
    active_l_id=[]
    pie_t_count ={} #For pie plot; key: Number of targets. 

    for l_id, t_id, l_fpts in zip(simi_l_ids, simi_l_targets, simi_l_fpts):
	if l_id not in skip_simi_ligands:
		#Append all target_ids to ligand_id, check if ligand occurs in multipe receptor
		if l_id in same_l_id:
			use_lid = str(l_id)+'_m_'+'_'.join(same_l_id[l_id])
		else:
			use_lid = str(l_id)+'_'+t_id

		#if l_id in craf_ids:
		#	use_lid = use_lid + '_craf'
		#Write SMILES for ligands that could be active against all input targets
		if len(simi_l_id_t[l_id]) == len(t_identifier):
			#Get the index of l_id in t_l_smi ligands list
			index_lid = t_l_smi[t_id][1].index(l_id)
			#Use this index to get the ligand SMILES
			l_smi = t_l_smi[t_id][0][index_lid]
			file_similar_smiles.write(l_smi+ ' ' + use_lid + '\n')
			count_lig+=1
			#Increment t_count
			if not t_count.has_key(t_id):
				t_count[t_id]=1
			else:
				t_count[t_id]+=1
			#Append ligand fingerprint
			active_l_all_t_fpts.append(l_fpts)
			active_l_id.append(str(l_id)+'_'+t_id)
		#Write target info for all similar ligands that do not need to be skipped
		file_similar_targets.write(use_lid)
		l_id_bar.append(str(l_id)+'_'+t_id)
		temp_count_t = 0 
		for id_target in simi_l_id_t[l_id]:
			file_similar_targets.write('\t'+id_target)
			temp_count_t+=1
		file_similar_targets.write('\n')
		l_id_bar_t_count.append(temp_count_t)
		if not pie_t_count.has_key(temp_count_t):
			pie_t_count[temp_count_t] = 1
		else:
			pie_t_count[temp_count_t]+=1
    file_similar_smiles.close()
    file_similar_targets.close()
    logging.info("%d similar ligands SMILES written to file Similar_Ligands_All.smi\n" % (count_lig,))

    if len(active_l_all_t_fpts) !=0:
	logging.info("Tanimoto plot for %d similar ligands: TM_Similar_Ligands_All.png\n" % (count_lig,))
	TC_mat_sim_all = CalcTan(active_l_all_t_fpts)
	#Plot the resulting matrix
	Plot_TM(np.triu(TC_mat_sim_all), 'TM_Similar_Ligands_All.png', [], active_l_id) #Take only upper half as matrix is symmteric
    
    t_sim_l_count=[] #for Bar plot
    t_id_sim_l_count=[] #for Bar plot
    no_l_t_sim_l=[] #for bar plot
    for key in t_count:
	logging.info("Number of similar ligands from target id %s: %d\n" % (key,t_count[key]))
	t_sim_l_count.append(t_count[key])
	t_id_sim_l_count.append(key)
	no_l_t_sim_l.append(len(t_l_smi[key][1]))

    #Make a BAR plot for number of ligands per target
    if len(t_sim_l_count) != 0:
	x_index = np.arange(len(t_id_sim_l_count)) #X-index for the t_ids
	b_width = 0.35 #Width of bars
    	fig, ax= plt.subplots()
    	bar1 = ax.bar(x_index, no_l_t_sim_l, b_width, color='b')
    	bar2 = ax.bar(x_index+b_width, t_sim_l_count, b_width, color='y')
    	ax.set_ylabel('Number of Ligands')
    	ax.set_xlabel('Target id')
    	ax.set_xticks(x_index+b_width)
    	ax.set_xticklabels(t_id_sim_l_count)
    	ax.legend((bar1[0], bar2[0]), ('Total', 'Selected Similar'))
    	plt.savefig('Bar_No_of_Ligands.png')
    	plt.close('all')

    #Make a BAR plot for number of targets per ligand
    if len(l_id_bar_t_count) != 0:
	x_index = np.arange(len(l_id_bar)) 
	fig, ax= plt.subplots(figsize=(24.0,8.0))
	fig.subplots_adjust(bottom=0.25)
	bar1 = ax.bar(x_index, l_id_bar_t_count, color='y')
	ax.set_ylabel('Number of Targets')
	ax.set_xlabel('Ligand id')
	ax.set_xticks(x_index)
	ax.set_yticks(range(len(t_identifier)+2))
	ax.set_xticklabels(l_id_bar, rotation='vertical', fontsize=2)
	plt.savefig('Bar_No_of_Targets.png', dpi=300)
	plt.close('all')

	#Make a pie plot of % number of targets per ligand 
	t_no_labels, t_values = pie_t_count.keys(), pie_t_count.values()
	fig= plt.subplots(figsize=(4.0,4.0))
	patches, texts, autotexts = plt.pie(t_values, labels=t_no_labels, autopct=make_autopct(t_values), pctdistance=0.75)
	for AT in autotexts:
		AT.set_color('white')
		AT.set_fontsize(8)
	plt.savefig('Pie_No_of_Targets.png', dpi=300)
	plt.close('all')
   	
