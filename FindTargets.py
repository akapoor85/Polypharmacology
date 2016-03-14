#!/usr/bin/python

from __future__ import division
import sys
import os, glob
import logging
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import pybel as pb
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols #For linear fingerprint (similar to Daylight)
from rdkit.Chem import AllChem #For cicular (ECFP4) fingerprint
import sanifix4
from operator import itemgetter


def Target_Ligands(t_iden_list, smiles_file_path):
    no_warn=0
    t_l_s={} #key is the target id
    for t_id in t_iden_list:
	try:
                t_smi_file = open(os.path.join(smiles_file_path, t_id+'.smi'), 'r')
        except:
		no_warn+=1
                logging.warning("Cannot locate file %s" % os.path.join(smiles_file_path, t_id+'.smi') + ".\n***Skipping target entry***")
		continue

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
    return t_l_s, no_warn

def CalcFpt(t_smi_list, FPType, usefp=''):
    fpts_list=[]
    l_in_fmt = 'smi' #Input format is SMILES
    if usefp != '': #Use OpenBabel
	for l_smi in t_smi_list:
		#First convert SMILES string to pybel mol
		py_mol = pb.readstring(l_in_fmt, l_smi)
		fpts_list.append(py_mol.calcfp(fptype=usefp))
    else: #Use RDKit
	count = 0 #To kep track of ligand index; Required if RDKit fails to generate molecule 
	for l_smi in t_smi_list:
		#First convert SMILES string to RDKit mol
		rd_mol = Chem.MolFromSmiles(l_smi, sanitize=False)
		#Sanitize mol
		try:
			Chem.SanitizeMol(rd_mol)
			count+=1
                        if FPType == 'linear':
                                fpts_list.append(FingerprintMols.FingerprintMol(rd_mol))
                        else:
                                fpts_list.append(AllChem.GetMorganFingerprintAsBitVect(rd_mol,2,nBits=16384))
		except ValueError as err:
			if "can't kekulize" in err.message.lower():
				#Attempt to sanitize it with sanifix4
				rd_mol_new = Chem.MolFromSmiles(l_smi, sanitize=False)
				rd_mol_fixed = sanifix4.AdjustAromaticNs(rd_mol_new)
				if rd_mol_fixed is not None: #sanitization was successful
					count+=1
					if FPType == 'linear':
						fpts_list.append(FingerprintMols.FingerprintMol(rd_mol_fixed))
					else:
						fpts_list.append(AllChem.GetMorganFingerprintAsBitVect(rd_mol_fixed,2,nBits=16384))
				else:
					return ('exit', count, err.message + '**LIGAND COULD NOT BE KEKULIZED WITH SANIFIX4 AS WELL**')
			else:
				#Exit if sanitization failed
				return ('exit', count, err.message)
    return fpts_list

def CalcTan(q_fpt, t_l_fpt_list, p_type):
    #Calculate Tanimoto similarity between q_ligand and t_ligands
    if p_type.lower() == 'openbabel':
	tc_array = np.array([q_fpt | l1 for l1 in t_l_fpt_list]) #this will be a length N array, where N = Total No of ligands in input target
    elif p_type.lower() == 'rdkit':
	tc_array = np.array([DataStructs.TanimotoSimilarity(q_fpt, l1) for l1 in t_l_fpt_list])
    else:
	sys.exit("***ERROR: Unrecognized Package Type when calculating Tanimoto similarity: %s***\nValid options: openbabel or rdkit" % p_type)
    return tc_array

def Find_Targets(target_file_path,smiles_file_path,query_smi_file,fp_type,package_type,tc_cutoff, K_KNN):
    no_warnings=0

    #Select fp2 fingerprint for OpenBabel 
    if package_type == 'openbabel':
	#Set fingerprint to fp2
	use_fp = 'fp2'
	
    #Open target info file
    try:
	t_i_file= open(target_file_path, 'r')
    except IOError:
       sys.exit("Cannot locate file %s" % target_file_path + ".\n***Invalid filename or path***")

    #Open query smi file
    try:
        q_smi_file= open(query_smi_file, 'r')
    except IOError:
       sys.exit("Cannot locate file %s" % query_smi_file + ".\n***Invalid filename or path***")

    #Read data from target_file_path to a dict. Key: target_id, value: [No_of_ligands, [target names]]
    t_info_dict = {}
    #Skip firstline, as it is header
    t_i_file.next()
    for line in t_i_file:
	line=line.strip('\n').split('\t')
	if not t_info_dict.has_key(line[0]):
		t_info_dict[line[0]] = [int(line[1]),line[2:]] 
	else:
		sys.exit('Multiple entry for target id: %s in file: %s' % (line[0], target_file_path)) 
    t_i_file.close()

    #Compile the list of ligands for all targets
    t_l_smi, no_warn= Target_Ligands(t_info_dict.keys(), smiles_file_path) #key: target id, value: [[ligand_SMILEs], [ligand_ids]]
    no_warnings+=no_warn

    #Get a list of sorted target ids
    sorted_tkeys_smi = sorted(t_l_smi, key= lambda k: int(k))

    #Calculate the fingerprints for ligands of all targets
    s_time = time.clock()
    t_l_fpts={} #key: target id, value: [ligand fingerprints]
    total_lig = 0
    if package_type == 'openbabel':
	logging.info("Using OpenBabel to compute linear (%s) fingerprint of compounds for all targets\n" % (use_fp))
	for t_id in sorted_tkeys_smi:
		t_l_fpts[t_id]= CalcFpt(t_l_smi[t_id][0], 'linear' ,use_fp) #Pass smiles list
		if len(t_l_fpts[t_id]) != len(t_l_smi[t_id][0]):#Although this will not happen as OpenBabel is not very picky about SMILES making chemical sense
			sys.exit("OpenBabel failed to generate fingerprint for one or more ligands.\nRemove problematic ligands and try again!")
		total_lig += len(t_l_fpts[t_id])
    elif package_type == 'rdkit':
	if fp_type == 'linear':
		logging.info("Using RDKit to compute linear (similar to daylight) fingerprint of compounds for all targets\n")
	else:
		logging.info("Using RDKit to compute circular (ECFP4 with bit length 16384) fingerprint of compounds for all targets\n")
        for t_id in sorted_tkeys_smi:
		if fp_type == 'linear':
			t_l_fpts[t_id]= CalcFpt(t_l_smi[t_id][0], 'linear') #Pass smiles list
		else:
			t_l_fpts[t_id]= CalcFpt(t_l_smi[t_id][0], 'circular') #Pass smiles list
		#Exit if problem generating fingerprints
		if t_l_fpts[t_id][0] == 'exit':
			f_index = t_l_fpts[t_id][1] #Index of ligand id for which RDKit failed to generate molecule
			print "RDKit failed to generate molecule from SMILES for ligand id %s in target id %s" % (t_l_smi[t_id][1][f_index], t_id)
			print "RDKit Error Message:\n%s" % t_l_fpts[t_id][2]
			sys.exit("\n***Fix this molecule and try again!***")
                total_lig += len(t_l_fpts[t_id])
    else:
	sys.exit("***ERROR: Unrecognized Package Type: %s***\nValid options: openbabel or rdkit" % package_type)
    e_time = time.clock()
    logging.info("Finished computing fingerprints.\n")
    logging.info("Fingerprints computed for %d ligands.\n" % (total_lig,))
    logging.info("Time spent computing fingerprints: %0.3f seconds process time\n" % (e_time-s_time,))

    #Read data from query_smi_file to a dict. Key: ligand_id, value: ligand SMILE
    q_smi_dict = {}
    for line in q_smi_file:
        line=line.strip('\n').split()
	if line: #if not an empty line
		if not q_smi_dict.has_key(line[1]):
			q_smi_dict[line[1]] = line[0] 
		else:
			no_warnings+=1
			logging.warning("Multiple entry for query ligand id: %s. Multiple Entry Skipped.\n" % line[1])
    q_smi_file.close()

    #Compile similarity data for each query SMILE 
    for q_id in q_smi_dict:
	#Calculate fingerprint of query ligand
        if package_type == 'openbabel':
                q_l_fpt = CalcFpt([q_smi_dict[q_id]], 'linear' ,use_fp) #q_l_fpt is a list of length 1
        elif package_type == 'rdkit':
                if fp_type == 'linear':
                        q_l_fpt = CalcFpt([q_smi_dict[q_id]], 'linear') #q_l_fpt is a list of length 1
                else:
                        q_l_fpt = CalcFpt([q_smi_dict[q_id]], 'circular') #q_l_fpt is a list of length 1
                #Skip if problem generating fingerprints
                if q_l_fpt[0] == 'exit':
                        logging.warning("Query Ligand Skipped. RDKit failed to generate molecule from SMILE for query ligand %s.\n\t\t**RDKit Error Message**:\n\t\t%s" % (q_id, q_l_fpt[2]))
                        continue

	#Create an output directory for q_id
	if not os.path.isdir(q_id):
		os.mkdir(q_id)
	else:
		if os.listdir(q_id) != []:
			smifiles= glob.glob(os.path.join(q_id, "*.smi"))
			for filename in smifiles:
				os.remove(filename)
	#Open a file to write information about identified targets
	iden_t_file = open(os.path.join(q_id, 'Targets_Identified.tsv'), 'w')
	iden_t_file.write('Target_ID\tNo_of_Ligands\tNo_Similar\tPercent_Similar\tMax_Score (1NN)\t'+ str(K_KNN)+'NN\tCentroid_Score\tTarget_Name(s)\n')
	#Create list to hold data for identfied targets
	iden_t_id = []
	per_l_iden_t = [] #Percent ligand in target similar to query ligand 
	iden_t_data = [] #Holds all information related to identifed target. Format: [(tid,No_of_lig,No_similar,Percent_Similar,Max_Score,K_KNN,Centroid_Score)] 
	#Go through each target
	for t_id in sorted_tkeys_smi:
		#Get the Tanimoto Similarity between q_ligand and ligands of target t_id
		if package_type == 'openbabel':
			tc_vals = CalcTan(q_l_fpt[0], t_l_fpts[t_id], 'openbabel')
		elif package_type == 'rdkit':
			tc_vals = CalcTan(q_l_fpt[0], t_l_fpts[t_id], 'rdkit')
		
		#Get the index of tc_vals where tanimoto >= cutoff; this will also be the index of ligand in t_l_smi[t_id] 
		index_tc = list(np.where(tc_vals >= tc_cutoff)[0])
		if index_tc: #if index_tc is not-empty
			#Update identified target list
			iden_t_id.append(t_id)
			Total_l_no = t_info_dict[t_id][0]
			per_l_iden_t.append(len(index_tc)/Total_l_no * 100)

			#Get the KNN scores
			CS = np.mean(tc_vals)
			NN_1 = tc_vals.max()
			NN_K = np.mean(sorted(tc_vals, reverse=True)[0:K_KNN])
			
			#Update iden_t_data
			iden_t_data.append((t_id,Total_l_no,len(index_tc),per_l_iden_t[-1],NN_1,NN_K,CS))
			#Write information about identified targets
			#iden_t_file.write(t_id + '\t' + str(Total_l_no) + '\t' + str(len(index_tc)) + '\t'+ str(per_l_iden_t[-1]) + '\t' + str(NN_1) + '\t' + str(NN_K) + '\t' + str(CS))
			#for t_name in t_info_dict[t_id][1]:
			#	iden_t_file.write('\t' + t_name)
			#iden_t_file.write('\n')

			#Open a file to write all similar ligands 
			sim_l_file = open(os.path.join(q_id, t_id+'.smi'), 'w')
			
			for ind in index_tc:
				tar_lig_id = t_l_smi[t_id][1][ind]
				tar_lig_smi = t_l_smi[t_id][0][ind] 
				tan_coff = str('%0.3f' % tc_vals[ind])
				sim_l_file.write(tar_lig_smi+ ' ' + str(tar_lig_id) +'_'+ t_id + '_' + tan_coff + '\n')
			sim_l_file.close()
	#Write the identified target data sorted by KNN score in descending order
	for entry in sorted(iden_t_data, key=itemgetter(5), reverse=True):
		#Write information about identified targets
		iden_t_file.write(entry[0] + '\t' + str(entry[1]) + '\t' + str(entry[2]) + '\t'+ str(entry[3]) + '\t' + str(entry[4]) + '\t' + str(entry[5]) + '\t' + str(entry[6]))
		for t_name in t_info_dict[entry[0]][1]:
			iden_t_file.write('\t' + t_name)
		iden_t_file.write('\n') 
	iden_t_file.close()
	#BAR plot of percent similar ligand for identified targets
	if iden_t_id:
		x_index = np.arange(len(iden_t_id))
		fig, ax= plt.subplots(figsize=(21.0,7.0))
		fig.subplots_adjust(bottom=0.25)
		bar1 = ax.bar(x_index, per_l_iden_t, color='y')
		ax.set_ylabel('Percent Similar Ligands')
		ax.set_xlabel('Target_Id')
		ax.set_xticks(x_index)
		ax.set_yticks(range(0,110,10))
		ax.set_xticklabels(iden_t_id, rotation='vertical', fontsize=7)
    		plt.savefig(os.path.join(q_id,'Percent_Similar_Ligands.png'), dpi=300)
    		plt.close('all')
 
    if no_warnings:
	print "Run completed with %d Warnings: Check log_runinfo for more details." % (no_warnings,)
		
