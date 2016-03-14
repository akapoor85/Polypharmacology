#!/usr/bin/python

import sys
import os, glob
import logging
import numpy as np
from rdkit import Chem
import sanifix4
from shutil import copy as cp

def Compile_Target_Info(t_i_fileObject):
    t_info_dict = {}
    #Skip firstline, as it is header
    t_i_fileObject.next()
    for line in t_i_fileObject:
        line=line.strip('\n').split('\t')
        if not t_info_dict.has_key(line[0]):
                t_info_dict[line[0]] = [int(line[1]),line[2:]]
        else:
                sys.exit('Multiple entry for target id: %s in file: %s' % (line[0], target_file_path))
    return t_info_dict

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

def Sanitize_SMI(t_l_smi_id_lists):
    problem_ligand = {} #Key: Ligand id. Value: problem_msz. 
    replace_SMILE = {} #Key: Ligand id. Value: Fixed_SMILE. For problem ligands that were successfully sanitized (Only kekulize cases)
    for l_smi, l_id in zip(t_l_smi_id_lists[0], t_l_smi_id_lists[1]):
	#Convert SMILES string to RDKit mol
	rd_mol = Chem.MolFromSmiles(l_smi, sanitize=False)
	#Sanitize mol
	try:
		Chem.SanitizeMol(rd_mol)
	except ValueError as err:
		if not problem_ligand.has_key(l_id):
			problem_ligand[l_id] = err.message
		if "can't kekulize" in err.message.lower():
			#Attempt to sanitize it with sanifix4
			rd_mol_new = Chem.MolFromSmiles(l_smi, sanitize=False)
			rd_mol_fixed = sanifix4.AdjustAromaticNs(rd_mol_new)
			#If sanitization was successful then send fixed SMILE to replace original SMILE
			if rd_mol_fixed is not None:
				if not replace_SMILE.has_key(l_id):
					replace_SMILE[l_id] = Chem.MolToSmiles(rd_mol_fixed)
    return problem_ligand, replace_SMILE

def main(argv):
    #######Input#######

    t_info_file_path = 'Targets_for_ligands.tsv'
    ligand_smi_folder = 'Target_Ligands'  #Path to folder containing ligand SMILES files
    MIN_LIGANDS = 10 #If any target has less than MIN_LIGANDS ligands after sanitization then warning msz will appear at the end of log file

    #######Configuration for logging#######
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s', filename='log_sanitize', filemode='w')

    #Open target info file
    try:
        t_i_file= open(t_info_file_path, 'r')
    except IOError:
       sys.exit("Cannot locate file %s" % t_info_file_path + ".\n***Invalid filename or path***")

    #Read data from target_info_file to a dict. Key: target_id, value: [No_of_ligands, [target names]]
    t_info_dict = Compile_Target_Info(t_i_file)
    t_i_file.close()
    
    #Compile the list of ligands for all targets
    t_l_smi, no_warn= Target_Ligands(t_info_dict.keys(), ligand_smi_folder) #key: target id, value: [[ligand_SMILEs], [ligand_ids]]
    if no_warn:
	print "One or more ligand SMILE files were not read. Check log_sanitize for more details."

    #Get a list of sorted target ids
    sorted_tkeys_smi = sorted(t_l_smi, key= lambda k: int(k))

    #Open a new file to write updated target information
    file_targets_info = open('Targets_for_ligands_sanitized.tsv', 'w') #For each .smi file write their corresponding target info
    file_targets_info.write('Target_ID'+ '\t' + 'No_of_Ligands' + '\t' + 'Target_Name(s)')

    #Create an Output Folder for sanitized smi files 
    if not os.path.isdir("Target_Ligands_Sanitized"):
        os.mkdir("Target_Ligands_Sanitized")
    else:
        if os.listdir("Target_Ligands_Sanitized") != []:
                smifiles= glob.glob(os.path.join("Target_Ligands_Sanitized", "*.smi"))
                for filename in smifiles:
                        os.remove(filename)

    #Go through each target_id smi file and sanitize (if required)
    uni_rem_ligs = set() #Keeps unique ligand id that are removed
    uni_fix_ligs =set() #Keeps unique ligand id that are fixed
    no_ligands_rem = 0
    target_warn = []
    logging.info("Sanitized Files Information: ")
    for t_id in sorted_tkeys_smi:
	problem_lig, replace_SMI = Sanitize_SMI(t_l_smi[t_id]) #pass [[ligand_SMILES], [ligand_ids]] for t_id
	#If not empty then ligands with problem were identified.
	if problem_lig: 
		logging.info("\n\t%s.smi: Number of problem ligands: %s; Number Fixed: %s" % (t_id, len(problem_lig), len(replace_SMI)))
		#Write a new smi file containing sanitized ligands only
		file_smi = open(os.path.join("Target_Ligands_Sanitized",str(t_id)+".smi"), 'w')

		for lig_smi, lig_id in zip(t_l_smi[t_id][0], t_l_smi[t_id][1]):
			if not problem_lig.has_key(lig_id):
				file_smi.write(lig_smi + ' ' + str(lig_id) + '\n')
			elif replace_SMI.has_key(lig_id):
				file_smi.write(replace_SMI[lig_id] + ' ' + str(lig_id) + '\n')
				logging.info("\t\tFIXED %s %s %s\n\t\t\tReplaced SMILE: %s" % (lig_id, lig_smi, problem_lig[lig_id], replace_SMI[lig_id]))
			else:
				logging.info("\t\tREMOVED %s %s %s:" % (lig_id, lig_smi, problem_lig[lig_id]))
		file_smi.close()
		#Determine unique ligand ids removed
		ligands_id_removed = set.difference( set(problem_lig.keys()), set(replace_SMI.keys()))
		no_ligands_rem+=len(ligands_id_removed)
		uni_rem_ligs = set.union(uni_rem_ligs, ligands_id_removed)
		uni_fix_ligs = set.union(uni_fix_ligs, set(replace_SMI.keys()))
		#Write to Target info file
		file_targets_info.write('\n' + str(t_id) + '\t' + str(t_info_dict[t_id][0] - len(ligands_id_removed)))
		for name in t_info_dict[t_id][1]:
			file_targets_info.write('\t' + name)
		if (t_info_dict[t_id][0] - len(ligands_id_removed)) < MIN_LIGANDS:
			target_warn.append(t_id)
	else:
		src = os.path.join(ligand_smi_folder,str(t_id)+".smi")
		dest = os.path.join("Target_Ligands_Sanitized",str(t_id)+".smi")
		cp(src,dest)
		#Write to Target info file
		file_targets_info.write('\n' + str(t_id) + '\t' + str(t_info_dict[t_id][0]))
		for name in t_info_dict[t_id][1]:
			file_targets_info.write('\t' + name)
    logging.info("Sanitization Complete!!!\nNumber of ligands removed (includes repetition): %s.\n\t\t Number of unique problematic ligands removed: %s" % (no_ligands_rem,len(uni_rem_ligs,)))
    logging.info("Number of unique problematic ligands fixed: %s" % (len(uni_fix_ligs)))
    logging.info("New SMILE files are in directory: Target_Ligands_Sanitized.\nNew Target Info file: Targets_for_ligands_sanitized.tsv")
    if target_warn:
	logging.warning("Following targets have less than %s ligands after sanitization" % (MIN_LIGANDS,))
	for t_id in target_warn:
		logging.warning("Target_id: %s" % (t_id,))

if __name__ == '__main__':
    sys.exit(main(sys.argv))

