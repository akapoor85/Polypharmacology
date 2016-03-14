#!/usr/bin/python

import sys
import logging
import numpy as np


def Filter_Sequence(reference,target_ids):
    try:
	ref_file= open(reference, 'r')
    except IOError:
       sys.exit("Cannot locate file %s" % reference + ".\n***Invalid filename or path***")

    try:
	query_target_ids= np.genfromtxt(target_ids, usecols=(0,),delimiter=',',dtype='string', skiprows=1)
    except IOError:
       sys.exit("Cannot locate file %s" % target_ids + ".\n***Invalid filename or path***")

    count_seq=0  #To keep track of number of sequences in reference file
    ref_seq_stat ={} #Key: Target_id, value: [[sequence1,sequence2,..],total_number_of_sequence_with_Target_id]
    no_warning=0 #Keep track of warning messages

    #Populate ref_seq_stat dict
    new_seq=False
    continue_seq=False
    for line in ref_file:
	if line[0]=='>':
		new_seq=True
		continue_seq=False
		line_split=line.split()
		ref_target_id=line_split[0][line_split[0].find('|')+1:].upper()   #This is the drugbank target id given in FASTA header
		headerline=line
		count_seq+=1
	elif new_seq:
		new_seq=False
		continue_seq=True
		if not ref_seq_stat.has_key(ref_target_id):
			ref_seq_stat[ref_target_id]=[[headerline+line],1] #line is the first line of sequence entry, 1 is total number of sequence
		else:                                         #If the target id is already present then a new sequence is started and the number of sequence is incremented
			ref_seq_stat[ref_target_id][0].append(headerline+line)
			#ref_seq_stat[ref_target_id][0].append(line)
			ref_seq_stat[ref_target_id][1]+=1
	elif continue_seq:
		ref_seq_stat[ref_target_id][0][ref_seq_stat[ref_target_id][1]-1]+=line   #append the sequence chunk in current line; (ref_seq_stat[DB_target_id][1]-1) will append to latest copy 

    #Reference Target Sequence: Update Log File
    logging.info("\tNumber of sequences in reference file: %d", count_seq)
    logging.info("\tNumber of UNIQUE target ids in reference file: %d", len(ref_seq_stat))
    for key in ref_seq_stat:
	if ref_seq_stat[key][1]>1:
		logging.warning("\tMultiple sequence entry for target_id: %s; Number of sequences: %d", key, ref_seq_stat[key][1])
		no_warning+=1

    #Target_ids from target_file: Update Log file
    logging.info("\tNumber of target ids in target identifier file: %d", len(query_target_ids))
    
    if len(query_target_ids)-len(set(query_target_ids)) > 0:
	logging.warning("\tNumber of duplicate target_ids in target identifier file: %d", len(query_target_ids)-len(set(query_target_ids)))
	unique=[]
	duplicate=[]
	for t_id in query_target_ids:
		if t_id.strip() not in unique:
			unique.append(t_id.strip())
		else:
			duplicate.append(t_id.strip())
	logging.warning("\tList of duplicate target_ids:\n"+str(duplicate) )
	no_warning+=1

    #Write sequences of targets in target_id from reference sequence file
    outfile = open('Filtered_Sequences.fasta', 'w')
    count_seq_written=0
    tid_no_seq=[]
    for t_id in unique:
	if ref_seq_stat.has_key(t_id):
		for line in ref_seq_stat[t_id][0]:
			outfile.write(line)
			if line[0]=='>':
				count_seq_written+=1
	else:
		tid_no_seq.append(t_id)
    logging.info("%d reference sequences written to file Filtered_Sequences.fasta", count_seq_written)
    logging.info("\tNumber of target ids in target identifier file with no sequence in reference sequence file: %d", len(tid_no_seq))
    logging.info("\nRun Completed with %d Warnings", no_warning)
    print "Run Completed with %d Warnings: Check logfile log_runinfo" % (no_warning,)


