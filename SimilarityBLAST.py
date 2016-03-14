#!/usr/bin/python

import sys
import os
import logging
import numpy as np
import subprocess


def Run_BLAST(blast_type,query_sequence,database_path,eval_cutoff,output_filename):
    '''
    BLAST OUTPUT FORMAT: Tabular
    COLUMN 1: Query Seq id
    COLUMN 2: Subject Seq id
    COLUMN 3: Percent Identity
    COLUMN 4: Alignment Length
    COLUMN 5: Query Seq length
    COLUMN 6: Subject Seq length
    COLUMN 7: Start of alignment in query
    COLUMN 8: End of alignment in query
    COLUMN 9: Start of alignment in subject
    COLUMN 10: End of alignment in subject
    COLUMN 11: Evalue
    COLUMN 12: Bitscore
    '''

    blast_output_format = '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore'
    proc=subprocess.Popen([blast_type, '-query', query_sequence, '-out', output_filename, '-db', database_path, '-evalue', eval_cutoff, '-outfmt', blast_output_format],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout_data, stderr_data)= proc.communicate()
    logging.info("BLAST Similarity Search output: %s\n", stdout_data)
    if stderr_data != '':
        logging.error("***BLAST Similarity Search ERROR***\n%s", stderr_data)
        sys.exit("***ERROR running BLAST similarity search***: Check log_runinfo")
    else:
	print "Success: Similarity Search Complete!!"
        logging.info("BLAST Similarity Search Output is in file: %s\n", output_filename)

def Extract_Hits(output_filename, database_path, dbtype):
    file_hits_temp = open('hits_temp', 'w')
    hitids = np.genfromtxt(output_filename, usecols=(1,), dtype='string')
    for name in set(hitids):    #set(hitnames) will remove duplicate entries
	file_hits_temp.write(name+'\n')
    file_hits_temp.close()

    #Extract hit sequences to file hits_output_filename
    proc1=subprocess.Popen(['blastdbcmd', '-db', database_path, '-dbtype', dbtype, '-entry_batch', 'hits_temp' ,'-outfmt','%f', '-out', 'hits_'+output_filename],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    (stdout_data_p1, stderr_data_p1)= proc1.communicate()
    if stderr_data_p1 != '':
        logging.error("***BLAST hit Sequence retrieval ERROR***\n%s", stderr_data_p1)
        sys.exit("***ERROR running BLAST command blastdbcmd***: Check log_runinfo")
    else:
        logging.info("Hit sequences written to file: %s\n", 'hits_'+output_filename)

def Extract_Hit_Names(output_filename):
    h_names={} 
    file_hits_name=open('hits_'+output_filename, 'r')
    for line in file_hits_name:
	if line[0]=='>':
		line=line.split()
		id_for_name=line[0][line[0].find('|')+1:]
		if not h_names.has_key(id_for_name):
			h_names[id_for_name]=[' '.join(line[3:])]
		else:
			h_names[id_for_name].append(' '.join(line[3:]))
    return h_names


def Similarity_Search(blast_type,query_sequence,database_path,eval_cutoff,output_filename,dbtype):
    #Call to Run_BLAST
    Run_BLAST(blast_type,query_sequence,database_path,eval_cutoff,output_filename)
    
    #Extract hits (entries of local BLAST database) name from BLAST output
    Extract_Hits(output_filename, database_path, dbtype)

    #Compile hit-query mapping
    hit_map={}
    output_data = np.genfromtxt(output_filename, dtype=None) #dtype=None will guess the datatype of each column
    for line in output_data:
	query_id = line[0]
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
	hit_coverage = '%.2f' % (100*(float((e_a_hit-s_a_hit+1))/hit_length))
	query_coverage = '%.2f' % (100*(float((e_a_query-s_a_query+1))/query_length))
	
	if not hit_map.has_key(hit_id):
		hit_map[hit_id]= [[query_id,p_identity,a_length,query_length,hit_length,s_a_hit,e_a_hit,hit_coverage,query_coverage,evalue]]
	else:
		hit_map[hit_id].append([query_id,p_identity,a_length,query_length,hit_length,s_a_hit,e_a_hit,hit_coverage,query_coverage,evalue])

    #Extract hit names from hits_output_filename
    hit_name=Extract_Hit_Names(output_filename)
    
    #OUTPUT:Write the hit-query mapping to a file
    file_map = open('hit_query_map_'+output_filename+'.tsv', 'w')
    file_map_90 = open('hit_query_map_'+output_filename+'_HitCoverage90.tsv', 'w')
    file_map.write('Hit_ID\tQuery_ID\tPercent_Identity\tAlignment_length\tQuery_length\tHit_length\tHit_alignment_start\tHit_alignment_end\tHit_coverage\tQuery_Coverage\tEvalue\n')
    file_map_90.write('Hit_ID\tQuery_ID\tPercent_Identity\tAlignment_length\tQuery_length\tHit_length\tHit_alignment_start\tHit_alignment_end\tHit_coverage\tQuery_Coverage\tEvalue\n')
    hits_90 =[]  #List of hits with coverage >=90
    for key in hit_map:
#	num_entry=1
	write_key=True    #Keeps track of when to write key in file_map_90
	file_map.write(key)
	for entries in hit_map[key]:
#		if num_entry>1:
#			file_map.write('\t')
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
#		num_entry+=1 

    #OUTPUT: Write the hit names to a file
    file_names=open('hits_name_'+output_filename+'.tsv', 'w')
    file_names_90=open('hits_name_'+output_filename+'_HitCoverage90.tsv', 'w')
    for key in hit_name:
	if key in set(hits_90):
		write_names_90=True
		file_names_90.write(key)
	else:
		write_names_90=False
	file_names.write(key)
	for val in hit_name[key]:     #To detect unlikely case of more than one entry of key (in that case output will have multiple tab seperated names)
		file_names.write('\t'+val)
		if write_names_90:
			file_names_90.write('\t'+val)
	file_names.write('\n')
	if write_names_90:
		file_names_90.write('\n')

    logging.info("Hit-query mapping written to file: %s\n", 'hit_query_map_'+output_filename+'.tsv')
    logging.info("Hit-query mapping for hits with hit coverage >=90.0 written to file: %s\n", 'hit_query_map_'+output_filename+'_HitCoverage90.tsv')
    logging.info("Hit names written to file: %s\n", 'hits_name_'+output_filename+'.tsv')
    logging.info("Hit names for hits with hit coverage >=90.0 written to file: %s\n", 'hits_name_'+output_filename+'_HitCoverage90.tsv')

    #Delete temporary file
    os.remove('hits_temp')

