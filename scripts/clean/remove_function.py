import pysam
import pandas as pd
import numpy as np
import re
import sys
import os
import logging
from datetime import datetime
import random

try:
	import global_para
except ImportError:
	from . import global_para

def bam_million(file_bam):
	bam_genome = pysam.AlignmentFile(file_bam,'rb')
	num_million = bam_genome.count(until_eof=True)/1000000
	bam_genome.close()
	return num_million



def bam_ispair():
    # define bam is paired or not;
	bam_genome = pysam.AlignmentFile(global_para.file_bam,'rb')
	for read in bam_genome:
		return(read.is_paired)
		break
		

def read_allreads(subdf_region, file_bam):
	readlist_raw = []
	readlist_pos_type = []
	readlist_region = []
	# read read list
	i = 0
	bam_genome = pysam.AlignmentFile(file_bam,'rb')
	for index, row in subdf_region.iterrows():
		i = i + 1;
		# logging.logger.info(i)
		tmp_chr = row['seqname']
		tmp_start = int(row['start']) 
		tmp_end = int(row['end'])
		# logging.logger.info(tmp_chr, tmp_start, tmp_end)
		tmp_region_name = tmp_chr + ':' + str(tmp_start) + '-' + str(tmp_end)
		# read all reads in region
		for read in bam_genome.fetch(tmp_chr,tmp_start,tmp_end):
			if not read.is_unmapped:
				readlist_raw.append(read)
				readlist_region.append(tmp_region_name)
				tmp_read_start = read.pos
				tmp_read_end = int(read.aend)
				# if read belonged to exon boundaries
				if tmp_read_start <= int(tmp_start) or tmp_read_end >= int(tmp_end):
					readlist_pos_type.append('edge')
				else:
					readlist_pos_type.append('between')
	bam_genome.close()
	return readlist_raw, readlist_pos_type, readlist_region


def read_nameinfo(readlist):
	# 1. give read name and if paired information
	namelist = [read.qname for read in readlist]
	isread1_list = [read.is_read1 for read in readlist]
	# df = pd.DataFrame(namelist, isread1_list)
	return(namelist, isread1_list)

def readlist_pair(readlist1, readlist2):
	# identify readlist which is paired
	# x, y = readlist_pair(gene_readlist_cDNA_merge, gene_readlist_uDNA_merge)
	readlist1_name, readlist1_isread1 = read_nameinfo(readlist1)
	readlist2_name, readlist2_isread1 = read_nameinfo(readlist2)
	tmp_list1_name_df = pd.DataFrame({"list1.name":readlist1_name, "list1.isread1":readlist1_isread1})
	tmp_list1_name_df['order_list1'] = np.array(tmp_list1_name_df.index).tolist();
	tmp_list2_name_df = pd.DataFrame({"list2.name":readlist2_name, "list2.isread1":readlist2_isread1})
	tmp_list2_name_df['order_list2'] = np.array(tmp_list2_name_df.index).tolist();
	tmp_list_name_df = tmp_list1_name_df.merge(tmp_list2_name_df, left_on = 'list1.name', right_on = 'list2.name')
	tmp_list_name_df_filter = tmp_list_name_df[tmp_list_name_df['list1.isread1']  != tmp_list_name_df['list2.isread1']]
	# remove
	readlist1_paired = [readlist1[i] for i in tmp_list_name_df_filter['order_list1'].tolist()]
	readlist2_paired = [readlist2[i] for i in tmp_list_name_df_filter['order_list2'].tolist()]
	return readlist1_paired, readlist2_paired

def diff_bam_readlist_out(filename1,filename_output, read_list2, all_df_region):
	read_list2_dict = dict()
	read_list2_region = dict()
	read_list2_seqname = [read.reference_name for read in read_list2]
	read_list2_seqname_unique = list(set(read_list2_seqname))
	for seqname in read_list2_seqname_unique:
		# logging.logger.info(seqname)
		read_list2_dict[seqname] = set([read for read in read_list2 if read.reference_name == seqname])
		tmp_region_list = [[all_df_region.iloc[i,]['start'] - 1000,all_df_region.iloc[i,]['end'] + 1000] for i in range(len(all_df_region)) if all_df_region.iloc[i,]['seqname'] == seqname]
		# merge regions
		tmp_region_list.sort(key=lambda interval: interval[0])
		merged = [tmp_region_list[0]]
		for current in tmp_region_list:
		    previous = merged[-1]
		    if current[0] <= previous[1]:
		        previous[1] = max(previous[1], current[1])
		    else:
		        merged.append(current)
		read_list2_region[seqname] = merged
	bam_genome1 = pysam.AlignmentFile(filename1,'rb')
	tmpfile = pysam.AlignmentFile(filename_output, "wb", template=bam_genome1)
	all_seqname_list = bam_genome1.references
	for seqname in all_seqname_list:
		print(seqname)
		if seqname in read_list2_seqname_unique:
			tmp_region_list = read_list2_region[seqname]
			tmp_read_list = list()
			for read1 in bam_genome1.fetch(seqname):
				if len(tmp_region_list) ==0:
					tmpfile.write(read1)
					continue
				if read1.pos < tmp_region_list[0][0]:
					tmpfile.write(read1)
				elif read1.pos >= tmp_region_list[0][0] and read1.pos <=tmp_region_list[0][1]:
					tmp_read_list.append(read1)
				else:
					tmp_region_list.remove(tmp_region_list[0])
					tmp_read_list_set = set(tmp_read_list)
					tmp_read_list_keep = list(tmp_read_list_set - read_list2_dict[seqname])
					tmp_read_list_keep = np.array(tmp_read_list_keep)[np.argsort([read.pos for read in tmp_read_list_keep])]
					read_list2_dict[seqname] = read_list2_dict[seqname] - tmp_read_list_set
					tmp_read_list = list()
					if len(tmp_read_list_keep) >0:
						for read in tmp_read_list_keep:
							tmpfile.write(read)
					if read1 not in read_list2_dict[seqname]:
						tmpfile.write(read1)
					else:
						read_list2_dict[seqname].remove(read1)
		else:
			for read in bam_genome1.fetch(seqname):
				tmpfile.write(read)
	tmpfile.close()
	bam_genome1.close()


def diff_list(list1, list2): 
    return (list(set(list1) - set(list2))) 


def check_bam_index(genome_bam_file):
    ## check index of bam file; if no, generate one.
    print('Checking index of input bam file')
    if os.path.exists(genome_bam_file + '.bai') or os.path.exists(re.sub('bam$','bai',genome_bam_file)):
        print('Index file exists')
    else:
        print('file is not indexed, now generating index')
        pysam.index(genome_bam_file)
        print('Index file created\n')
        return
        
def create_dir(output_dir):
	if not os.path.isdir(output_dir):
		print("directory %s is not exists, create it"%output_dir)
		os.makedirs(output_dir)
	else:
		pass


def get_output_filename_bam(file_bam, output_dir):
	file_bam_basename = os.path.basename(file_bam)
	file_bam_clean_basename = re.sub('.bam$','.clean.bam',file_bam_basename)
	filename_output = os.path.join(output_dir, file_bam_clean_basename)
	return filename_output

def get_output_filename_region(file_region, output_dir):
	file_region_basename = os.path.basename(file_region)
	file_region_clean_basename = file_region_basename + ".clean"
	filename_output = os.path.join(output_dir, file_region_clean_basename)
	return filename_output




def reads_remove_automatic(exon_readlist_raw,exon_readlist_pos_type,exon_readlist_cDNA, exon_readlist_gDNA,exon_readlist_uDNA):
	# logging.logger.info('clean method:automatic')
	exon_readlist_edge = [exon_readlist_raw[i] for i in range(len(exon_readlist_raw)) if exon_readlist_pos_type[i]=="edge"]
	exon_readlist_cDNA_edge = list(set(exon_readlist_cDNA).intersection(set(exon_readlist_edge)))
	exon_readlist_gDNA_edge = list(set(exon_readlist_gDNA).intersection(set(exon_readlist_edge)))
	exon_readlist_uDNA_edge = list(set(exon_readlist_uDNA).intersection(set(exon_readlist_edge)))
	n_cDNA_edge = len(exon_readlist_cDNA_edge)
	n_gDNA_edge = len(exon_readlist_gDNA_edge)
	n_uDNA_edge = len(exon_readlist_uDNA_edge)
	n_cDNA_exon = len(exon_readlist_cDNA)
	n_gDNA_exon = len(exon_readlist_gDNA)
	n_uDNA_exon = len(exon_readlist_uDNA)
	n_exon_raw = len(exon_readlist_raw)
	if n_gDNA_exon == 0:
		n_uDNA_remove = n_uDNA_exon
		exon_readlist_remove = exon_readlist_raw
	else:
		# if (n_gDNA_edge + n_cDNA_edge) - n_cDNA_exon == 0:
		if any([(n_gDNA_edge + n_cDNA_edge) - n_cDNA_exon == 0, n_gDNA_edge + n_cDNA_edge == 0]):
			n_uDNA_remove = 0
		else:
			n_uDNA_remove = round(n_exon_raw*n_cDNA_edge/(n_gDNA_edge + n_cDNA_edge) - n_cDNA_exon)
		if n_uDNA_remove <=0:
			exon_readlist_remove =  exon_readlist_cDNA
		elif n_uDNA_remove>=len(exon_readlist_uDNA):
			exon_readlist_remove = exon_readlist_uDNA + exon_readlist_cDNA
		else:
			exon_readlist_remove_uDNA = random.sample(exon_readlist_uDNA, k = n_uDNA_remove)
			exon_readlist_remove = exon_readlist_remove_uDNA + exon_readlist_cDNA
	n_exon_keep = n_exon_raw - len(exon_readlist_remove)
	return(exon_readlist_remove, n_exon_raw, n_exon_keep)


def reads_remove_rpm(exon_readlist_raw,exon_readlist_rpm_number,bl_gDNA_remove, exon_readlist_cDNA, exon_readlist_gDNA,exon_readlist_uDNA):
	# logging.logger.info('clean method: reads per million')
	n_cDNA_exon = len(exon_readlist_cDNA)
	n_gDNA_exon = len(exon_readlist_gDNA)
	n_uDNA_exon = len(exon_readlist_uDNA)
	n_exon_raw = len(exon_readlist_raw)
	n_exon_keep_all = exon_readlist_rpm_number
	n_uDNA_remove = n_exon_raw - n_exon_keep_all - n_cDNA_exon
	if n_uDNA_remove <=0:
		exon_readlist_remove = exon_readlist_cDNA
	elif n_uDNA_remove >n_uDNA_exon:
		if bl_gDNA_remove == True:
			n_gDNA_remove = n_uDNA_remove - n_uDNA_exon
			exon_readlist_remove_gDNA = random.sample(exon_readlist_gDNA, k = n_gDNA_remove)
			exon_readlist_remove = exon_readlist_cDNA + exon_readlist_uDNA + exon_readlist_remove_gDNA
		else:
				exon_readlist_remove = exon_readlist_cDNA + exon_readlist_uDNA
	else:
		exon_readlist_remove_uDNA = random.sample(exon_readlist_uDNA, k = n_uDNA_remove)
		exon_readlist_remove = exon_readlist_remove_uDNA + exon_readlist_cDNA
	return exon_readlist_remove

def reads_remove_fraction(exon_readlist_raw,exon_readlist_fraction,bl_gDNA_remove, exon_readlist_cDNA, exon_readlist_gDNA,exon_readlist_uDNA):
	n_cDNA_exon = len(exon_readlist_cDNA)
	n_gDNA_exon = len(exon_readlist_gDNA)
	n_uDNA_exon = len(exon_readlist_uDNA)
	n_exon_raw = len(exon_readlist_raw)
	n_exon_keep_all = round(exon_readlist_fraction * n_exon_raw)
	n_uDNA_remove = n_exon_raw - n_exon_keep_all - n_cDNA_exon
	if n_uDNA_remove <=0:
		exon_readlist_remove = exon_readlist_cDNA
	elif n_uDNA_remove >n_uDNA_exon:
		if bl_gDNA_remove == True:
			n_gDNA_remove = n_uDNA_remove - n_uDNA_exon
			exon_readlist_remove_gDNA = random.sample(exon_readlist_gDNA, k = n_gDNA_remove)
			exon_readlist_remove = exon_readlist_cDNA + exon_readlist_uDNA + exon_readlist_remove_gDNA
		else:
				exon_readlist_remove = exon_readlist_cDNA + exon_readlist_uDNA
	else:
		exon_readlist_remove_uDNA = random.sample(exon_readlist_uDNA, k = n_uDNA_remove)
		exon_readlist_remove = exon_readlist_remove_uDNA + exon_readlist_cDNA
	return exon_readlist_remove
