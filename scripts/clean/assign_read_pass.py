import pysam
import pandas as pd
import numpy as np
import re
import sys
import random
import os
import logging
from itertools import compress
from datetime import datetime

try:
	from .remove_function import *
except:
	from remove_function import *

def read_cDNA(exon_readlist_raw, exon_readlist_pos_type, tmp_chr, tmp_start, tmp_end):
	# logging.logger.info('define if reads belong to cDNA')
	# 1. start, start_S
	# 2. end, end_S;
	# 3. start, start_consensus_seq
	# 4. end, end_consequence
	# begin
	readlist_cDNA = list()
	# logging.logger.info('define if reads belong to cDNA')
	# 1. start, start_S
	len_list= len(exon_readlist_raw);
	list_nuc_start = list();
	for i in range(len_list):
		tmp_read = exon_readlist_raw[i]
		tmp_cigar = tmp_read.cigar[0][0]
		tmp_cigar_len = tmp_read.cigar[0][1]
		tmp_seq = tmp_read.to_dict()['seq']
		# if re.search("^[0-9]*S",tmp_cigar) and exon_readlist_raw[i].pos == tmp_start:
		if tmp_cigar == 4 and tmp_read.pos == tmp_start:
			readlist_cDNA.append(tmp_read)
			list_nuc_start.append(tmp_seq[0:int(tmp_cigar_len)][::-1])
		elif tmp_cigar == 5 and tmp_read.pos == tmp_start:
			readlist_cDNA.append(tmp_read)
		# 2. start, start_consensus_seq
	consensus_seq_start = str();
	list_nuc_start.sort(key = len)
	if len(list_nuc_start)==0:
		pass;
	else:
		for i in range(0,len(list_nuc_start[-1])):
			tmp = [x[0] for x in list_nuc_start if x]
			tmp_count = [ tmp.count(n)  for n in list('ACGT') ]
			tmp_sum = sum(tmp_count)+ 0.01
			tmp_freq = [float(x)/tmp_sum for x in tmp_count]
			list_nuc_start = [x[1:] for x in list_nuc_start if x]
			if tmp_count[tmp_freq.index(max(tmp_freq))]>=2 and max(tmp_freq)>=0.8:
				consensus_seq_start = consensus_seq_start +  ['ACGT'][0][tmp_freq.index(max(tmp_freq))]
			else:
				break
		for read in exon_readlist_raw:
			if read.pos < tmp_start:
				tmp_seq = read.to_dict()['seq']   
				tmp_seq_hang = tmp_seq[0:(len(tmp_seq) - len([tmp_pos for tmp_pos in read.positions if tmp_pos>=tmp_start]))][::-1]
				if any([re.match(tmp_seq_hang,consensus_seq_start),re.match(consensus_seq_start,tmp_seq_hang)]) and len(tmp_seq_hang)>0:
					readlist_cDNA.append(read)
		# 2. end, end_S;
	list_nuc_end = list();
	for i in range(len_list):
		tmp_read = exon_readlist_raw[i]
		tmp_cigar = tmp_read.cigar[-1][0]
		tmp_cigar_len = tmp_read.cigar[-1][1]
		tmp_seq = tmp_read.to_dict()['seq']
		if tmp_cigar == 4 and tmp_read.aend == tmp_end:
			readlist_cDNA.append(tmp_read)
			list_nuc_end.append(tmp_seq[-int(tmp_cigar_len):])
		elif tmp_cigar == 5 and tmp_read.aend == tmp_end:
			readlist_cDNA.append(tmp_read)
		# 4. end, end_consequence
	consensus_seq_end = str();
	list_nuc_end.sort(key = len)
	if len(list_nuc_end)==0:
		pass
	else:
		for i in range(0,len(list_nuc_end[-1])):
			tmp = [x[0] for x in list_nuc_end if x]
			tmp_count = [ tmp.count(n) for n in list('ACGT') ]
			tmp_sum = sum(tmp_count)+ 0.01
			tmp_freq = [float(x)/tmp_sum for x in tmp_count]
			list_nuc_end = [x[1:] for x in list_nuc_end if x]
			if tmp_count[tmp_freq.index(max(tmp_freq))]>=2 and max(tmp_freq)>=0.8:
				consensus_seq_end = consensus_seq_end +  ['ACGT'][0][tmp_freq.index(max(tmp_freq))]
			else:
				break
		for read in exon_readlist_raw:
			if read.aend>tmp_end:
				tmp_seq = read.to_dict()['seq']
				tmp_seq_hang = tmp_seq[::-1][0:int(len(tmp_seq) - len([tmp_pos for tmp_pos in read.positions if tmp_pos<tmp_end]))][::-1]
				if any([re.match(tmp_seq_hang,consensus_seq_end),re.match(consensus_seq_end,tmp_seq_hang)]) and len(tmp_seq_hang)>0:
					readlist_cDNA.append(read)
	readlist_cDNA = list(set(readlist_cDNA))
	return(readlist_cDNA)


def read_gDNA(exon_readlist_raw, exon_readlist_pos_type, tmp_chr, tmp_start, tmp_end,exon_readlist_cDNA):
	# logging.logger.info('read_gDNA')
	# 1. start, span the exon regions;
	# 2. end, span the exon regions;
	readlist_gDNA = list();
	len_list= len(exon_readlist_raw);
	for i in range(len_list):
		# logging.logger.info(i)
		tmp_read = exon_readlist_raw[i]
		tmp_read_pos_start = tmp_read.pos
		tmp_read_pos_end = tmp_read.aend
		if exon_readlist_pos_type[i] == 'edge':
			if (int(tmp_read_pos_start) - int(tmp_start))*(tmp_read_pos_end- int(tmp_start)) < 0:
				readlist_gDNA.append(tmp_read)
			elif (int(tmp_read_pos_start) - int(tmp_end))*(tmp_read_pos_end - int(tmp_end))< 0:
				readlist_gDNA.append(tmp_read)	
	readlist_gDNA = list(set(readlist_gDNA) - set(exon_readlist_cDNA))
	
	return readlist_gDNA


def read_unclass(readlist_raw, readlist_cDNA, readlist_gDNA):
	# logging.logger.info('read_unclass')
	readlist_uDNA = list(set(readlist_raw) - set(readlist_cDNA) - set(readlist_gDNA))
	return readlist_uDNA




def cdna_exon_read_assign(tmp_region, exon_df_region, file_bam):
	tmp_chr = exon_df_region.iloc[0]['seqname']
	tmp_start = exon_df_region.iloc[0]['start']
	tmp_end = exon_df_region.iloc[0]['end']
	exon_readlist_raw, exon_readlist_pos_type, exon_readlist_region = read_allreads(exon_df_region, file_bam)
	## class the raw read list from exon region
	exon_readlist_cDNA = read_cDNA(exon_readlist_raw, exon_readlist_pos_type, tmp_chr, tmp_start, tmp_end)
	exon_readlist_gDNA = read_gDNA(exon_readlist_raw, exon_readlist_pos_type, tmp_chr, tmp_start, tmp_end,exon_readlist_cDNA)
	exon_readlist_uDNA = read_unclass(exon_readlist_raw, exon_readlist_cDNA, exon_readlist_gDNA)
	return exon_readlist_raw, exon_readlist_pos_type, exon_readlist_region, exon_readlist_cDNA, exon_readlist_gDNA, exon_readlist_uDNA




def cdna_gene_read_assign_single_end(gene_df_region,file_bam):
	gene_readlist_raw = list()
	gene_readlist_pos_type = list()
	gene_readlist_region = list()
	gene_readlist_cDNA = list()
	gene_readlist_gDNA = list()
	gene_readlist_uDNA = list()
	gene_readlist_raw_merge = list()
	gene_readlist_pos_type_merge = list()
	gene_readlist_region_merge = list()
	gene_readlist_cDNA_merge = list()
	gene_readlist_gDNA_merge = list()
	gene_readlist_uDNA_merge = list()
	for tmp_region, exon_df_region in gene_df_region.groupby('tmp_region'):
		exon_readlist_raw, exon_readlist_pos_type, exon_readlist_region, exon_readlist_cDNA, exon_readlist_gDNA, exon_readlist_uDNA = cdna_exon_read_assign(tmp_region, exon_df_region, file_bam)
		gene_readlist_raw.append(exon_readlist_raw) 
		gene_readlist_pos_type.append(exon_readlist_pos_type)
		gene_readlist_region.append(exon_readlist_region)
		gene_readlist_cDNA.append(exon_readlist_cDNA)
		gene_readlist_gDNA.append(exon_readlist_gDNA)
		gene_readlist_uDNA.append(exon_readlist_uDNA)
	return(gene_readlist_raw, gene_readlist_pos_type, gene_readlist_cDNA, gene_readlist_gDNA, gene_readlist_uDNA)



def cdna_gene_read_assign_paired_end(gene_df_region, file_bam):
	gene_readlist_raw = list()
	gene_readlist_pos_type = list()
	gene_readlist_region = list()
	gene_readlist_cDNA = list()
	gene_readlist_gDNA = list()
	gene_readlist_uDNA = list()
	gene_readlist_raw_merge = list()
	gene_readlist_pos_type_merge = list()
	gene_readlist_region_merge = list()
	gene_readlist_cDNA_merge = list()
	gene_readlist_gDNA_merge = list()
	gene_readlist_uDNA_merge = list()
	for tmp_region, exon_df_region in gene_df_region.groupby('tmp_region'):
		exon_readlist_raw, exon_readlist_pos_type, exon_readlist_region, exon_readlist_cDNA, exon_readlist_gDNA, exon_readlist_uDNA = cdna_exon_read_assign(tmp_region, exon_df_region, file_bam)
		gene_readlist_raw.append(exon_readlist_raw) 
		gene_readlist_pos_type.append(exon_readlist_pos_type)
		gene_readlist_region.append(exon_readlist_region)
		gene_readlist_cDNA.append(exon_readlist_cDNA)
		gene_readlist_gDNA.append(exon_readlist_gDNA)
		gene_readlist_uDNA.append(exon_readlist_uDNA)
		# merge all readlist information
		gene_readlist_raw_merge.extend(exon_readlist_raw)
		gene_readlist_pos_type_merge.extend(exon_readlist_pos_type)
		gene_readlist_region_merge.extend(exon_readlist_region)
		gene_readlist_cDNA_merge.extend(exon_readlist_cDNA)
		gene_readlist_gDNA_merge.extend(exon_readlist_gDNA)
		gene_readlist_uDNA_merge.extend(exon_readlist_uDNA)		
	gene_readlist_cDNA_merge_tmp_twoexon = list()
	for i in range(len(gene_readlist_raw)):
		tmp_gene_readlist_uDNA_list1 = gene_readlist_uDNA[i]
		tmp_gene_readlist_uDNA_list2 = list(set(exon_readlist_uDNA) - set(tmp_gene_readlist_uDNA_list1))
		tmp_gene_readlist1_cDNA, tmp_gene_readlist2_cDNA = readlist_pair(tmp_gene_readlist_uDNA_list1, tmp_gene_readlist_uDNA_list2)
		tmp_gene_readlist1_cDNA.extend(tmp_gene_readlist2_cDNA)
		gene_readlist_cDNA_merge_tmp_twoexon.extend(list(set(tmp_gene_readlist1_cDNA) - set(gene_readlist_gDNA_merge)))
	gene_readlist_cDNA_merge_tmp_twoexon = list(set(gene_readlist_cDNA_merge_tmp_twoexon))
	gene_readlist_uDNA_merge = diff_list(gene_readlist_uDNA_merge, gene_readlist_cDNA_merge_tmp_twoexon)
	# 2. one of reads are belonged to cds reads
	gene_readlist_cDNA_merge_tmp_1cdna = list()
	tmp_gene_readlist_uDNA_list1 = gene_readlist_uDNA_merge
	tmp_gene_readlist_uDNA_list2 = list(set(gene_readlist_cDNA_merge_tmp_twoexon + gene_readlist_cDNA_merge ))
	tmp_gene_readlist1_cDNA, tmp_gene_readlist2_cDNA = readlist_pair(tmp_gene_readlist_uDNA_list1, tmp_gene_readlist_uDNA_list2)
	gene_readlist_cDNA_merge_tmp_1cdna = list(set(tmp_gene_readlist1_cDNA))
	gene_readlist_cDNA_merge = list(set(gene_readlist_cDNA_merge + gene_readlist_cDNA_merge_tmp_1cdna + gene_readlist_cDNA_merge_tmp_twoexon))
	gene_readlist_uDNA_merge = diff_list(gene_readlist_uDNA_merge, gene_readlist_cDNA_merge_tmp_1cdna)
	# 2. gDNA
	# one of reads are located outside exon regions
	tmp_gene_readlist_uDNA_list1 = list(compress(gene_readlist_uDNA_merge,[read.is_proper_pair for read in gene_readlist_uDNA_merge]))
	tmp_gene_readlist_uDNA_list2 = gene_readlist_raw_merge
	tmp_gene_readlist1_uDNA, tmp_gene_readlist2_uDNA = readlist_pair(tmp_gene_readlist_uDNA_list1, tmp_gene_readlist_uDNA_list2)
	gene_readlist_gDNA_merge_tmp_1noexon = diff_list(tmp_gene_readlist_uDNA_list1, tmp_gene_readlist1_uDNA)
	# one of reads are belonged to gDNA
	tmp_gene_readlist_uDNA_list1 = list(compress(gene_readlist_uDNA_merge,[read.is_proper_pair for read in gene_readlist_uDNA_merge]))
	tmp_gene_readlist_uDNA_list2 = gene_readlist_gDNA_merge
	tmp_gene_readlist1_cDNA, tmp_gene_readlist2_cDNA = readlist_pair(tmp_gene_readlist_uDNA_list1, tmp_gene_readlist_uDNA_list2)
	gene_readlist_gDNA_merge_tmp_1gDNA = tmp_gene_readlist1_cDNA
	gene_readlist_cDNA_merge_new = list(set(gene_readlist_cDNA_merge + gene_readlist_cDNA_merge_tmp_1cdna + gene_readlist_cDNA_merge_tmp_twoexon ))
	gene_readlist_gDNA_merge_new = list(set(gene_readlist_gDNA_merge +  gene_readlist_gDNA_merge_tmp_1gDNA + gene_readlist_gDNA_merge_tmp_1noexon))
	gene_readlist_uDNA_merge_new = diff_list(gene_readlist_uDNA_merge, gene_readlist_cDNA_merge_new + gene_readlist_gDNA_merge_new)
	# return gene_readlist_cDNA_merge_new, gene_readlist_gDNA_merge_new, gene_readlist_uDNA_merge_new
	# reassign genelist to results
	gene_readlist_cDNA = [list.clear() for list in gene_readlist_cDNA]
	gene_readlist_gDNA = [list.clear() for list in gene_readlist_gDNA]
	gene_readlist_uDNA = [list.clear() for list in gene_readlist_uDNA]
	for i in range(len(gene_readlist_raw)):
		tmp_gene_readlist_raw = set(gene_readlist_raw[i])
		gene_readlist_cDNA[i] = list(set(gene_readlist_cDNA_merge_new).intersection(tmp_gene_readlist_raw))
		gene_readlist_gDNA[i] = list(set(gene_readlist_gDNA_merge_new).intersection(tmp_gene_readlist_raw))
		gene_readlist_uDNA[i] = list(set(gene_readlist_uDNA_merge_new).intersection(tmp_gene_readlist_raw))
	return(gene_readlist_raw,gene_readlist_pos_type, gene_readlist_cDNA, gene_readlist_gDNA, gene_readlist_uDNA)

