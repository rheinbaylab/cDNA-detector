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

try:
	from .set_log import *
except:
	from set_log import *

try:
	from .assign_read_pass import *
except:
	from assign_read_pass import *

try:
	import global_para
except ImportError:
	from . import global_para


def cdna_remove_region(global_para):
	file_region = global_para.file_region
	file_bam = global_para.file_bam
	bl_gDNA_remove = global_para.bl_gDNA_remove
	# bl_paired = global_para.bl_paired
	bl_paired = bam_ispair()
	clean_way = global_para.clean_way
	output_dir = global_para.output_dir
	file_bam_clean = global_para.file_bam_clean
	file_region_clean = global_para.file_region_clean
	create_dir(output_dir)
	set_logger(global_para.out_log)
	check_bam_index(file_bam)
	logging.logger.info("clean region file is: %s"%file_region_clean)
	logging.logger.info('uncleaned bam file is: %s'%file_bam)
	# check file and check directory
	# 1. decode whole region file
	all_df_region = pd.read_csv(file_region, sep="\t",header = 0)
	all_df_region.seqname = all_df_region.seqname.astype('str')
	all_df_region['tmp_region'] = all_df_region.apply(lambda x:str(x.seqname) + ":" + str(x.start) + "-" + str(x.end), axis = 1)


	# 2. class all reads in exon regions into a big dictionary
	all_dict_region = dict()
	for gene_name, gene_df_region in all_df_region.groupby('gene_name'):
		logging.logger.info("assign reads for gene: %s"%gene_name)
		gene_dict = dict()
		if bl_paired == True:
			gene_readlist_raw, gene_readlist_pos_type, gene_readlist_cDNA, gene_readlist_gDNA, gene_readlist_uDNA = cdna_gene_read_assign_paired_end(gene_df_region,file_bam)
		else:
			gene_readlist_raw, gene_readlist_pos_type, gene_readlist_cDNA, gene_readlist_gDNA, gene_readlist_uDNA = cdna_gene_read_assign_single_end(gene_df_region,file_bam)
		for col in gene_df_region:
			gene_dict[col] = gene_df_region[col].tolist()
		gene_dict['gene_readlist_raw'] = gene_readlist_raw
		gene_dict['gene_readlist_cDNA'] = gene_readlist_cDNA
		gene_dict['gene_readlist_gDNA'] = gene_readlist_gDNA
		gene_dict['gene_readlist_uDNA'] = gene_readlist_uDNA
		gene_dict['gene_readlist_pos_type'] = gene_readlist_pos_type
		all_dict_region[gene_name] = gene_dict


	# 3. remove all reads in exon regions
	all_readlist_remove = list();
	num_bam_M = bam_million(file_bam)
	for gene_name in all_dict_region:
		gene_dict = all_dict_region[gene_name]
		if clean_way == "automatic":
			gene_dict['rpm'] = [float(0)] * len(gene_dict['seqname'])
			gene_dict['fraction'] = [float(0)] *len(gene_dict['seqname'])
		for i in range(len(gene_dict['gene_readlist_raw'])):
			exon_readlist_raw 	   = 	gene_dict['gene_readlist_raw'][i]
			exon_readlist_pos_type = 	gene_dict['gene_readlist_pos_type'][i]
			exon_readlist_cDNA     = 	gene_dict['gene_readlist_cDNA'][i]
			exon_readlist_gDNA     = 	gene_dict['gene_readlist_gDNA'][i]
			exon_readlist_uDNA     = 	gene_dict['gene_readlist_uDNA'][i]
			if clean_way == "fraction":
				exon_readlist_fraction = float(gene_dict['fraction'][i])
				exon_readlist_remove = reads_remove_fraction(exon_readlist_raw,exon_readlist_fraction,bl_gDNA_remove, exon_readlist_cDNA, exon_readlist_gDNA,exon_readlist_uDNA)
			elif clean_way == "rpm":
				exon_readlist_rpm_number = round(num_bam_M * float(gene_dict['rpm'][i]))
				exon_readlist_remove = reads_remove_rpm(exon_readlist_raw,exon_readlist_rpm_number,bl_gDNA_remove, exon_readlist_cDNA, exon_readlist_gDNA,exon_readlist_uDNA)
			elif  clean_way == "automatic":
				exon_readlist_remove, n_exon_raw, n_exon_keep = reads_remove_automatic(exon_readlist_raw,exon_readlist_pos_type,exon_readlist_cDNA, exon_readlist_gDNA,exon_readlist_uDNA)
				value_rpm_exon = n_exon_keep/num_bam_M
				if n_exon_raw == 0:
					value_fraction_exon = 0
				else:
					value_fraction_exon = n_exon_keep/n_exon_raw
				gene_dict['rpm'][i] = value_rpm_exon
				gene_dict['fraction'][i] = value_fraction_exon

			all_readlist_remove.extend(exon_readlist_remove)

	logging.logger.info('generate clean bam file')
	diff_bam_readlist_out(file_bam,file_bam_clean, all_readlist_remove, all_df_region)
	pysam.index(file_bam_clean)


	if clean_way == "automatic":
		logging.logger.info("generate statistic file")
		columns_list = ['seqname', 'start', 'end', 'gene_name','rpm', 'fraction']
		all_df_region = all_df_region.filter(columns_list)
		for columns in columns_list:
			x = list()
			[x.extend(all_dict_region[gene][columns]) for gene in all_dict_region.keys()]
			all_df_region[columns] = x
		all_df_region.to_csv(file_region_clean, header = True,sep = '\t', index = False)


	logging.logger.info("Program finished successfully")

class para_process:
	def __init__(self, args):
		self.file_region = args.region
		self.file_bam = args.bam
		self.bl_gDNA_remove = args.gDNA_remove
		self.bl_paired = args.paired
		self.clean_way = args.clean_way
		self.output_dir = args.output_dir
		self.file_bam_clean = file_bam_clean
		self.file_region_clean = file_region_clean
		self.file_bam_clean = get_output_filename_bam(self.file_bam, self.output_dir)
		self.file_region_clean = get_output_filename_region(self.file_region, self.output_dir)


if __name__ == '__main__':
	class global_para:
		pass
	cdna_remove_region(global_para)





