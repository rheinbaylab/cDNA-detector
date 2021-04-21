## import modules
import pysam
import pandas as pd
import numpy as np
import re
import os
import sys
import gffpandas.gffpandas as gffpd
from Bio import SeqIO
import logging 
pd.options.mode.chained_assignment = None  # default='warn', this is used to remove warning from pandas chain assignment
np.seterr(divide = 'ignore') 


try:
    from . import global_para
except ImportError:
    import global_para

try:
    from .module_import_gtf import *
except ImportError:
    from module_import_gtf import *



def convert_gff2saf(global_para):
    # 1. check if output directory exist
    if os.path.isdir(global_para.output_dir):
        pass
    else:
        os.makedirs(global_para.output_dir)

    # 2. generate output filenames
    global_para.out_gtf_file = filename_gff2saf(global_para.gtf_file,global_para.output_dir)
    global_para.out_log = re.sub('.saf','.prepare.log',global_para.out_gtf_file)

    # 3. set the recording log
    global logger
    logger = set_logger(global_para.out_log)
    logger.info("input file is %s"%(global_para.gtf_file))
    logger.info("output file is %s"%global_para.out_gtf_file)
    logger.info('begin ...............')

    # 4. convert annotation file into dataframe with correct format
    logger.info('loading annotation files: be patient, a little slow')

    if global_para.format == "gtf":
        df_gff_cds_expand = f_import_transform_gtf(global_para)
    elif global_para.format == "bed":
        df_gff_cds_expand = f_import_transform_bed(global_para)
    else:
        logger.error("wrong input format, please use correct one: gtf or bed")
        exit(1)

    logger.info('finished!!!')

    # calculate chromosome information
    list_chr_all = df_gff_cds_expand.seq_id.unique().tolist()
    list_chr_keep = get_chromosome_keep(list_chr_all, global_para.chr_exclude, global_para.chr_include)


    ## 5. expand the gff file and only select bestref/best annotation for exon regions;
    logger.info("filter chromosomes")
    df_gff_cds_expand_filter = filter_chromosome(df_gff_cds_expand, list_chr_keep)
    ## 6. add exon boundary information and seq information
    logger.info("extract exon boundary sequences")
    genome_dict               = f_read_genome_seq(global_para.genome_fasta)
    df_gff_cds_expand_filter_1 = f_exon_boundary_info(df_gff_cds_expand_filter, genome_dict)
    df_gff_cds_expand_filter_2 = f_exon_flank_seq(df_gff_cds_expand_filter_1, genome_dict)

    ## 7. output results
    df_gtf_unique             = f_gff_exon_to_saf(df_gff_cds_expand_filter_2)
    list_column_upper         = ["exon_flank_start20", "exon_flank_end20", "exon_boundary_start_nearseq20", "exon_boundary_end_nearseq20"]
    df_gtf_unique             = f_sequence_upper_column(df_gtf_unique, list_column_upper)
    df_gtf_unique.to_csv(global_para.out_gtf_file, sep = '\t', index = False)

    ## send messages to users
    n_line = len(df_gtf_unique)
    logger.info("gene model files haved been generated")
    logger.info("%d coding exons have been used"%n_line)
    logger.info("please use the following gene model file in the next subcommand detect:")
    logger.info("%s"%global_para.out_gtf_file)

if __name__  == '__main__':
    global_para = global_class
