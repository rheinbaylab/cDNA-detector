import pysam
import pandas as pd
import numpy as np
import re
import os
import sys
import collections
import scipy
from scipy import stats
import statsmodels


def f_get_input_files(args, global_para):
    genome_list = ['hg19','hg38','mm10','hg19_nochr','hg38_nochr']
    gene_model_list = [
    os.path.join(global_para.script_path,'data/gene_model/hg19.saf'),
    os.path.join(global_para.script_path,'data/gene_model/hg38.saf'),
    os.path.join(global_para.script_path,'data/gene_model/mm10.saf'),
    os.path.join(global_para.script_path,'data/gene_model/hg19_nochr.saf'),
    os.path.join(global_para.script_path,'data/gene_model/hg38_nochr.saf')
    ]
    # global_para.script_path = script_path
    blastn_database_list = [
    os.path.join(global_para.script_path,'data/blastn_database/total_blast.hg19.fa'),
    os.path.join(global_para.script_path,'data/blastn_database/total_blast.hg38.fa'),
    os.path.join(global_para.script_path,'data/blastn_database/total_blast.mm10.fa'),
    os.path.join(global_para.script_path,'data/blastn_database/total_blast.hg19_nochr.fa'),
    os.path.join(global_para.script_path,'data/blastn_database/total_blast.hg38_nochr.fa')
    ]
    genome_human_list = ['hg19','hg38','hg19_nochr','hg38_nochr']
    human_source_known = os.path.join(global_para.script_path, 'data/source_inference_additional/table_known_source.tsv')
    gene_model_dict = {genome_list[i]:gene_model_list[i] for i in range(len(genome_list))}
    blastn_databases_dict = {genome_list[i]:blastn_database_list[i] for i in range(len(genome_list))}
    source_known_dict = {genome_human_list[i]:human_source_known for i in range(len(genome_human_list))}
    # gene_model file
    if args.gene_model in genome_list:
        out_gene_model = gene_model_dict[args.gene_model]
    else:
        out_gene_model = args.gene_model
    # blastn databases file
    if (args.blastn_database == "") and (args.gene_model in genome_list):
        out_blastn_database = blastn_databases_dict[args.gene_model]
    elif args.blastn_database != "":
        out_blastn_database = args.blastn_database
    else:
        out_blastn_database = os.path.join(global_para.script_path, 'data/blastn_database/total_blast.fa')
    # known source file
    if args.file_source_known == "" and (args.gene_model in genome_human_list):
        out_source_known = source_known_dict[args.gene_model]
    elif args.file_source_known !="":
        out_source_known = args.file_source_known
    else:
        out_source_known = ""
    return out_gene_model, out_blastn_database, out_source_known


