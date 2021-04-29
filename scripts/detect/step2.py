#!/usr/bin/env python
# coding: utf-8

# import modules
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
import subprocess
import logging
from statsmodels.stats.multitest import fdrcorrection
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Seq


pd.options.mode.chained_assignment = None  # default='warn', this is used to remove warning from pandas chain assignment
np.seterr(divide = 'ignore') 
# np.seterr(divide = 'warn') 

try:
    from . import global_para
except ImportError:
    import global_para

try:
    from .consensus_seq import *
except ImportError:
    from consensus_seq import *

try:
    from .math_stat import *
except ImportError:
    from math_stat import *

try:
    from .unaligned_reads import *
except ImportError:
    from unaligned_reads import *

try:
    from .multi_processor import *
except ImportError:
    from multi_processor import *

try:
    from .set_log import *
except ImportError:
    from set_log import *

try:
    from .clipped_seq_source_inference import *
except ImportError:
    from clipped_seq_source_inference import *



def detect_cdna(global_para):
    # 1. preparation for loading files and directory
    ### check index file, if no, create it.
    check_bam_index(global_para.genome_bam_file)
    global_para.genome_bam_paired = bam_ispair()
    ### check output directory exist, if not, create it
    if os.path.isdir(global_para.output_dir):
        pass
    else:
        os.makedirs(global_para.output_dir)
    ### create logger file
    global_para.logger = set_logger()
    global_para.logger.info('Begin to process sample: %s'%global_para.genome_bam_file)

    # 2. load gene model
    ### load transcript-based gene model file
    df_transcript_exon = read_gene_model(global_para.gtf_gene_unique_file)
    df_transcript_exon = df_transcript_exon.drop_duplicates()
    global_para.exon_distance = 10
    df_transcript_exon = f_close_exon_merge(df_transcript_exon)
    df_transcript_exon['region'] = f_df_add_region(df_transcript_exon)
    ### create gene-based exon file
    list_attribute_gene = ['seqname','start','end','gene_id','gene_name','region']
    df_gene_exon = df_transcript_exon.filter(list_attribute_gene).drop_duplicates()
    ### check if chromosome names are consistent between bam and gene model
    df_gene_exon = f_overlap_reference(global_para.genome_bam_file,df_gene_exon)



    # 3. scan whole exon regions and get a initial selection for possible cDNAs
    ### split gene model and scan it 
    global_para.logger.info('Start scanning exon regions...')
    df_gene_exon = f_exon_group_by_distance(df_gene_exon)
    x = apply_parallel(df_gene_exon.groupby('read_group'),f_find_unalign_readlist_multi)
    ### get background statistic information
    global_para.num_total_read_exon = sum([i[0] for i in x]) 
    global_para.num_unalign_read = sum([i[1] for i in x])
    ### get information of each exons
    df_gene_exon['num_unalign_read'] = [i for sublist in x for i in sublist[2]]
    df_gene_exon['num_unalign_read_exon_edge'] = [i for sublist in x for i in sublist[3]]
    df_gene_exon['num_all_read_exon_edge'] = [i for sublist in x for i in sublist[4]]
    ### calculate ratio of unaligned reads and remove extremely high coverage unaligned reads from background information
    p_unalign_read = f_get_ratio_unaligned(global_para)
    if p_unalign_read >= global_para.ratio_ehc and  global_para.exclude_ehc == True: 
        f_adjust_bg_unaligned_reads(df_gene_exon)
    ### set backgrounad information to df, and calculate intial p-value for exons.
    df_gene_exon['num_bg_total_read_exon'] = global_para.num_total_read_exon
    df_gene_exon['num_bg_unalign_read_exon'] = global_para.num_unalign_read
    print("total reads in exon regions is: %d "%global_para.num_total_read_exon)
    print("clipped reads in exon regions is: %d "%global_para.num_unalign_read)
    df_gene_exon['pvalue'] = df_gene_exon.apply(lambda x:beta_binomial_significance(x.num_unalign_read_exon_edge, x.num_all_read_exon_edge, x.num_bg_unalign_read_exon, x.num_bg_total_read_exon),axis = 1)
    ### filter specific gene list
    tmp_gene_name_filter_step1 = f_filter_genelist_step1(df_transcript_exon, df_gene_exon)
    ### get short transcript-exon and gene-exon information, and add column 'region' for connection between transcript-exon and gene-exon
    df_transcript_exon_1 = df_transcript_exon.query('gene_name in @tmp_gene_name_filter_step1')
    df_gene_exon_1 = df_gene_exon.query('gene_name in @tmp_gene_name_filter_step1')
    df_region = f_df_collapse_seq_boundary_info(df_transcript_exon_1)
    df_gene_exon_1 = df_gene_exon_1.merge(df_region, on = 'region')
    df_gene_exon_filter = df_gene_exon_1.copy()


    # 4. stat reads around exon edges for each regions.
    global_para.logger.info('Evaluating %d potential candidate cDNA(s) ......'%len(df_gene_exon_filter.gene_name.unique()))
    list_stat_all_region_fdr = list()
    for gene, df_gene in df_gene_exon_filter.groupby(['gene_name']):
        print(gene)
        global_para.logger.info(gene)
        df_gene = df_gene_exon_filter.query('gene_name==@gene')
        list_stat_start = list()
        list_stat_end = list();
        ### recal information of each boundary
        for region, df_exon in df_gene.groupby('region'):
            list_stat_start.append(f_df_exon_start_stat(df_exon))
            list_stat_end.append(f_df_exon_end_stat(df_exon))
        df_stat_start = pd.concat(list_stat_start)
        df_stat_end = pd.concat(list_stat_end)
        if len(df_stat_start)==0 or len(df_stat_end)==0:
            continue
        ### assign specific clipped position to each boundary
        df_stat_start = f_df_gene_start_stat_remove_dup(df_stat_start)
        df_stat_end = f_df_gene_end_stat_remove_dup(df_stat_end)
        ### merge reads in start/end to whole exons
        df_gene = df_gene.merge(df_stat_start, left_on = "region",right_on = "region_start" ,suffixes = ("","_start"), how = "left")
        df_gene = df_gene.merge(df_stat_end, left_on = "region",right_on = "region_end", suffixes = ("","_end"), how = "left")
        ### remove duplication information
        df_gene = f_df_gene_merge_removedup(df_gene)
        if len(df_gene)==0:
            continue
        ### adjust most-left position
        df_exon_start = df_gene.query('is_exon_boundary_start.str.contains("1")', engine = 'python')
        for tmp_index_start in df_exon_start.index.tolist():
            tmp_df_exon_start = df_gene.loc[tmp_index_start]
            if tmp_df_exon_start.is_exon_match_end and tmp_df_exon_start.bbinom_p_end <= global_para.cutoff_pvalue and tmp_df_exon_start.bbinom_p_start >=0.05:
                df_gene.loc[tmp_index_start] = exon_start_adj(tmp_df_exon_start)
        ### adjust for exon end position
        df_exon_end = df_gene.query('is_exon_boundary_end.str.contains("1")', engine = 'python')
        for tmp_index_end in df_exon_end.index.tolist():
            tmp_df_exon_end = df_gene.loc[tmp_index_end]
            if tmp_df_exon_end.is_exon_match_start and tmp_df_exon_end.bbinom_p_start <= global_para.cutoff_pvalue and tmp_df_exon_end.bbinom_p_end >=0.05:
                df_gene.loc[tmp_index_end] = exon_end_adj(tmp_df_exon_end)
        ### generate clipped sequences for gene
        df_gene['combine_pvalue'] = df_gene.apply(lambda x:f_combin_p([x.bbinom_p_start,x.bbinom_p_end]),axis = 1)
        df_gene['consensus_seq_start'] = df_gene.consensus_seq_start.fillna('')
        df_gene['consensus_seq_end'] = df_gene.consensus_seq_end.fillna('')
        df_gene['clipped_seq_start'] = df_gene.apply(lambda x:x.consensus_seq_start[::-1],axis = 1)
        df_gene['clipped_seq_end'] = df_gene.consensus_seq_end
        list_stat_all_region_fdr.append(df_gene)


    f_if_0cdna(list_stat_all_region_fdr)
    df_gene_fdr = pd.concat(list_stat_all_region_fdr)
    df_gene_fdr['combine_qvalue'] = fdrcorrection(df_gene_fdr['combine_pvalue'])[1]


    # 5. utilize different information to define if exon might be cDNA-enriched
    list_stat_all_region = []
    for gene_name, df_gene in df_gene_fdr.groupby('gene_name'):
        df_gene['is_exon_keep'] = df_gene.apply(lambda x:f_exon_filter(x, global_para),axis = 1)
        df_gene_filter = df_gene.query('is_exon_keep==True')
        if len(df_gene_filter)>=1:
            tmp1= df_gene.query('is_exon_boundary_start=="1"').query('is_exon_match_end==True').query('bbinom_p_end<@global_para.cutoff_pvalue')
            tmp2= df_gene.query('is_exon_boundary_end=="1"').query('is_exon_match_start==True').query('bbinom_p_start<@global_para.cutoff_pvalue')
            df_gene_filter = pd.concat([df_gene_filter,tmp1,tmp2])
            df_gene_filter = df_gene_filter.drop_duplicates()
        list_stat_all_region.append(df_gene_filter)

    ### merge all info for all genes
    f_if_0cdna(list_stat_all_region)
    df_stat_all_region = pd.concat(list_stat_all_region)
    f_if_0cdna(df_stat_all_region)


    # 6. generate output files
    ### 6.1 exon stat
    df_stat_all_region.reset_index(inplace = True,drop = True)
    df_stat_all_region_full = df_stat_all_region.copy()
    df_stat_all_region = df_stat_all_region.rename(columns = {'num_bad_reads_ajust_start':'num_clipped_start', 'num_all_reads_ajust_start':'num_total_start', 'bbinom_p_start':'bbinom_pvalue_start', 'num_bad_reads_ajust_end':'num_clipped_end', 'num_all_reads_ajust_end':'num_total_end', 'bbinom_p_end':'bbinom_pvalue_end', 'bg_unalign_start':'num_bg_clipped', 'bg_total_start':'num_bg_total'})
    select_columns_rename = ['seqname', 'exon_start', 'exon_end', 'pos_start','pos_end','gene_id', 'gene_name','transcript_id','num_clipped_start', 'num_total_start',  'bbinom_pvalue_start','num_clipped_end', 'num_total_end',  'bbinom_pvalue_end','num_bg_clipped', 'num_bg_total','combine_pvalue', 'combine_qvalue']
    # df_stat_all_region = df_stat_all_region.astype({'pos_start':int,'pos_end':int,'num_clipped_start':int,'num_total_start':int,'num_clipped_end':int,'num_total_end':int,'num_bg_clipped':int,'num_bg_total':int}, errors = 'ignore')
    df_stat_all_region = df_stat_all_region.astype({'exon_start':"Int64",'exon_end':"Int64", 'pos_start':"Int64",'pos_end':'Int64','num_clipped_start':'Int64','num_total_start':'Int64','num_clipped_end':'Int64','num_total_end':'Int64','num_bg_clipped':'Int64','num_bg_total':'Int64'})
    df_stat_all_region = df_stat_all_region.filter(select_columns_rename)
    df_stat_all_region['exon_start'] = df_stat_all_region['exon_start'] -1
    df_stat_all_region['pos_start'] = df_stat_all_region['pos_start'] -1
    df_stat_all_region.to_csv(global_para.out_exon_stat, sep = "\t", index = False)



    ## 6.2 gene stat for each transcript; then select most possible one. then get the source of each transcript information.
    df_stat_all_transcript = f_stat_most_possible_transcript(df_stat_all_region, df_gene_exon_1, df_gene_fdr, df_gene_exon_filter)

    ## recheck the exon;
    df_stat_all_transcript['num_exon_detected'] = df_stat_all_transcript.apply(lambda x:f_recheck_exon(x,df_stat_all_region_full),axis =1 )
    df_stat_all_transcript['ratio'] = df_stat_all_transcript.apply(lambda x:x.num_exon_detected/x.num_exon_transcript, axis = 1)
    df_stat_all_transcript['combined_pvalue'] = df_stat_all_transcript.apply(lambda x:f_combine_pvalue_transcript(x.transcript_id, x.gene_name, x.num_exon_transcript, df_stat_all_region), axis = 1)
    df_stat_all_transcript['avg_log10pvalue'] = df_stat_all_transcript.apply(lambda x:f_avg_minus_log10P(x.transcript_id, x.gene_name, x.num_exon_transcript, df_stat_all_region), axis = 1)
    df_stat_all_transcript['avg_pvalue'] = df_stat_all_transcript.apply(lambda x:f_avg_pvalue(x.transcript_id, x.gene_name, x.num_exon_transcript, df_stat_all_region), axis = 1)
    df_stat_all_transcript['avg_cDNA'] = df_stat_all_transcript.apply(lambda x:f_avg_cDNA(x.transcript_id, x.gene_name, x.num_exon_transcript, df_stat_all_region), axis = 1)
    df_stat_all_transcript['median_cDNA'] = df_stat_all_transcript.apply(lambda x:f_median_cDNA(x.transcript_id, x.gene_name, x.num_exon_transcript, df_stat_all_region), axis = 1)




    # get source information
    list_select_transcript = df_stat_all_transcript.transcript_id.unique().tolist()
    df_stat_all_region_full_reinfo = df_stat_all_region_full.merge(df_transcript_exon_1, on = "region", suffixes = ['_1',""]).query('transcript_id in @list_select_transcript')
    table_blastn = f_get_blast_4gene(df_stat_all_region_full_reinfo)
    table_blastn_source = f_source_inference(df_stat_all_region_full_reinfo, table_blastn)
    table_source_known = f_df_source_known(global_para.file_source_known)
    table_blastn_source['source_known_databases'] = table_blastn_source.apply(lambda x:f_source_known(x.gene_name, table_source_known), axis = 1)
    table_blastn_source.to_csv(global_para.out_blastn_seq_source,sep = '\t', index = False)
    df_table_blastn_source = table_blastn_source.filter(['gene_name','source_inference','source_known_databases'])


    df_stat_all_transcript = df_stat_all_transcript.merge(df_table_blastn_source, how = 'left', on = "gene_name")
    df_stat_all_transcript = df_stat_all_transcript.sort_values(by = ['median_cDNA','avg_pvalue'], ascending=[False, True])
    df_stat_all_transcript.to_csv(global_para.out_gene_stat,index = False, sep = "\t")



    ## 6.3 filter gene and exon information
    df_stat_all_transcript_filter = df_stat_all_transcript.query('ratio>=@global_para.cutoff_ratio_gene').query('median_cDNA>=@global_para.cutoff_num_exon_unaligned_reads')
    f_if_0cdna(df_stat_all_transcript_filter)
    list_transcript_filter = df_stat_all_transcript_filter.transcript_id.unique().tolist()
    df_stat_all_region_filter = df_stat_all_region[df_stat_all_region.apply(lambda x:str_list_compare(list_transcript_filter, x.transcript_id),axis = 1)]
    df_stat_all_region_filter.to_csv(global_para.out_exon_stat_filter, sep = '\t', index = False)
    df_stat_all_transcript_filter = df_stat_all_transcript_filter.sort_values(by = ['median_cDNA','avg_pvalue'], ascending=[False, True])
    df_stat_all_transcript_filter.to_csv(global_para.out_gene_stat_filter,index = False, sep = "\t")


    ## code for filter source inference
    f_if_0cdna(df_stat_all_transcript_filter)
    df_stat_all_transcript_filter_source_filter = f_filter_by_source(df_stat_all_transcript_filter)
    df_stat_all_transcript_filter_source_filter = df_stat_all_transcript_filter_source_filter.sort_values(by = ['median_cDNA','avg_pvalue'], ascending=[False, True])
    df_stat_all_transcript_filter_source_filter.to_csv(global_para.out_gene_stat_filter_source,index = False, sep = "\t")
    f_if_0cdna(df_stat_all_transcript_filter_source_filter)





    ## 6.4 generate region information: need to merge some regions
    select_columns = ['seqname','pos_start','pos_end','gene_id','gene_name','region']
    select_columns_rename = ['seqname','start','end','gene_id','gene_name','region','filter']
    df_stat_all_region_filter['region'] = df_stat_all_region_full.loc[df_stat_all_region_filter.index].region.tolist()
    list_gene_filter = df_stat_all_transcript_filter_source_filter.gene_name.tolist()
    df_region_stat_bed = df_stat_all_region_filter.filter(select_columns).query('gene_name in @list_gene_filter')
    ## 6.4.1 merge interval
    list_merge = list()
    for gene_name, sub_df_pass in df_region_stat_bed.groupby('gene_name'):
        print(gene_name)
        sub_df_pass = df_region_stat_bed.query('gene_name==@gene_name').dropna()
        tmp_transcript = df_stat_all_transcript_filter.query('gene_name==@gene_name').iloc[0].transcript_id
        sub_tmp = df_transcript_exon_1.query('transcript_id==@tmp_transcript')
        sub_tmp['num_pass'] = sub_tmp.apply(lambda x:str_list_num(sub_df_pass.region.tolist(), x.region),axis = 1)
        sub_df_fail = sub_tmp.query('num_pass==0').filter(['seqname','start','end','gene_id','gene_name','region'])
        if len(sub_df_fail) >0:
            sub_df_fail = f_recal_pos_exon_fail(sub_df_fail,tmp_transcript,  df_gene_fdr,df_transcript_exon_1)
        sub_df_fail.columns = select_columns
        sub_df_fail['filter'] = 'fail'
        sub_df_pass['filter'] = 'pass'
        sub_df = pd.concat([sub_df_fail, sub_df_pass])
        interval_sub_df = [[sub_df.iloc[i].pos_start,sub_df.iloc[i].pos_end,sub_df.iloc[i]['filter']] for i in range(len(sub_df))]
        merge_interval_list = merge_interval(interval_sub_df)
        sub_df_new = sub_df.iloc[0:len(merge_interval_list)].copy()
        sub_df_new.pos_start = [i[0] for i in merge_interval_list]
        sub_df_new.pos_end = [i[1] for i in merge_interval_list]
        sub_df_new['filter'] = [i[2] for i in merge_interval_list]
        sub_df_new['region'] = sub_df_new.apply(lambda x:x['gene_name'] + "|" + x['seqname'] + ':' + str(x['pos_start']) + "-" + str(x['pos_end']), axis = 1)
        list_merge.append(sub_df_new)
        del(sub_df_new)


    f_if_0cdna(list_merge)
    df_region_stat_bed_merge = pd.concat(list_merge)
    df_region_stat_bed_merge.columns = select_columns_rename
    f_generate_bed(df_region_stat_bed_merge)
    del(df_region_stat_bed_merge['region'])
    df_region_stat_bed_merge.to_csv(global_para.out_bed_merge, sep = "\t", index = False)


    # 7. tell users information
    global_para.logger.info("Program finished successfully")
    global_para.logger.info('%d cDNA detected'%(len(df_stat_all_transcript_filter_source_filter.gene_name.unique())))



if __name__ == "__main__":
    pass
