import pandas as pd
import numpy as np
import pysam
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Seq
import subprocess


try:
    from . import global_para
except ImportError:
    import global_para
try:
    from .unaligned_reads import *
except ImportError:
    from unaligned_reads import *

np.seterr(divide = 'ignore') 
pd.options.mode.chained_assignment = None  # default='warn', this is used to remove warning from pandas chain assignment



def bam_ispair():
    # define bam is paired or not;
	bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
	for read in bam_genome:
		return(read.is_paired)
		break

def exon_start_adj_get_unalign(df_exon):
    unalign_readlist_true = list();
    bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
    for read in bam_genome.fetch(df_exon.seqname,int(df_exon.pos_end)-1,int(df_exon.pos_end)+1):  
        if read.mapping_quality >=global_para.min_quality  and read.is_secondary == False and read.is_unmapped == False:
            tmp_cigar = read.to_dict()['cigar']
            tmp_nm = read.get_tag('NM')
            if tmp_cigar == '[0-9]*M' and tmp_nm == 0:
                pass;
            if  re.search("[0-9]*S$",tmp_cigar):
                if read.aend == df_exon.pos_end:
                    unalign_readlist_true.append(read);
                    continue
    return(unalign_readlist_true)



def exon_end_adj_get_unalign(df_exon):
    unalign_readlist_true = list();
    bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
    for read in bam_genome.fetch(df_exon.seqname,int(df_exon.pos_start)-1,int(df_exon.pos_start)+1): 
        if read.mapping_quality >=global_para.min_quality  and read.is_secondary == False and read.is_unmapped == False:
            tmp_cigar = read.to_dict()['cigar']
            tmp_nm = read.get_tag('NM')
            if tmp_cigar == '[0-9]*M' and tmp_nm == 0:
                pass;
            if  re.search("^[0-9]*S",tmp_cigar):
                if read.pos == (df_exon.pos_start -1 ):
                    unalign_readlist_true.append(read);
                    continue
    return(unalign_readlist_true)


def exon_adj_get_mate(read_list):
    bamFile_handle = pysam.Samfile(global_para.genome_bam_file,'rb')
    unalign_readlist_true_mate = list()
    read_list_new = [read for read in read_list if read.is_proper_pair]
    # for read in read_list:
    for read in read_list_new: # add
        try:
            read_mate = bamFile_handle.mate(read)
            unalign_readlist_true_mate.append(read_mate)
        except:
            pass
    return(unalign_readlist_true_mate)



def exon_adj_start_exchange(df_exon, df_exon_adj):
    df_exon.pos_start = df_exon_adj.pos_start
    df_exon.count_unalign_start = df_exon_adj.count_unalign_start
    df_exon.count_exon_match_start = df_exon_adj.count_exon_match_start
    df_exon.is_exon_match_start = df_exon_adj.is_exon_match_start
    df_exon.consensus_seq_start = df_exon_adj.consensus_seq_start
    df_exon.num_bad_reads_start = df_exon_adj.num_bad_reads_start
    df_exon.num_bad_reads_ajust_start = df_exon_adj.num_bad_reads_ajust_start
    df_exon.num_all_reads_start = df_exon_adj.num_all_reads_start
    df_exon.num_all_reads_ajust_start = df_exon_adj.num_all_reads_ajust_start
    df_exon.region_start = df_exon_adj.region_start
    df_exon.bg_unalign_start = df_exon_adj.bg_unalign_start
    df_exon.bg_total_start = df_exon_adj.bg_total_start
    df_exon.bbinom_p_start = df_exon_adj.bbinom_p_start
    df_exon.q_value_start = df_exon_adj.q_value_start
    return(df_exon)


def exon_adj_end_exchange(df_exon, df_exon_adj):
    df_exon.pos_end = df_exon_adj.pos_end
    df_exon.count_unalign_end = df_exon_adj.count_unalign_end
    df_exon.count_exon_match_end = df_exon_adj.count_exon_match_end
    df_exon.is_exon_match_end = df_exon_adj.is_exon_match_end
    df_exon.consensus_seq_end = df_exon_adj.consensus_seq_end
    df_exon.num_bad_reads_end = df_exon_adj.num_bad_reads_end
    df_exon.num_bad_reads_ajust_end = df_exon_adj.num_bad_reads_ajust_end
    df_exon.num_all_reads_end = df_exon_adj.num_all_reads_end
    df_exon.num_all_reads_ajust_end = df_exon_adj.num_all_reads_ajust_end
    df_exon.region_end = df_exon_adj.region_end
    df_exon.bg_unalign_end = df_exon_adj.bg_unalign_end
    df_exon.bg_total_end = df_exon_adj.bg_total_end
    df_exon.bbinom_p_end = df_exon_adj.bbinom_p_end
    df_exon.q_value_end = df_exon_adj.q_value_end
    return(df_exon)




def exon_start_adj(df_exon):
    if global_para.genome_bam_paired == True:
        df_exon_return = df_exon.copy()
        # get clipped reads in exon end position;
        unalign_readlist_true = exon_start_adj_get_unalign(df_exon)
        # get clipped mates read list;
        unalign_readlist_mate = exon_adj_get_mate(unalign_readlist_true)
        # get another boundary information
        tmp_read_df = f_readlist_unalign_boundary_todf(unalign_readlist_mate,'start', df_exon.start)
        tmp_read_df = tmp_read_df.query('pos<=@df_exon.start')
        if len(tmp_read_df) == 0:
            unalign_readlist = unalign_readlist_true +  unalign_readlist_mate
            unalign_readlist_pos = [read.pos for read in unalign_readlist]
            unalign_readlist_pos = [x for x in unalign_readlist_pos if x <df_exon.start]
            unalign_readlist_pos.sort()
            if len(unalign_readlist_pos)>=2:
                df_exon_return.pos_start = unalign_readlist_pos[1]
                df_exon_return = f_recount_unalign(df_exon_return,"start")
            df_exon_return.abs_d_start = abs(df_exon_return.pos_start - df_exon_return.start)
        else:
          tmp_pos = tmp_read_df.pos.value_counts().sort_values(ascending = False).index[0]
          df_exon_sub = df_exon.copy()
          df_exon_sub.start = tmp_pos
          df_exon_sub_t = df_exon_sub.to_frame().T
          df_exon_sub_t_adj = f_df_exon_start_stat(df_exon_sub_t)
          if len(df_exon_sub_t_adj) >=1:
              df_exon_sub_t_adj = f_df_gene_start_stat_remove_dup(df_exon_sub_t_adj)
              df_exon_sub_adj = df_exon_sub_t_adj.iloc[0]
              if df_exon_sub_adj.bbinom_p_start < df_exon.bbinom_p_start and df_exon_sub_adj.bbinom_p_start <0.01:
                  df_exon_return = exon_adj_start_exchange(df_exon,df_exon_sub_adj)
              else:
                  unalign_readlist = unalign_readlist_true +  unalign_readlist_mate
                  unalign_readlist_pos = [read.pos for read in unalign_readlist]
                  unalign_readlist_pos = [x for x in unalign_readlist_pos if x <df_exon.start]
                  unalign_readlist_pos.sort()
                  if len(unalign_readlist_pos)>=2:
                    df_exon_return.pos_start = unalign_readlist_pos[1]
                    df_exon_return = f_recount_unalign(df_exon_return,"start")
                  df_exon_return.abs_d_start = abs(df_exon_return.pos_start - df_exon_return.start)
    else:
        df_exon_return = df_exon.copy()
        # get clipped reads in exon end position;
        unalign_readlist_true = exon_start_adj_get_unalign(df_exon)
        # get another boundary information
        tmp_read_df = f_readlist_unalign_boundary_todf(unalign_readlist_true,'start', df_exon.start)
        tmp_read_df = tmp_read_df.query('pos<=@df_exon.start')
        if len(tmp_read_df) == 0:
            unalign_readlist = unalign_readlist_true
            unalign_readlist_pos = [read.pos for read in unalign_readlist]
            unalign_readlist_pos = [x for x in unalign_readlist_pos if x <df_exon.start]
            unalign_readlist_pos.sort()
            if len(unalign_readlist_pos)>=2:
                df_exon_return.pos_start = unalign_readlist_pos[1]
                df_exon_return = f_recount_unalign(df_exon_return,"start")
            df_exon_return.abs_d_start = abs(df_exon_return.pos_start - df_exon_return.start)
        else:
          tmp_pos = tmp_read_df.pos.value_counts().sort_values(ascending = False).index[0]
          df_exon_sub = df_exon.copy()
          df_exon_sub.start = tmp_pos
          df_exon_sub_t = df_exon_sub.to_frame().T
          df_exon_sub_t_adj = f_df_exon_start_stat(df_exon_sub_t)
          if len(df_exon_sub_t_adj)>=1:
              df_exon_sub_t_adj = f_df_gene_start_stat_remove_dup(df_exon_sub_t_adj)
              df_exon_sub_adj = df_exon_sub_t_adj.iloc[0]
              if df_exon_sub_adj.bbinom_p_start < df_exon.bbinom_p_start and df_exon_sub_adj.bbinom_p_start <0.01:
                  df_exon_return = exon_adj_start_exchange(df_exon,df_exon_sub_adj)
              else:
                  unalign_readlist = unalign_readlist_true
                  unalign_readlist_pos = [read.pos for read in unalign_readlist]
                  unalign_readlist_pos = [x for x in unalign_readlist_pos if x <df_exon.start]
                  unalign_readlist_pos.sort()
                  if len(unalign_readlist_pos)>=2:
                    df_exon_return.pos_start = unalign_readlist_pos[1]
                    df_exon_return = f_recount_unalign(df_exon_return,"start")
                  df_exon_return.abs_d_start = abs(df_exon_return.pos_start - df_exon_return.start)
    return(df_exon_return)


def exon_end_adj(df_exon):
    if global_para.genome_bam_paired == True:
        df_exon_return = df_exon.copy()
        # get clipped reads in exon end position;
        unalign_readlist_true = exon_end_adj_get_unalign(df_exon)
        # get clipped mates read list;
        unalign_readlist_mate = exon_adj_get_mate(unalign_readlist_true)
        # get another boundary information
        tmp_read_df = f_readlist_unalign_boundary_todf(unalign_readlist_mate,'end', df_exon.end)
        tmp_read_df = tmp_read_df.query('pos>=@df_exon.end')
        if len(tmp_read_df)>=1:
            tmp_pos = tmp_read_df.pos.value_counts().sort_values(ascending = False).index[0]
            df_exon_sub = df_exon.copy()
            df_exon_sub.end = tmp_pos
            df_exon_sub_t = df_exon_sub.to_frame().T
            df_exon_sub_t_adj = f_df_exon_end_stat(df_exon_sub_t)
            if len(df_exon_sub_t_adj)>=1:
                df_exon_sub_t_adj = f_df_gene_end_stat_remove_dup(df_exon_sub_t_adj)
                df_exon_sub_adj = df_exon_sub_t_adj.iloc[0]
                if df_exon_sub_adj.bbinom_p_end < df_exon.bbinom_p_end and df_exon_sub_adj.bbinom_p_end <0.01:
                    df_exon_return = exon_adj_end_exchange(df_exon,df_exon_sub_adj)
                else:
                    unalign_readlist = unalign_readlist_true +  unalign_readlist_mate
                    unalign_readlist_pos = [read.pos for read in unalign_readlist]
                    unalign_readlist_pos = [x for x in unalign_readlist_pos if x >df_exon.end]
                    unalign_readlist_pos.sort()
                    if len(unalign_readlist_pos)>=2:
                        df_exon_return.pos_end = unalign_readlist_pos[-2]
                        df_exon_return = f_recount_unalign(df_exon_return,"end")
                        df_exon_return.abs_d_end = abs(df_exon_return.pos_end - df_exon_return.end)
        else:
            unalign_readlist = unalign_readlist_true +  unalign_readlist_mate
            unalign_readlist_pos = [read.aend for read in unalign_readlist]
            unalign_readlist_pos = [x for x in unalign_readlist_pos if x >df_exon.end]
            unalign_readlist_pos.sort()
            if len(unalign_readlist_pos)>=2:
                df_exon_return.pos_end = unalign_readlist_pos[-2]
                df_exon_return = f_recount_unalign(df_exon_return,"end")
            df_exon_return.abs_d_end = abs(df_exon_return.pos_end - df_exon_return.end)
    else:
        df_exon_return = df_exon.copy()
        # get clipped reads in exon end position;
        unalign_readlist_true = exon_end_adj_get_unalign(df_exon)
        # get another boundary information
        tmp_read_df = f_readlist_unalign_boundary_todf(unalign_readlist_true,'end', df_exon.end)
        tmp_read_df = tmp_read_df.query('pos>=@df_exon.end')
        if len(tmp_read_df)>=1:
            tmp_pos = tmp_read_df.pos.value_counts().sort_values(ascending = False).index[0]
            df_exon_sub = df_exon.copy()
            df_exon_sub.end = tmp_pos
            df_exon_sub_t = df_exon_sub.to_frame().T
            df_exon_sub_t_adj = f_df_exon_end_stat(df_exon_sub_t)
            if len(df_exon_sub_t_adj)>=1:
                df_exon_sub_t_adj = f_df_gene_end_stat_remove_dup(df_exon_sub_t_adj)
                df_exon_sub_adj = df_exon_sub_t_adj.iloc[0]
                if df_exon_sub_adj.bbinom_p_end < df_exon.bbinom_p_end and df_exon_sub_adj.bbinom_p_end <0.01:
                    df_exon_return = exon_adj_end_exchange(df_exon,df_exon_sub_adj)
                else:
                    unalign_readlist = unalign_readlist_true
                    unalign_readlist_pos = [read.pos for read in unalign_readlist]
                    unalign_readlist_pos = [x for x in unalign_readlist_pos if x >df_exon.end]
                    unalign_readlist_pos.sort()
                    if len(unalign_readlist_pos)>=2:
                        df_exon_return.pos_end = unalign_readlist_pos[-2]
                        df_exon_return = f_recount_unalign(df_exon_return,"end")
                        df_exon_return.abs_d_end = abs(df_exon_return.pos_end - df_exon_return.end)
        else:
            unalign_readlist = unalign_readlist_true
            unalign_readlist_pos = [read.aend for read in unalign_readlist]
            unalign_readlist_pos = [x for x in unalign_readlist_pos if x >df_exon.end]
            unalign_readlist_pos.sort()
            if len(unalign_readlist_pos)>=2:
                df_exon_return.pos_end = unalign_readlist_pos[-2]
                df_exon_return = f_recount_unalign(df_exon_return,"end")
            df_exon_return.abs_d_end = abs(df_exon_return.pos_end - df_exon_return.end)
    return(df_exon_return)



def f_len_list(list_input):
    len_list = list()
    for x in list_input:
        if type(x) == str:
            len_list.append(max([len(x_single) for x_single in x.split(',')]))
        else:
            len_list.append(0)
    return(len_list)





def f_source_inference(df_stat_all_region, table_blastn):
    column_names = ["gene_name", "clipped_seq_left", "clipped_seq_right", "clipped_left_CDS_distance", "clipped_right_CDS_distance", "clipped_left_length", "clipped_right_length", "clipped_left_blastn", "clipped_right_blastn", "source_inference","first_choice"]
    genelist = df_stat_all_region.gene_name.unique().tolist()
    df_genelist = pd.DataFrame({'gene_name':genelist})
    table_blastn_source = pd.DataFrame(columns = column_names,index = genelist);
    table_blastn_source.gene_name = genelist
    tmp_df_start = df_stat_all_region.query('is_exon_boundary_start=="1"')
    tmp_df_end = df_stat_all_region.query('is_exon_boundary_end=="1"')
    table_blastn_source.clipped_seq_left = df_genelist.merge(tmp_df_start.filter(['gene_name','clipped_seq_start']).groupby('gene_name')['clipped_seq_start'].apply(','.join).reset_index(), how = "left").clipped_seq_start.tolist()
    table_blastn_source.clipped_seq_right = df_genelist.merge(tmp_df_end.filter(['gene_name','clipped_seq_end']).groupby('gene_name')['clipped_seq_end'].apply(','.join).reset_index(), how = "left").clipped_seq_end.tolist()
    table_blastn_source.clipped_left_CDS_distance = df_genelist.merge(tmp_df_start.filter(['gene_name','abs_d_start']).groupby('gene_name')['abs_d_start'].apply(max).reset_index(), how = "left").abs_d_start.tolist()
    table_blastn_source.clipped_right_CDS_distance = df_genelist.merge(tmp_df_end.filter(['gene_name','abs_d_end']).groupby('gene_name')['abs_d_end'].apply(max).reset_index(), how = "left").abs_d_end.tolist()
    table_blastn_source.clipped_left_length = f_len_list(table_blastn_source.clipped_seq_left.tolist())
    table_blastn_source.clipped_right_length = f_len_list(table_blastn_source.clipped_seq_right.tolist())
    if len(table_blastn)==0:
        pass;
    else:
        table_blastn = table_blastn.sort_values(['query_acc.ver','bit_score'],ascending = False).drop_duplicates(subset = 'query_acc.ver')
        table_blastn = table_blastn.query('pct_identity>=95')
        if len(table_blastn)==0:
            pass
        else:
            table_blastn['gene_name'] = table_blastn.apply(lambda x:x["query_acc.ver"].split(":")[0],axis = 1)
            table_blastn['which_boundary'] = table_blastn.apply(lambda x:x["query_acc.ver"].split(":")[1],axis = 1)
            table_blastn['blastn_database'] = table_blastn.apply(lambda x:x["subject_acc.ver"].split(":")[0],axis = 1)
            table_blastn['blastn_database_detail'] = table_blastn.apply(lambda x:x["subject_acc.ver"].split(":")[1],axis = 1)
            table_blastn = pd.concat([table_blastn.query('blastn_database!="UTR"'),table_blastn.query('blastn_database=="UTR"').query('blastn_database_detail==gene_name')],sort = False)
            table_blastn_source.clipped_left_blastn = df_genelist.merge(table_blastn.query('which_boundary=="start"').filter(['gene_name','blastn_database']),how = "left").blastn_database.tolist()
            table_blastn_source.clipped_right_blastn = df_genelist.merge(table_blastn.query('which_boundary=="end"').filter(['gene_name','blastn_database']),how = "left").blastn_database.tolist()
            table_blastn_source.first_choice = df_genelist.merge(table_blastn.sort_values(['gene_name',"bit_score"],ascending = False).drop_duplicates(subset = "gene_name"), how = "left").blastn_database.tolist()
    table_blastn_source.source_inference = table_blastn_source.apply(lambda x:f_source_inference_single(x),axis = 1)
    del(table_blastn_source['first_choice'])
    return(table_blastn_source)




def f_source_inference_single(df_blastn_gene):
    tmp_inference = "unknown"
    tmp_blastn_database_lable = ['vector', 'retro_element', 'UTR']
    tmp_blastn_database_result = ['vector', 'retrocopy', 'retrocopy']
    tmp_blastn_database_dict = {tmp_blastn_database_lable[i]:tmp_blastn_database_result[i] for i in range(len(tmp_blastn_database_lable))}
    if np.any([pd.isna(df_blastn_gene.clipped_left_CDS_distance),pd.isna(df_blastn_gene.clipped_right_CDS_distance)]):
        if df_blastn_gene.clipped_left_CDS_distance<=5 or df_blastn_gene.clipped_right_CDS_distance <=5:
            if df_blastn_gene.first_choice == "retro_element":
               tmp_inference = "unknown"
            elif df_blastn_gene.first_choice in tmp_blastn_database_lable:
                tmp_inference = tmp_blastn_database_dict[df_blastn_gene.first_choice]
            elif pd.isna(df_blastn_gene.first_choice) == False and df_blastn_gene.first_choice not in tmp_blastn_database_lable:
                tmp_inference = df_blastn_gene.first_choice
            else:
                tmp_inference = "unknown"
        elif df_blastn_gene.clipped_left_CDS_distance>5 or df_blastn_gene.clipped_right_CDS_distance >5:
            if df_blastn_gene.first_choice in tmp_blastn_database_lable:
                tmp_inference = tmp_blastn_database_dict[df_blastn_gene.first_choice]
            elif pd.isna(df_blastn_gene.first_choice) == False and df_blastn_gene.first_choice not in tmp_blastn_database_lable:
                tmp_inference = df_blastn_gene.first_choice
            elif df_blastn_gene.clipped_left_CDS_distance>10 or df_blastn_gene.clipped_right_CDS_distance >10:
                tmp_inference = "retrocopy-likely"
            else:
                tmp_inference = "unknown"
        else:
            tmp_inference = "unknown"
    elif np.all([pd.isna(df_blastn_gene.clipped_left_CDS_distance)==False,pd.isna(df_blastn_gene.clipped_right_CDS_distance)==False]):
        if df_blastn_gene.clipped_left_CDS_distance<=5 and df_blastn_gene.clipped_right_CDS_distance <=5:
            if df_blastn_gene.first_choice == "retro_element":
               tmp_inference = "unknown"
            elif df_blastn_gene.first_choice in tmp_blastn_database_lable:
                tmp_inference = tmp_blastn_database_dict[df_blastn_gene.first_choice]
            elif pd.isna(df_blastn_gene.first_choice) == False and df_blastn_gene.first_choice not in tmp_blastn_database_lable:
                tmp_inference = df_blastn_gene.first_choice
            else:
                tmp_inference = "vector-likely"
        elif df_blastn_gene.clipped_left_CDS_distance>10 and df_blastn_gene.clipped_right_CDS_distance >10:
            if df_blastn_gene.first_choice in tmp_blastn_database_lable:
                tmp_inference = tmp_blastn_database_dict[df_blastn_gene.first_choice]
            elif (pd.isna(df_blastn_gene.first_choice) == False) and (df_blastn_gene.first_choice not in tmp_blastn_database_lable):
                tmp_inference = df_blastn_gene.first_choice
            else:
                tmp_inference = "retrocopy-likely"
        else:
            if df_blastn_gene.first_choice in tmp_blastn_database_lable:
                tmp_inference = tmp_blastn_database_dict[df_blastn_gene.first_choice]
            elif pd.isna(df_blastn_gene.first_choice) == False and df_blastn_gene.first_choice not in tmp_blastn_database_lable:
                tmp_inference = df_blastn_gene.first_choice
            else:
                tmp_inference = "unknown"
    else:
            tmp_inference = "unknown"
    return(tmp_inference)


def f_recount_unalign(df_exon_return,type):
  if type=="start":
    bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
    unalign_readlist = f_find_unalign_readlist_boundary(bam_genome,df_exon_return.seqname, df_exon_return.pos_start-1,df_exon_return.pos_start+1,"start")
    unalign_readlist_filter = [read for read in unalign_readlist if read.pos==(df_exon_return.pos_start-1)]
    df_exon_return.num_bad_reads_start = len(unalign_readlist_filter)
    df_exon_return.num_bad_reads_ajust_start = df_exon_return.num_bad_reads_start
    num_all_reads = len(f_find_unalign_readlist(bam_genome,df_exon_return.seqname, df_exon_return.pos_start-1,df_exon_return.pos_start)[0])
    df_exon_return.num_all_reads_start = num_all_reads
    df_exon_return.num_all_reads_ajust_start = num_all_reads
    df_exon_return.bbinom_p_start = beta_binomial_significance(df_exon_return.num_bad_reads_ajust_start, df_exon_return.num_all_reads_ajust_start, df_exon_return.bg_unalign_start, df_exon_return.bg_total_start)
  if type == "end":
    bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
    unalign_readlist = f_find_unalign_readlist_boundary(bam_genome,df_exon_return.seqname, df_exon_return.pos_end-1,df_exon_return.pos_end+1,"end")
    unalign_readlist_filter = [read for read in unalign_readlist if read.aend==(df_exon_return.pos_end)]
    df_exon_return.num_bad_reads_end = len(unalign_readlist_filter)
    df_exon_return.num_bad_reads_ajust_end = df_exon_return.num_bad_reads_end
    num_all_reads = len(f_find_unalign_readlist(bam_genome,df_exon_return.seqname, df_exon_return.pos_end-1,df_exon_return.pos_end)[0])
    df_exon_return.num_all_reads_end = num_all_reads
    df_exon_return.num_all_reads_ajust_end = num_all_reads
    df_exon_return.bbinom_p_end = beta_binomial_significance(df_exon_return.num_bad_reads_ajust_end, df_exon_return.num_all_reads_ajust_end, df_exon_return.bg_unalign_end, df_exon_return.bg_total_end)
  return(df_exon_return)




def f_source_known(gene_name, table_source_known):
    tmp_df = table_source_known.query('gene_name==@gene_name')
    if len(tmp_df)>0:
        source_known = ",".join(tmp_df['source'].tolist())
    else:
        source_known = "None"
    return(source_known)


def f_df_source_known(file_source_known):
    if file_source_known == "":
        df = pd.DataFrame(columns = ['gene_name',"source"])
    else:
        df = pd.read_table(file_source_known, header = 0, sep = '\t')
    return df





def f_get_blast_4gene(df_stat_all_region_full_reinfo):
    # identify source;
    records = list()
    for gene_name, df_gene in df_stat_all_region_full_reinfo.groupby('gene_name'):
        tmp_df_seq_start = df_gene.query('is_exon_boundary_start=="1"').query('q_value_start<@global_para.cutoff_pvalue').filter(['gene_name','clipped_seq_start'])
        if len(tmp_df_seq_start)>=1:
            tmp_df_seq_start['name'] = tmp_df_seq_start.apply(lambda x:x.gene_name + ":start",axis = 1)
            tmp_df_seq_start.columns = ['gene_name','clipped_seq','name']
        tmp_df_seq_end = df_gene.query('is_exon_boundary_end=="1"').query('q_value_end<@global_para.cutoff_pvalue').filter(['gene_name','clipped_seq_end'])
        if len(tmp_df_seq_end)>=1:
            tmp_df_seq_end['name'] = tmp_df_seq_end.apply(lambda x:x.gene_name + ":end",axis = 1)
            tmp_df_seq_end.columns = ['gene_name','clipped_seq','name']
        tmp_df_seq = pd.concat([tmp_df_seq_start, tmp_df_seq_end], sort = False)
        if len(tmp_df_seq)>=1:
            for i in range(len(tmp_df_seq)):
              if len(tmp_df_seq.iloc[i]['clipped_seq']) >=6:
                record = SeqRecord(Seq.Seq(tmp_df_seq.iloc[i]['clipped_seq']),id=tmp_df_seq.iloc[i]['name'])
                records.append(record)
        del(tmp_df_seq_start)
        del(tmp_df_seq_end)
        del(tmp_df_seq)
    list_blastn_columns = ["query_acc.ver", "subject_acc.ver", "pct_identity", "alignment_length", "mismatches", "gap_opens", "q._start", "q._end", "s._start", "s._end", "evalue", "bit_score"]
    if len(records)>=1:
        SeqIO.write(records, global_para.out_blastn_seq, "fasta")
        blastn_command = "bash " + os.path.join(global_para.script_path, "scripts/detect/f_blastn_source_inference.sh") + " " +  global_para.out_blastn_seq + " " + global_para.blastn_database
        global_para.logger.info("Running blastn for inferring sources")
        subprocess.call(blastn_command,shell = True)
        if os.path.isfile(global_para.out_blastn_seq_table) and os.path.getsize(global_para.out_blastn_seq_table):
            table_blastn = pd.read_csv(global_para.out_blastn_seq_table, header = 0, sep = '\t')
            table_blastn.columns = list_blastn_columns
        else:
            table_blastn = pd.DataFrame(columns = list_blastn_columns)
    else:
        table_blastn = pd.DataFrame(columns = list_blastn_columns)
    return table_blastn




def f_filter_by_source(df_stat_all_transcript_filter):
    if global_para.source_inference_include!="" and global_para.source_inference_exclude == "":
        if global_para.source_inference_include == "all":
            tmp_list_boolean_inference = [True for x in df_stat_all_transcript_filter.source_inference.tolist()]
        else:
            list_source_inference_include = global_para.source_inference_include.split(',')
            tmp_list_boolean_inference = [True if str_list_compare(list_source_inference_include, str(x)) else False for x in df_stat_all_transcript_filter.source_inference.tolist()]
    elif global_para.source_inference_include=="" and global_para.source_inference_exclude != "":
        list_source_inference_exclude = global_para.source_inference_exclude.split(',')
        tmp_list_boolean_inference = [False if str_list_compare(list_source_inference_exclude, str(x)) else True for x in df_stat_all_transcript_filter.source_inference.tolist()]
    else:
        global_para.logger.error('error parameter input')
        exit(1)
    if global_para.source_known_databases_include!="" and global_para.source_known_databases_exclude=="":
        if global_para.source_known_databases_include == "all":
            tmp_list_boolean_known = [True for x in df_stat_all_transcript_filter.source_known_databases]
        else:
            list_source_known_include = global_para.source_known_databases_include.split(',')
            tmp_list_boolean_known = [True if str_list_compare(list_source_known_include, str(x)) else False for x in df_stat_all_transcript_filter.source_known_databases.tolist()]
    elif global_para.source_known_databases_include == "" and global_para.source_known_databases_exclude != "":
        list_source_known_exclude = global_para.source_known_databases_exclude.split(',')
        tmp_list_boolean_known = [False if str_list_compare(list_source_known_exclude, str(x)) else True for x in df_stat_all_transcript_filter.source_known_databases.tolist()]
    else:
        global_para.logger.error('error parameter input')
        exit(1)
    df_stat_all_transcript_filter_source_filter = df_stat_all_transcript_filter[np.array([(a and b) for a,b in zip(tmp_list_boolean_inference,tmp_list_boolean_known)])]
    return df_stat_all_transcript_filter_source_filter





def f_recal_pos_exon_fail(sub_df_fail, tmp_transcript, df_gene_fdr, df_transcript_exon_1):
    for region in sub_df_fail.region:
        # print(region)
        tmp_df_gene_fdr_region = df_gene_fdr.query('region==@region').dropna()
        index_value = sub_df_fail.query('region==@region').index
        if (len(tmp_df_gene_fdr_region)>0 and tmp_df_gene_fdr_region.iloc[0]['is_exon_match_start']):
            sub_df_fail.loc[index_value,'start'] = int(tmp_df_gene_fdr_region.iloc[0]['pos_start']) -1 
        else:
            tmp_df_transcript = df_transcript_exon_1.query('transcript_id==@tmp_transcript').query('region==@region')
            if tmp_df_transcript.is_exon_boundary_start.iloc[0] == "0":
                n_distance = f_sequence_pair_compare_n(tmp_df_transcript.exon_flank_start20.iloc[0], tmp_df_transcript.exon_boundary_start_nearseq20.iloc[0],'start')
                tmp_start_new = int(tmp_df_transcript.start.iloc[0]) -1 - n_distance
            else:
                bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
                tmp_start = tmp_df_transcript.iloc[0].start -1 
                all_readlist,unalign_readlist,n_unalign_readlist_boundary, n_all_readlist_boundary = f_find_unalign_readlist(bam_genome,tmp_df_transcript.iloc[0].seqname,tmp_df_transcript.iloc[0].start-10, tmp_df_transcript.iloc[0].start)
                tmp_start_new = int(tmp_df_transcript.start.iloc[0]) -1
                if len(unalign_readlist)>0:
                    tmp_read_df = f_readlist_unalign_boundary_todf(unalign_readlist,'start', tmp_df_transcript.iloc[0].start)
                    tmp_read_df_filter = tmp_read_df.query('pos>=(@tmp_start-10) and pos<=@tmp_start')
                    if len(tmp_read_df_filter)>=1:
                        tmp_start_new = pd.DataFrame(tmp_read_df_filter.pos.value_counts()).reset_index().sort_values(by  = ['pos','index'], ascending = [0,1]).iloc[0]['index'] - 1
            sub_df_fail.loc[index_value,'start'] = tmp_start_new
        if len(tmp_df_gene_fdr_region)>0 and tmp_df_gene_fdr_region.iloc[0]['is_exon_match_end']:
            sub_df_fail.loc[index_value,'end'] = int(tmp_df_gene_fdr_region.iloc[0]['pos_end'])
        else:
            tmp_df_transcript = df_transcript_exon_1.query('transcript_id==@tmp_transcript').query('region==@region')
            if tmp_df_transcript.is_exon_boundary_end.iloc[0] == "0":
                n_distance = f_sequence_pair_compare_n(tmp_df_transcript.exon_flank_end20.iloc[0], tmp_df_transcript.exon_boundary_end_nearseq20.iloc[0],'end')
                tmp_end_new = int(tmp_df_transcript.end.iloc[0]) +  n_distance
            else:
                bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
                tmp_end = tmp_df_transcript.iloc[0].end
                all_readlist,unalign_readlist,n_unalign_readlist_boundary, n_all_readlist_boundary = f_find_unalign_readlist(bam_genome,tmp_df_transcript.iloc[0].seqname,tmp_end , tmp_end + 10)
                tmp_end_new = int(tmp_df_transcript.end.iloc[0])
                if len(unalign_readlist)>0:
                    tmp_read_df = f_readlist_unalign_boundary_todf(unalign_readlist,'end', tmp_df_transcript.iloc[0].end)
                    tmp_read_df_filter = tmp_read_df.query('pos>=(@tmp_end) and pos<=@tmp_end + 10')
                    if len(tmp_read_df_filter)>=1:
                        tmp_end_new = pd.DataFrame(tmp_read_df_filter.pos.value_counts()).reset_index().sort_values(by  = ['pos','index'], ascending = [0,1]).iloc[0]['index']
            sub_df_fail.loc[index_value,'end'] = tmp_end_new
    return sub_df_fail


