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
from statsmodels.stats.multitest import fdrcorrection
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

def f_0cdna():
    # if 0 cdna detected, report messages to users
    global_para.logger.info("Program finished successfully")
    global_para.logger.info("No cDNA detected. Exiting.")
    exit(0)

def f_if_0cdna(obj):
    if len(obj) == 0:
        f_0cdna()



def read_gene_model(gtf_gene_unique_file):
    # load gene model into a dataframe
    print('Loading gene model table')
    dict_type = {
    "seqname":"str",
    "start":"int64",
    "end":"int64",
    "gene_id":"str",
    "gene_name":"str",
    "transcript_id":"str",
    "exon_flank_start20":"str",
    "exon_flank_end20":"str",
    "is_exon_boundary_start":"str",
    "is_exon_boundary_end":"str",
    "exon_boundary_start_nearseq20":"str",
    "exon_boundary_end_nearseq20":"str"}
    df_gene_exon_unique = pd.read_csv(gtf_gene_unique_file, sep = '\t',header = 0)
    df_gene_exon_unique = df_gene_exon_unique.astype(dict_type)
    # convert all sequences to uppercase
    df_gene_exon_unique['exon_flank_start20'] = df_gene_exon_unique['exon_flank_start20'].str.upper()
    df_gene_exon_unique['exon_flank_end20'] = df_gene_exon_unique['exon_flank_end20'].str.upper()
    df_gene_exon_unique['exon_boundary_start_nearseq20'] = df_gene_exon_unique['exon_boundary_start_nearseq20'].str.upper()
    df_gene_exon_unique['exon_boundary_end_nearseq20'] = df_gene_exon_unique['exon_boundary_end_nearseq20'].str.upper()
    df_gene_exon_unique = df_gene_exon_unique.fillna('')
    print('Loaded %d exons\n'%(len(df_gene_exon_unique)))
    return(df_gene_exon_unique)

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

def f_overlap_reference(genome_bam_file,df_gene_exon):
    # overlap reference for input bam and gene model
    bam_genome = pysam.AlignmentFile(genome_bam_file,'rb')
    reference_bam = bam_genome.references
    bam_genome.close()
    reference_exon = df_gene_exon.seqname.unique().tolist()
    overlap_reference = [x for x in reference_bam if x in reference_exon]
    if len(overlap_reference)==0: global_para.logger.error('chromosome names are not matched between gene model and bam file'); exit(1)
    df_gene_exon = df_gene_exon.query('seqname in @overlap_reference')
    return df_gene_exon
    

def f_close_exon_merge(df_transcript_exon):
    df_transcript_exon = df_transcript_exon.sort_values(['transcript_id','start'])
    df_transcript_exon = df_transcript_exon.reset_index(drop = True)
    df_transcript_exon['start_next'] = df_transcript_exon.groupby(['transcript_id'])['start'].shift(-1)
    df_transcript_exon['dis_exon'] = abs(df_transcript_exon['end'] - df_transcript_exon['start_next'])
    df_transcript_exon_close = df_transcript_exon.query('dis_exon<@global_para.exon_distance')
    list_transcript = df_transcript_exon_close.transcript_id.unique().tolist()
    if len(list_transcript) >0:
        list_df_transcript_merge = []
        for transcript_id in list_transcript:
            sub_df = df_transcript_exon.query('transcript_id==@transcript_id')
            sub_df_new = f_df_1transcript_merge(sub_df)
            list_df_transcript_merge.append(sub_df_new)
        df_transcript_exon_close_new = pd.concat(list_df_transcript_merge)
        df_transcript_exon_noclose = df_transcript_exon.query('transcript_id not in @list_transcript')
        df_transcript_exon_new = pd.concat([df_transcript_exon_close_new, df_transcript_exon_noclose])
    else:
        df_transcript_exon_new = df_transcript_exon.copy()
    del df_transcript_exon_new['start_next']
    del df_transcript_exon_new['dis_exon']
    return df_transcript_exon_new


def f_df_1transcript_merge(sub_df):
    list_new = []
    list_iskeep = [True]*len(sub_df)
    for i in range(len(sub_df)):
        # print("line %s"%i)
        if list_iskeep[i] and sub_df.iloc[i].dis_exon <global_para.exon_distance:
            j = 1
            up_list = sub_df.iloc[[i]]
            while up_list.iloc[0].dis_exon < global_para.exon_distance:
                down_list = sub_df.iloc[[i+j]]
                up_list['end'] = down_list.iloc[0]['end']
                up_list['exon_flank_end20'] = down_list.iloc[0]['exon_flank_end20']
                up_list['is_exon_boundary_end'] = down_list.iloc[0]['is_exon_boundary_end']
                up_list['exon_boundary_end_nearseq20'] = down_list.iloc[0]['exon_boundary_end_nearseq20']
                up_list['start_next'] = down_list.iloc[0]['start_next']
                up_list['dis_exon'] = down_list.iloc[0]['dis_exon']
                list_iskeep[i+j] = False
                j = j+1
            list_new.append(up_list)
        elif list_iskeep[i] == False:
            continue
        else:
            list_new.append(sub_df.iloc[[i]])
    sub_df_new = pd.concat(list_new)
    return sub_df_new


def f_exon_group_by_distance(df_gene_exon_unique):
    # split gene model dataframe with distance <=150 into one new columns : read_group
    list_df = list()
    df_gene_exon_unique = df_gene_exon_unique.sort_values(['seqname','start','end'])
    for name, sub_df in df_gene_exon_unique.groupby('seqname'):
        sub_df['diff_start'] = (sub_df.start - sub_df.end.shift(1)).tolist()
        sub_df['diff_end'] = (sub_df.start.shift(-1) - sub_df.end).tolist()
        list_df.append(sub_df)
    df_gene_exon_unique_merge = pd.concat(list_df)
    df_gene_exon_unique_merge['if_merge'] = True
    df_gene_exon_unique_merge.loc[df_gene_exon_unique_merge.query('diff_start>150 and diff_end>150').index,'if_merge'] = False
    list_group = list()
    i = -1
    tmp_value_previous = ''
    for value in df_gene_exon_unique_merge.if_merge:
        if tmp_value_previous == value:
            list_group.append(i)
        else:
            i = i+1;
            list_group.append(i)
            tmp_value_previous = value
    df_gene_exon_unique_merge['read_group_1'] = list_group
    df_gene_exon_unique_merge['read_group'] = df_gene_exon_unique_merge.read_group_1//500
    df_gene_exon_unique_merge_new = pd.concat([sub_df for name, sub_df in df_gene_exon_unique_merge.groupby('read_group')])
    del df_gene_exon_unique_merge_new['read_group_1']
    del df_gene_exon_unique_merge_new['diff_start']
    del df_gene_exon_unique_merge_new['diff_end']
    return(df_gene_exon_unique_merge_new)



def f_df_add_region(df):
    df_s_region =  df.gene_name + "|" + df.seqname + ":" + df.start.astype('str') + "-" + df.end.astype('str')
    return df_s_region


def f_get_ratio_unaligned(global_para):
    if global_para.num_total_read_exon>0:
        p_unalign_read = global_para.num_unalign_read/global_para.num_total_read_exon
    else:
        f_0cdna()

    if p_unalign_read>=0.2: 
        global_para.logger.warning("Warning: Proportion of clipped reads in exon regions is high! Pay attention to detected cDNAs (if detected).")
    global_para.logger.info('Proportion of soft-clipped reads is: %0.3f\n'%(p_unalign_read))
    return p_unalign_read





def f_filter_genelist_step1(df_transcript_exon, df_gene_exon):
    df_gene_exon_expand = df_transcript_exon.copy()
    df_gene_exon['pvalue'] = df_gene_exon.pvalue.fillna(1)
    df_gene_exon_expand = df_gene_exon_expand.merge(df_gene_exon, on = 'region', suffixes = ['','y'])
    df_gene_exon_expand['pvalue_transcript'] = df_gene_exon_expand.groupby(['transcript_id'])['pvalue'].transform(f_combin_p)
    df_gene_exon_expand['fdr_transcript'] = fdrcorrection(df_gene_exon_expand['pvalue_transcript'])[1]
    tmp_genelist_filter1 = df_gene_exon_expand.query('pvalue_transcript<=@global_para.cutoff_pvalue').gene_name.unique().tolist()
    if len(tmp_genelist_filter1) <= 100:
        tmp_genelist_filter = tmp_genelist_filter1
    else:
        tmp_genelist_filter1 = df_gene_exon_expand.query('fdr_transcript<=@global_para.cutoff_pvalue').gene_name.unique().tolist()
        if len(tmp_genelist_filter1)<=500:
            tmp_genelist_filter = tmp_genelist_filter1
        else:
            df_gene_exon_expand_filter = df_gene_exon_expand.query('gene_name in @tmp_genelist_filter1') 
            df_gene_exon_expand_filter['n_transcript'] = df_gene_exon_expand_filter.groupby(['transcript_id'])['transcript_id'].transform('count')
            df_gene_exon_expand_filter['n_transcript_detect'] = df_gene_exon_expand_filter.query('pvalue<=@global_para.cutoff_pvalue').groupby(['transcript_id'])['transcript_id'].transform('count')
            df_gene_exon_expand_filter['n_transcript_detect'] = df_gene_exon_expand_filter.n_transcript_detect.fillna(0)
            df_gene_exon_expand_filter['ratio'] = df_gene_exon_expand_filter.apply(lambda x:x.n_transcript_detect/x.n_transcript,axis = 1)
            df_gene_exon_expand_filter = df_gene_exon_expand_filter.query('ratio>=@global_para.cutoff_ratio_gene')
            tmp_genelist_filter = df_gene_exon_expand_filter.gene_name.unique().tolist()
            if len(df_gene_exon_expand_filter.gene_name.unique())>=500:
                tmp_genelist_filter = df_gene_exon_expand_filter.sort_values('fdr_transcript')['gene_name'].unique()[0:500].tolist()
    f_if_0cdna(tmp_genelist_filter)
    return tmp_genelist_filter

def f_df_collapse_seq_boundary_info(df_gene_exon_unique_transcript):
    list_region = []
    list_is_exon_boundary_start = []
    list_is_exon_boundary_end = []
    list_exon_boundary_start_nearseq20 = []
    list_exon_boundary_end_nearseq20 = []
    list_transcript = []
    for region, sub_df in df_gene_exon_unique_transcript.groupby('region'):
        list_region.append(region)
        list_is_exon_boundary_start.append(','.join(set(sub_df.is_exon_boundary_start.tolist())))
        list_is_exon_boundary_end.append(','.join(set(sub_df.is_exon_boundary_end.tolist())))
        list_exon_boundary_start_nearseq20.append(','.join(set(sub_df.exon_boundary_start_nearseq20.tolist())))
        list_exon_boundary_end_nearseq20.append(','.join(set(sub_df.exon_boundary_end_nearseq20.tolist())))
        list_transcript.append(','.join(set(sub_df.transcript_id.tolist())))
    dict_region = {
    'region': list_region, 
    'is_exon_boundary_start': list_is_exon_boundary_start, 
    'is_exon_boundary_end': list_is_exon_boundary_end, 
    'exon_boundary_start_nearseq20': list_exon_boundary_start_nearseq20, 
    'exon_boundary_end_nearseq20': list_exon_boundary_end_nearseq20, 
    'transcript_id': list_transcript
    }
    df_region = pd.DataFrame(dict_region)
    return(df_region)





def f_adjust_bg_unaligned_reads(df_gene_exon_unique):
    tmp_unalign_count_each_gene = df_gene_exon_unique.groupby('gene_name').num_unalign_read_exon_edge.sum()
    tmp_unalign_count_each_gene_filter = tmp_unalign_count_each_gene[tmp_unalign_count_each_gene >= max(global_para.ratio_ehc*global_para.num_total_read_exon, global_para.count_ehc)]
    if len(tmp_unalign_count_each_gene_filter) >=1:
        tmp_ehc_gene = ",".join(tmp_unalign_count_each_gene_filter.index.tolist())
        global_para.logger.info("genes with extremely high coverage of clipped reads are detected, corresponding reads will not be account for background calculation")
        global_para.logger.info("%s"%tmp_ehc_gene)
        global_para.num_unalign_read = global_para.num_unalign_read - sum(tmp_unalign_count_each_gene_filter.values.tolist())



def f_expand_by_transcript(df_gene_exon_unique):
    # expand transcript by seperation ','
    new_df = pd.DataFrame(df_gene_exon_unique.transcript_id.str.split(',').tolist(),index = df_gene_exon_unique.index).stack()
    new_df.index = new_df.index.droplevel(-1)
    new_df.name = 'transcript_split'
    df_gene_exon_unique_expand = df_gene_exon_unique.join(new_df)
    n_transcript = df_gene_exon_unique_expand.groupby(['transcript_split'])['transcript_split'].count()
    n_transcript.name = 'n_transcript'
    sum_transcript = df_gene_exon_unique_expand.groupby(['transcript_split'])['num_unalign_read_exon_edge'].sum()
    sum_transcript.name = 'read_transcript'
    df_transcript = pd.concat([n_transcript,sum_transcript],axis = 1)
    df_gene_exon_unique_expand = df_gene_exon_unique_expand.merge(df_transcript,left_on = 'transcript_split',right_index = True)
    df_gene_exon_unique_expand['cutoff_transcript_read'] = global_para.cutoff_num_exon_unaligned_reads*df_gene_exon_unique_expand.n_transcript
    return(df_gene_exon_unique_expand)



## get unaligned sequences from one read
def f_get_unalign_seq(read,pos,type):
    if type == 'start':
        tmp_seq_ref = read.seq[0:(pos -1 - read.pos )]
        tmp_seq_align = tmp_seq_ref[::-1]
    else:
        tmp_seq = read.seq[::-1][0:(read.aend - pos )]
        tmp_seq_ref = tmp_seq[::-1]
        tmp_seq_align = tmp_seq_ref
    return([tmp_seq_ref, tmp_seq_align])


# get read id from one read
def f_readid(readlist):
    readid_list = list();
    for read in readlist:
        tmp_readid = read.to_dict()['name'] + "_" +  str(int(read.is_read1))
        readid_list.append(tmp_readid)
    return(readid_list)


# count clipped reads number for multiple exon regions
# This is used in the background exon scanning
def f_find_unalign_readlist_multi(df_sub):
    bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
    df_sub_2list = df_sub.apply(lambda x:f_find_unalign_readlist(bam_genome,x.seqname,x.start-1, x.end),axis = 1)
    df_sub['list'] = df_sub_2list
    list2_1 = df_sub.query('if_merge==False').list.tolist()
    list1_1 = df_sub.query('if_merge==True').list.tolist()
    tmp_num_total_read_exon = sum([len([read for sub_read in list2_1 for read in sub_read[0]]),len(set([read for sub_read in list1_1 for read in sub_read[0]]))])
    tmp_num_unalign_read = sum([len([read for sub_read in list2_1 for read in sub_read[1]]),len(set([read for sub_read in list1_1 for read in sub_read[1]]))])
    tmp_num_unalign_read_exon = [len(i[1]) for i in df_sub_2list]
    tmp_num_unalign_read_exon_edge = [i[2] for i in df_sub_2list]
    tmp_num_all_read_exon_edge = [i[3] for i in df_sub_2list]
    bam_genome.close()
    return([tmp_num_total_read_exon,tmp_num_unalign_read,tmp_num_unalign_read_exon, tmp_num_unalign_read_exon_edge,tmp_num_all_read_exon_edge,df_sub['gene_name'].tolist()])


# count clipped reads for one exon regions
# just used for background exon scanning
def f_find_unalign_readlist(bam_genome,chr, start, end):
    unalign_readlist = list();
    all_readlist = list();
    unalign_readlist_boundary = list()
    n_unalign_readlist_boundary = 0
    n_all_readlist_boundary = 0
    for read in bam_genome.fetch(chr,start,end):       
        if read.mapping_quality >=global_para.min_quality and read.is_secondary == False and read.is_unmapped == False:
            all_readlist.append(read)
    if len(all_readlist)>0:
        all_readlist_cigar = [read for read in all_readlist if (read.cigar[0][0]==4 or read.cigar[-1][0]==4)]
        if read.has_tag('MD'):
            all_readlist_md = [read for read in all_readlist if read.has_tag('MD') and re.search("^0[ACGT]|[ACGT]0$",read.get_tag('MD'))]
            unalign_readlist = list(set(all_readlist_md+all_readlist_cigar))    
        else:
            unalign_readlist = all_readlist_cigar
        n_unalign_readlist_boundary = len([read for read in unalign_readlist if read.pos <= (start +5 ) or read.aend >=(end-5)])
        n_all_readlist_boundary = len([read for read in all_readlist if read.pos <= (start +5 ) or read.aend >=(end-5)])
    return(all_readlist,unalign_readlist,n_unalign_readlist_boundary, n_all_readlist_boundary)


# convert clipped reads into data.frame
def f_read_unalign_type_todf(readlist,type):
    readlist_readid = list();
    readlist_type = list();
    readlist_seq_ref = list();
    readlist_seq_align = list();
    readlist_seqname  = list();
    readlist_pos  = list();
    for read in readlist:
        tmp_readid = read.to_dict()['name'] + "_" +  str(int(read.is_read1))
        tmp_seqname = read.to_dict()['ref_name']
        tmp_read_pos_start = int(read.to_dict()['ref_pos'])
        tmp_read_pos_end =  int(read.aend)
        if type == 'start':
            tmp_cigar = read.cigar[0][0]
            if tmp_cigar == 4:
                tmp_type = 'start'
                tmp_seq_n = int(read.cigar[0][1])
                tmp_seq_ref = read.seq[0:tmp_seq_n]
                tmp_seq_align = tmp_seq_ref[::-1]
                readlist_readid.append(tmp_readid)
                readlist_seqname.append(tmp_seqname)
                readlist_pos.append(tmp_read_pos_start)
                readlist_type.append(tmp_type)
                readlist_seq_ref.append(tmp_seq_ref)
                readlist_seq_align.append(tmp_seq_align)
            if read.has_tag('MD') and re.search("^0[ACGT]",read.get_tag('MD')):
                tmp_type = 'start'
                tmp_seq_n = 1
                tmp_seq_ref = read.seq[0:tmp_seq_n]
                tmp_seq_align = tmp_seq_ref[::-1]
                readlist_readid.append(tmp_readid)
                readlist_seqname.append(tmp_seqname)
                readlist_pos.append(tmp_read_pos_start + 1)
                readlist_type.append(tmp_type)
                readlist_seq_ref.append(tmp_seq_ref)
                readlist_seq_align.append(tmp_seq_align)
        if type == 'end':
            tmp_cigar = read.cigar[-1][0]
            if tmp_cigar == 4:
                tmp_type = 'end'
                tmp_seq_n = int(read.cigar[-1][1])
                tmp_seq = read.seq[::-1][0:tmp_seq_n]
                tmp_seq_ref = tmp_seq[::-1]
                tmp_seq_align = tmp_seq_ref
                readlist_readid.append(tmp_readid)
                readlist_seqname.append(tmp_seqname)
                readlist_pos.append(tmp_read_pos_end)
                readlist_type.append(tmp_type)
                readlist_seq_ref.append(tmp_seq_ref)
                readlist_seq_align.append(tmp_seq_align)
            if read.has_tag('MD') and re.search("[ACGT]0$",read.get_tag('MD')):
                tmp_type = 'end'
                tmp_seq_n = 1
                tmp_seq = read.seq[::-1][0:tmp_seq_n]
                tmp_seq_ref = tmp_seq[::-1]
                tmp_seq_align = tmp_seq_ref
                readlist_readid.append(tmp_readid)
                readlist_seqname.append(tmp_seqname)
                readlist_pos.append(tmp_read_pos_end - 1)
                readlist_type.append(tmp_type)
                readlist_seq_ref.append(tmp_seq_ref)
                readlist_seq_align.append(tmp_seq_align)
            
    df = pd.DataFrame({'readid':readlist_readid,'seqname':readlist_seqname,'pos':readlist_pos,'type':readlist_type,'breakpoint_seq_ref':readlist_seq_ref,'breakpoint_seq_align': readlist_seq_align})
    return(df)

    

def f_cal_bad_read_all_reads_specific_pos(bam_genome, df_exon, pos,type, unalign_readlist_df):
    # correct clipped reads number with a specific position
    chr = df_exon.iloc[0]['seqname']
    unalign_readlist_df_filter = unalign_readlist_df.query('pos==@pos')
    consensus_seq = f_get_consensus_sequence(unalign_readlist_df_filter['breakpoint_seq_align'].tolist())
    num_bad_reads = len(unalign_readlist_df_filter)
    all_reads = f_find_unalign_readlist(bam_genome,chr,pos-1,pos)[0]
    num_all_reads = len(all_reads)
    all_reads_readid = f_readid(all_reads)
    tmp_all_reads_isin_unalign = list(~pd.Series(all_reads_readid).isin(unalign_readlist_df_filter['readid'].tolist()))
    all_reads_diff = [all_reads[i] for i in range(0,len(tmp_all_reads_isin_unalign)) if tmp_all_reads_isin_unalign[i]]
    all_reads_diff_align_seq = [f_get_unalign_seq(read,pos,type)[1] for read in all_reads_diff]
    num_all_reads_ajust = num_all_reads - len([seq for seq in all_reads_diff_align_seq if len(seq)==0])
    all_reads_diff_align_seq = [seq for seq in all_reads_diff_align_seq if len(seq)>0]
    num_all_reads_consensus = f_consensus_seq_count_exact(consensus_seq, all_reads_diff_align_seq)[consensus_seq]
    num_bad_reads_ajust = num_bad_reads + num_all_reads_consensus
    return([consensus_seq, num_bad_reads, num_bad_reads_ajust,num_all_reads,num_all_reads_ajust])



def f_find_unalign_readlist_boundary(bam_genome,chr,start,end,type):
    unalign_readlist = list();
    if type == 'start':
        for read in bam_genome.fetch(chr,start,end):
            if read.mapping_quality >=global_para.min_quality  and read.is_secondary == False and read.is_unmapped == False:
                tmp_cigar = read.cigar[0][0]
                tmp_nm = read.get_tag('NM')
                if tmp_cigar == 0 and tmp_nm == 0:
                    pass;
                elif tmp_cigar == 4:
                    unalign_readlist.append(read)
                    continue;
                elif read.has_tag('MD') and re.search("^0[ACGT]",read.get_tag('MD')):
                    unalign_readlist.append(read)
                    continue;
    elif type == 'end':
        for read in bam_genome.fetch(chr,start,end):  
            if read.mapping_quality >=global_para.min_quality  and read.is_secondary == False and read.is_unmapped == False:
                tmp_cigar = read.cigar[-1][0]
                tmp_nm = read.get_tag('NM')
                if tmp_cigar == 0 and tmp_nm == 0:
                    pass;
                elif  tmp_cigar == 4:
                        unalign_readlist.append(read);
                        continue;
                elif read.has_tag('MD') and re.search("[ACGT]0$",read.get_tag('MD')):
                        unalign_readlist.append(read)
    else:
        print('error input') 
    return(unalign_readlist)





def f_readlist_unalign_boundary_todf(readlist,boundary_type, exon_pos):
    readlist_readid = list();
    readlist_boundary_type = list();
    readlist_seq_ref = list();
    readlist_seq_align = list();
    readlist_seqname  = list();
    readlist_pos  = list();
    readlist_seq_ref_full =  list();
    for read in readlist:
        tmp_readid = read.to_dict()['name'] + "_" +  str(int(read.is_read1))
        tmp_seqname = read.to_dict()['ref_name']
        tmp_read_pos_start = int(read.to_dict()['ref_pos'])
        tmp_read_pos_end =  int(read.aend)
        if boundary_type == 'start':
            tmp_cigar = read.cigar[0][0]
            if tmp_cigar == 4:
                tmp_boundary_type = 'start'
                tmp_seq_n = int(read.cigar[0][1])
                tmp_seq_ref = read.seq[0:tmp_seq_n]
                tmp_seq_align = tmp_seq_ref[::-1]
                if exon_pos-1 >= read.pos:
                    tmp_seq_ref_full  = read.seq[0:(tmp_seq_n + read.get_overlap(read.pos ,exon_pos-1))]
                else:
                    tmp_seq_ref_full = ''
                readlist_readid.append(tmp_readid)
                readlist_seqname.append(tmp_seqname)
                readlist_pos.append(tmp_read_pos_start)
                readlist_boundary_type.append(tmp_boundary_type)
                readlist_seq_ref.append(tmp_seq_ref)
                readlist_seq_align.append(tmp_seq_align)
                readlist_seq_ref_full.append(tmp_seq_ref_full)
            elif read.has_tag('MD') and re.search("^0[ACGT]",read.get_tag('MD')):
                tmp_boundary_type = 'start'
                tmp_seq_n = 1
                tmp_seq_ref = read.seq[0:tmp_seq_n]
                tmp_seq_align = tmp_seq_ref[::-1]
                tmp_seq_ref_full = read.seq[0:read.get_overlap(read.pos, exon_pos -1)]
                readlist_readid.append(tmp_readid)
                readlist_seqname.append(tmp_seqname)
                readlist_pos.append(tmp_read_pos_start + 1)
                readlist_boundary_type.append(tmp_boundary_type)
                readlist_seq_ref.append(tmp_seq_ref)
                readlist_seq_align.append(tmp_seq_align)
                readlist_seq_ref_full.append(tmp_seq_ref_full)
        if boundary_type == 'end':
            tmp_cigar = read.cigar[-1][0]
            if tmp_cigar == 4:
                tmp_boundary_type = 'end'
                tmp_seq_n = int(read.cigar[-1][1])
                tmp_seq = read.seq[::-1][0:tmp_seq_n]
                tmp_seq_ref = tmp_seq[::-1]
                tmp_seq_align = tmp_seq_ref
                if read.aend >= exon_pos:
                    tmp_seq_ref_full = read.seq[::-1][0:(tmp_seq_n + read.get_overlap(exon_pos, read.aend) )][::-1]
                else:
                    tmp_seq_ref_full = ''
                readlist_readid.append(tmp_readid)
                readlist_seqname.append(tmp_seqname)
                readlist_pos.append(tmp_read_pos_end)
                readlist_boundary_type.append(tmp_boundary_type)
                readlist_seq_ref.append(tmp_seq_ref)
                readlist_seq_ref_full.append(tmp_seq_ref_full)
                readlist_seq_align.append(tmp_seq_align)
            elif read.has_tag('MD') and re.search("[ACGT]0$",read.get_tag('MD')):
                tmp_boundary_type = 'end'
                tmp_seq_n = 1
                tmp_seq = read.seq[::-1][0:tmp_seq_n]
                tmp_seq_ref = tmp_seq[::-1]
                tmp_seq_align = tmp_seq_ref
                tmp_seq_ref_full = read.seq[::-1][0:(read.get_overlap(exon_pos, read.aend) )][::-1]
                readlist_readid.append(tmp_readid)
                readlist_seqname.append(tmp_seqname)
                readlist_pos.append(tmp_read_pos_end - 1)
                readlist_boundary_type.append(tmp_boundary_type)
                readlist_seq_ref.append(tmp_seq_ref)
                readlist_seq_align.append(tmp_seq_align)            
                readlist_seq_ref_full.append(tmp_seq_ref_full)
    df = pd.DataFrame({'readid':readlist_readid,'seqname':readlist_seqname,'pos':readlist_pos,'boundary_type':readlist_boundary_type,'breakpoint_seq_ref':readlist_seq_ref,'breakpoint_seq_align': readlist_seq_align,'breakpoint_seq_ref_full':readlist_seq_ref_full})
    return(df)





def merge_interval(list_interval):
    # interval_sub_df = [[sub_df.iloc[i].pos_start,sub_df.iloc[i].pos_end,sub_df.iloc[i]['filter']] for i in range(len(sub_df))]
    # list_interval = interval_sub_df.copy()
    list_interval.sort(key = lambda interval: interval[0])
    list_merge_interval = list()
    list_merge_interval.append(list_interval[0])
    for interval in list_interval:
        if interval[0] <= list_merge_interval[-1][1]:
            list_merge_interval[-1][1] = max(list_merge_interval[-1][1], interval[1])
            list_merge_interval[-1][2] = "fail" if np.any([i == "fail" for i in [list_merge_interval[-1][2], interval[2]]]) else "pass"
        else:
            list_merge_interval.append(interval)
    return(list_merge_interval)





def f_startswith_str_2(str1, str2):
    if len(str1)>=len(str2):
        is_starts = str1.startswith(str2)
    else:
        is_starts = str2.startswith(str1)
    return is_starts

def f_endswith_str_2(str1, str2):
    if len(str1)>=len(str2):
        is_starts = str1.endswith(str2)
    else:
        is_starts = str2.endswith(str1)
    return is_starts



def f_df_exon_start_stat(df_exon):
    # i = 1
    # df_exon = df_gene.iloc[[i]]
    # 1. get bad readlist here: start:
    def f_get_exon_match_start(x):
        return(any([(f_endswith_str_2(exon_last_element, x['breakpoint_seq_ref_full']) and  x['breakpoint_seq_ref_full']!="") for exon_last_element in exon_last_seqlist]) )
    tmp_start = df_exon.iloc[0]['start']
    tmp_end = df_exon.iloc[0]['end']
    tmp_seqname = df_exon.iloc[0]['seqname']
    region = df_exon.iloc[0]['region']
    bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
    unaligned_readlist_start = f_find_unalign_readlist_boundary(bam_genome,tmp_seqname,max(0,tmp_start-5),tmp_start + 5,'start')
    if len(unaligned_readlist_start)==0:
        tmp_pos_df_filter_stat = pd.DataFrame()
    else:
        tmp_read_df = f_readlist_unalign_boundary_todf(unaligned_readlist_start,'start', tmp_start)
        exon_last_seqlist = [x.upper() for x in df_exon.iloc[0]['exon_boundary_start_nearseq20'].split(',')]
        tmp_read_df['exon_seq_match_start'] = False
        if exon_last_seqlist ==['']:
            pass
        else:
            tmp_read_df['exon_seq_match_start']   = tmp_read_df.apply(lambda x:f_get_exon_match_start(x),axis = 1)
        # 1. get exon position
        dict_pos = dict()
        for pos, tmp_read_df_sub in tmp_read_df.groupby('pos'):
            dict_sub_pos = dict()
            dict_sub_pos['count_unalign'] = len(tmp_read_df_sub)
            dict_sub_pos['count_exon_match'] = sum(tmp_read_df_sub['exon_seq_match_start'])
            dict_sub_pos['is_exon_match'] = dict_sub_pos['count_exon_match']>0
            dict_pos[pos] = dict_sub_pos
        # filter of exons
        tmp_pos_df = pd.DataFrame.from_dict(dict_pos,orient = 'index').reset_index()
        tmp_pos_df.columns = ['pos_start','count_unalign','count_exon_match','is_exon_match']
        if any(tmp_pos_df['is_exon_match']):
            tmp_pos_df_filter = tmp_pos_df.query('count_exon_match>0')
            tmp_list = tmp_pos_df['count_unalign']
            tmp_pos_df_filter = tmp_pos_df[tmp_list>max(tmp_list)/2]
        else:
            tmp_list = tmp_pos_df['count_unalign']
            tmp_pos_df_filter = tmp_pos_df[tmp_list>max(tmp_list)/2]
        # 3. make statistics
        y = tmp_pos_df_filter.apply(lambda x:f_cal_bad_read_all_reads_specific_pos(bam_genome, df_exon,x.pos_start ,'start', tmp_read_df),axis = 1)
        tmp_pos_df_filter_stat = tmp_pos_df_filter.merge(pd.DataFrame(y.tolist(),columns = ['consensus_seq','num_bad_reads','num_bad_reads_ajust','num_all_reads','num_all_reads_ajust'],index = y.index),left_index = True, right_index = True)
        tmp_pos_df_filter_stat['region'] = region
        tmp_pos_df_filter_stat['exon_start'] = tmp_start
        bam_genome.close()
    return(tmp_pos_df_filter_stat)



def f_df_exon_end_stat(df_exon):
    def f_get_exon_match_end(x):
        return(any([(f_startswith_str_2(exon_last_element, x['breakpoint_seq_ref_full']) and  x['breakpoint_seq_ref_full']!="") for exon_last_element in exon_last_seqlist]) )
    tmp_start = df_exon.iloc[0]['start']
    tmp_end = df_exon.iloc[0]['end']
    tmp_seqname = df_exon.iloc[0]['seqname']
    tmp_type = 'end'
    tmp_exon_pos = tmp_end
    region = df_exon.iloc[0]['region']
    bam_genome = pysam.AlignmentFile(global_para.genome_bam_file,'rb')
    unaligned_readlist = f_find_unalign_readlist_boundary(bam_genome,tmp_seqname,max(tmp_exon_pos-5,0),tmp_exon_pos + 5,tmp_type)
    if len(unaligned_readlist)==0:
        tmp_pos_df_filter_stat = pd.DataFrame()
    else:
        tmp_read_df = f_readlist_unalign_boundary_todf(unaligned_readlist,tmp_type, tmp_exon_pos)
        exon_last_seqlist = [x.upper() for x in df_exon.iloc[0]['exon_boundary_end_nearseq20'].split(',')]
        tmp_read_df['exon_seq_match_end'] = False
        if exon_last_seqlist ==['']:
            pass
        else:
            tmp_read_df['exon_seq_match_end']   = tmp_read_df.apply(lambda x:f_get_exon_match_end(x),axis = 1)
        dict_pos = dict()
        for pos, tmp_read_df_sub in tmp_read_df.groupby('pos'):
            dict_sub_pos = dict()
            dict_sub_pos['count_unalign'] = len(tmp_read_df_sub)
            dict_sub_pos['count_exon_match'] = sum(tmp_read_df_sub['exon_seq_match_end'])
            dict_sub_pos['is_exon_match'] = dict_sub_pos['count_exon_match']>0
            dict_pos[pos] = dict_sub_pos
        tmp_pos_df = pd.DataFrame.from_dict(dict_pos,orient = 'index').reset_index()
        tmp_pos_df.columns = ['pos_end','count_unalign','count_exon_match','is_exon_match']
        if any(tmp_pos_df['is_exon_match']):
            tmp_pos_df_filter = tmp_pos_df.query('count_exon_match>0')
            tmp_list = tmp_pos_df_filter['count_unalign']
            tmp_pos_df_filter = tmp_pos_df_filter[tmp_list>max(tmp_list)/2]
        else:
            tmp_list = tmp_pos_df['count_unalign']
            tmp_pos_df_filter = tmp_pos_df[tmp_list>max(tmp_list)/2]
        # 3. statistics
        y = tmp_pos_df_filter.apply(lambda x:f_cal_bad_read_all_reads_specific_pos(bam_genome, df_exon,x.pos_end ,tmp_type, tmp_read_df),axis = 1)
        tmp_pos_df_filter_stat = tmp_pos_df_filter.merge(pd.DataFrame(y.tolist(),columns = ['consensus_seq','num_bad_reads','num_bad_reads_ajust','num_all_reads','num_all_reads_ajust'],index = y.index),left_index = True, right_index = True)
        tmp_pos_df_filter_stat['region'] = region
        tmp_pos_df_filter_stat['exon_end'] = tmp_end
        bam_genome.close()
    return(tmp_pos_df_filter_stat)




def f_df_gene_start_stat_remove_dup(df_stat_start):
    df_stat_start['abs_d'] = df_stat_start.apply(lambda x:abs(x.pos_start - x.exon_start ), axis = 1)
    df_stat_start['bg_unalign'] = global_para.num_unalign_read
    df_stat_start['bg_total'] = global_para.num_total_read_exon
    df_stat_start['bbinom_p'] = df_stat_start.apply(lambda x:beta_binomial_significance(x.num_bad_reads_ajust, x.num_all_reads_ajust, x.bg_unalign, x.bg_total), axis=1)
    df_stat_start['q_value'] = fdrcorrection(df_stat_start['bbinom_p'])[1]
    df_stat_start.columns = [x  if x.endswith('_start') else (x + '_start') for x in df_stat_start.columns]
    # remove duplication:
    # 1. one exon have multiple postions;
    df_stat_start = df_stat_start.reset_index()
    df_stat_start_tmpdup = df_stat_start[df_stat_start.duplicated(['region_start'],keep = False)]
    if len(df_stat_start_tmpdup)>=1:
        list_drop = list()
        for tmp_pos, sub_df in df_stat_start_tmpdup.groupby(['region_start']):
            sub_df_is_exon_true = sub_df[sub_df['is_exon_match_start']];
            sub_df_is_exon_false= sub_df[~sub_df['is_exon_match_start']]
            if len(sub_df_is_exon_true)>0 and len(sub_df_is_exon_false)>0:
                list_drop.append(sub_df_is_exon_false)
            elif len(sub_df_is_exon_true)==0 and len(sub_df_is_exon_false)>0:
                sub_df_is_exon_false_sort = sub_df_is_exon_false.sort_values(by = ['q_value_start','abs_d_start'])
                sub_df_is_exon_false_drop = sub_df_is_exon_false_sort[sub_df_is_exon_false_sort.duplicated('region_start',keep = "first")]
                list_drop.append(sub_df_is_exon_false_drop)
            else:
                pass
        if len(list_drop)>0:
            df_stat_start = df_stat_start.drop(pd.concat(list_drop).index)
    return(df_stat_start)

def f_df_gene_end_stat_remove_dup(df_stat_end):
    df_stat_end['abs_d'] = df_stat_end.apply(lambda x:abs(x.pos_end - x.exon_end), axis = 1)
    df_stat_end['bg_unalign'] = global_para.num_unalign_read
    df_stat_end['bg_total'] = global_para.num_total_read_exon
    df_stat_end['bbinom_p'] = df_stat_end.apply(lambda x:beta_binomial_significance(x.num_bad_reads_ajust, x.num_all_reads_ajust, x.bg_unalign, x.bg_total), axis=1)
    df_stat_end['q_value'] = fdrcorrection(df_stat_end['bbinom_p'])[1]
    df_stat_end.columns = [x  if x.endswith('_end') else (x + '_end') for x in df_stat_end.columns]
    # remove duplicate
    df_stat_end = df_stat_end.reset_index()
    df_stat_end_tmpdup = df_stat_end[df_stat_end.duplicated(['region_end'],keep = False)]
    if len(df_stat_end_tmpdup)>=1:
        list_drop = list()
        for tmp_pos, sub_df in df_stat_end_tmpdup.groupby(['region_end']):
            sub_df_is_exon_true = sub_df[sub_df['is_exon_match_end']];
            sub_df_is_exon_false = sub_df[~sub_df['is_exon_match_end']]
            if len(sub_df_is_exon_true)>0 and len(sub_df_is_exon_false)>0:
                list_drop.append(sub_df_is_exon_false)
            elif len(sub_df_is_exon_true)==0 and len(sub_df_is_exon_false)>0:
                # sub_df_is_exon_false_drop = sub_df_is_exon_false[sub_df_is_exon_false.sort_values(by = ['q_value_end','abs_d_end']).duplicated('region_end',keep = "first")]
                sub_df_is_exon_false_sort = sub_df_is_exon_false.sort_values(by = ['q_value_end','abs_d_end'])
                sub_df_is_exon_false_drop = sub_df_is_exon_false_sort[sub_df_is_exon_false_sort.duplicated('region_end',keep = "first")]
                list_drop.append(sub_df_is_exon_false_drop)
            else:
                pass
        if len(list_drop)>0:
            df_stat_end = df_stat_end.drop(pd.concat(list_drop).index)
    return(df_stat_end)




def f_df_gene_merge_removedup(df_gene):
    df_gene_full = df_gene.copy()
    # df_gene = df_gene_full.copy()
    df_gene = df_gene.dropna(subset = ['pos_start','pos_end'], how = "any", axis = 0)
    # 2. remove duplication:
    # 2. one position have multiple exons
    ## drop start
    df_gene_start_tmpdup = df_gene[df_gene.pos_start.duplicated(keep = False)]
    if len(df_gene_start_tmpdup)>0:
        list_drop = list()
        for tmp_pos, sub_df in df_gene_start_tmpdup.groupby(['pos_start']):
            sub_df_is_exon_true = sub_df[sub_df['is_exon_match_end']==True];
            if len(sub_df_is_exon_true) ==0 and len(sub_df.query('is_exon_boundary_start.str.contains("1")',engine = 'python'))>=1:
                sub_df_is_exon_true = sub_df.query('is_exon_boundary_start.str.contains("1")', engine = 'python')
            sub_df_is_exon_false = sub_df.drop(sub_df_is_exon_true.index)
            if len(sub_df_is_exon_true)>0 and len(sub_df_is_exon_false)>0:
                list_drop.append(sub_df_is_exon_false)
            elif len(sub_df_is_exon_true)==0 and len(sub_df_is_exon_false)>0:
                sub_df_is_exon_false_sort = sub_df_is_exon_false.sort_values(by = ['q_value_end','abs_d_end'])
                sub_df_is_exon_false_drop = sub_df_is_exon_false_sort[sub_df_is_exon_false_sort.duplicated('pos_start',keep = "first")]
                list_drop.append(sub_df_is_exon_false_drop)
            else:
                pass
        if len(list_drop)>0:
            df_gene = df_gene.drop(pd.concat(list_drop).index)
    ## drop end
    df_gene_end_tmpdup = df_gene[df_gene.pos_end.duplicated(keep = False)]
    if len(df_gene_end_tmpdup)>0:
        list_drop = list()
        for tmp_pos, sub_df in df_gene_end_tmpdup.groupby(['pos_end']):
            sub_df_is_exon_true = sub_df[sub_df['is_exon_match_start']==True];
            if len(sub_df_is_exon_true) ==0 and len(sub_df.query('is_exon_boundary_end.str.contains("1")', engine = 'python'))>=1:
                sub_df_is_exon_true = sub_df.query('is_exon_boundary_end.str.contains("1")', engine = 'python')
            sub_df_is_exon_false = sub_df.drop(sub_df_is_exon_true.index)   
            if len(sub_df_is_exon_true)>0 and len(sub_df_is_exon_false)>0:
                list_drop.append(sub_df_is_exon_false)
            elif len(sub_df_is_exon_true)==0 and len(sub_df_is_exon_false)>0:
                sub_df_is_exon_false_sort = sub_df_is_exon_false.sort_values(by = ['q_value_start','abs_d_start'])
                sub_df_is_exon_false_drop = sub_df_is_exon_false_sort[sub_df_is_exon_false_sort.duplicated('pos_end',keep = "first")]
                list_drop.append(sub_df_is_exon_false_drop)
            else:
                pass
        if len(list_drop)>0:
            df_gene = df_gene.drop(pd.concat(list_drop).index)
    # start na selection
    na_start = df_gene_full.bbinom_p_start.isna().tolist()
    na_end   = df_gene_full.bbinom_p_end.isna().tolist()
    df_gene_start_na_end = df_gene_full[[ (not a) and b for a,b in zip(na_start, na_end)]]
    df_gene_end_na_start = df_gene_full[[a and (not b) for a,b in zip(na_start, na_end)]]
    df_gene_start_na_end_sort = pd.DataFrame(columns = df_gene_start_na_end.columns)
    df_gene_end_na_start_sort = pd.DataFrame(columns = df_gene_start_na_end.columns)
    if len(df_gene_start_na_end) >=1:
        region_is_in_list = [a in df_gene.region.tolist() for a in df_gene_start_na_end.region.tolist()]
        pos_start_is_in_list = [a in df_gene.pos_start.tolist() for a in df_gene_start_na_end.pos_start.tolist()]
        df_gene_start_na_end = df_gene_start_na_end[[(not a) and (not b) for (a,b) in zip(region_is_in_list, pos_start_is_in_list)]]
        if df_gene_start_na_end.shape[0]>=1:
            df_gene_start_na_end = df_gene_start_na_end.query('is_exon_match_start == True')
            if df_gene_start_na_end.shape[0]>=1:
                df_gene_start_na_end_sort = df_gene_start_na_end.sort_values(by = ['q_value_start','abs_d_start'])
                df_gene_start_na_end_sort = df_gene_start_na_end_sort.drop_duplicates('region',keep = "first")
                df_gene_start_na_end_sort['q_value_end'] = 1
                df_gene_start_na_end_sort['bbinom_p_end'] = 1
                df_gene_start_na_end_sort['is_na'] = 1
    # end na selection
    if len(df_gene_end_na_start)>=1:
        region_is_in_list = [a in df_gene.region.tolist() for a in df_gene_end_na_start.region.tolist()]
        pos_end_is_in_list = [a in df_gene.pos_end.tolist() for a in df_gene_end_na_start.pos_end.tolist()]
        df_gene_end_na_start = df_gene_end_na_start[[(not a) and (not b) for (a,b) in zip(region_is_in_list, pos_end_is_in_list)]]
        if df_gene_end_na_start.shape[0]>=1:
            df_gene_end_na_start = df_gene_end_na_start.query('is_exon_match_end == True')
            if df_gene_end_na_start.shape[0]>=1:
                df_gene_end_na_start_sort = df_gene_end_na_start.sort_values(by = ['q_value_end','abs_d_end'])
                df_gene_end_na_start_sort = df_gene_end_na_start_sort.drop_duplicates('region',keep = "first")
                df_gene_end_na_start_sort['q_value_start'] = 1
                df_gene_end_na_start_sort['bbinom_p_start'] = 1
                df_gene_end_na_start_sort['is_na'] = 1
    # df_gene merge
    df_gene['is_na'] = 0
    df_gene_full = pd.concat([df_gene,df_gene_start_na_end_sort, df_gene_end_na_start_sort])
    return(df_gene_full)





def f_exon_filter(df_exon, global_para):
    if all([re.search("1", df_exon.is_exon_boundary_start), re.search("1", df_exon.is_exon_boundary_end)]):
        is_exon_match_region = True
        is_pvalue = all([df_exon.combine_qvalue <= global_para.cutoff_pvalue,df_exon.bbinom_p_start <= global_para.cutoff_pvalue, df_exon.bbinom_p_end <= global_para.cutoff_pvalue])
    else:
        is_pvalue = df_exon.combine_qvalue <= global_para.cutoff_pvalue
        is_exon_match_region = any([df_exon.is_exon_match_start and len(df_exon.clipped_seq_start) >=1 , df_exon.is_exon_match_end and len(df_exon.clipped_seq_end) >=1])
    is_keep = all([is_pvalue,is_exon_match_region])
    return(is_keep)






def f_transcript_stat(transcript_id,df_stat_all_region, df_gene_exon_unique):
    sub_df_stat_all_region = df_stat_all_region[df_stat_all_region.apply(lambda x: transcript_id in x.transcript_id.split(',') , axis = 1)]
    tmp_gene_id = sub_df_stat_all_region.gene_id.iloc[0]
    tmp_gene_name = sub_df_stat_all_region.gene_name.iloc[0]
    tmp_num_exon_detected = sub_df_stat_all_region.filter(['seqname','exon_start','exon_end']).drop_duplicates().shape[0]
    tmp_num_exon_transcript = df_gene_exon_unique[df_gene_exon_unique.apply(lambda x: transcript_id in x.transcript_id.split(',') , axis = 1)].query('gene_id==@tmp_gene_id').shape[0]
    return([tmp_gene_id, tmp_gene_name,transcript_id,tmp_num_exon_detected,tmp_num_exon_transcript])

def str_list_compare(list_str, str_input):
    # to define if one of list_str can be part of str_input;
    result = [x for x in list_str if x in str_input]
    if len(result)>0:
        return(True)
    else:
        return(False)




def str_list_num(list_str, str_input):
    result = len([x for x in list_str if str_input in x])
    return(result)


def df_max(df, col_name):
    max_value = df.sort_values(col_name, ascending = False)[col_name].iloc[0]
    df_max_value = df[df[col_name]==max_value]
    return(df_max_value)




def f_transcript_source(genename,df_stat_all_transcript_filter, df_gene_fdr, df_gene_exon_unique_filter_exon):
    tmp_df_transcript_filter = df_stat_all_transcript_filter.query("gene_name==@genename")
    tmp_df_select_max = df_max(tmp_df_transcript_filter, 'num_exon_detected')
    tmp_transcript = ''
    if len(tmp_df_select_max) == 1:
        tmp_transcript = tmp_df_select_max.iloc[0].transcript_id
    elif len(tmp_df_select_max) >1:
        tmp_df_gene_fdr = df_gene_fdr.query("gene_name==@genename").query('combine_qvalue<0.05')
        tmp_df_select_max['num_tmp_fdr'] = tmp_df_select_max.apply(lambda x:str_list_num(tmp_df_gene_fdr.transcript_id.tolist(), x.transcript_id),axis = 1)
        tmp_df_select_max = df_max(tmp_df_select_max,'num_tmp_fdr')
        tmp_transcript = tmp_df_select_max.iloc[0].transcript_id
        if len(tmp_df_select_max)>1:
            tmp_df_gene_exon_unique = df_gene_exon_unique_filter_exon.query("gene_name==@genename")
            tmp_df_select_max['num_tmp_uniq'] = tmp_df_select_max.apply(lambda x:str_list_num(tmp_df_gene_exon_unique.query('pvalue<0.05').transcript_id.tolist(), x.transcript_id),axis = 1)
            tmp_df_select_max = df_max(tmp_df_select_max,'num_tmp_uniq')
            tmp_transcript = tmp_df_select_max.iloc[0].transcript_id
            if len(tmp_df_select_max)>1:
                tmp_df_select_max = df_max(tmp_df_select_max,'ratio')
                tmp_transcript = tmp_df_select_max.iloc[0].transcript_id
                if len(tmp_df_select_max)>1:
                    tmp_df_select_max = tmp_df_select_max.iloc[[0]]
                    tmp_transcript = tmp_df_select_max.iloc[0].transcript_id
    return(tmp_transcript)



def f_generate_bed(df_region_stat_bed_merge):
    global_para.out_bed_merge_color = os.path.join(global_para.output_dir, (global_para.sampleid + '.merge_region.bed'))
    dict_color = dict({"pass":"0,153,153", "fail":"128,128,128"})
    df_bed = df_region_stat_bed_merge.iloc[:,0:3].copy()
    df_bed['name'] = df_region_stat_bed_merge.iloc[:,5]
    df_bed['score'] = '.'
    df_bed['strand'] = '.'
    df_bed['thickStart'] = df_bed['start']
    df_bed['thickEnd'] = df_bed['end']
    df_bed['itemRgb'] = df_region_stat_bed_merge.apply(lambda x:dict_color[x['filter']],axis = 1)
    f = open(global_para.out_bed_merge_color, 'w')
    tmp_comment_line = 'track name="'  + global_para.sampleid + '" description="cDNA regions" itemRgb="On"\n'
    f.write(tmp_comment_line)
    df_bed.to_csv(f, sep = "\t", index = False, header = False)
    f.close()




def f_stat_most_possible_transcript(df_stat_all_region, df_gene_exon_1, df_gene_fdr, df_gene_exon_filter):
    list_transcript = list(set([x for sublist in df_stat_all_region.transcript_id.tolist() for x in sublist.split(',')]))
    df_stat_all_transcript = pd.DataFrame({"transcript_id":list_transcript}).apply(lambda x:f_transcript_stat(x.transcript_id,df_stat_all_region,df_gene_exon_1),axis = 1,result_type = "expand")
    df_stat_all_transcript.columns = ['gene_id','gene_name','transcript_id','num_exon_detected','num_exon_transcript']
    df_stat_all_transcript['ratio'] = df_stat_all_transcript.apply(lambda x: round(x.num_exon_detected/x.num_exon_transcript, 3), axis=1)
    df_stat_possible_transcript = pd.DataFrame({'gene_name':df_stat_all_transcript.gene_name.unique().tolist(), 'possible_transcript':"unknown"})
    df_stat_possible_transcript['possible_transcript'] = df_stat_possible_transcript.apply(lambda x:f_transcript_source(x.gene_name, df_stat_all_transcript, df_gene_fdr, df_gene_exon_filter), axis = 1)
    df_stat_all_transcript = df_stat_all_transcript.merge(df_stat_possible_transcript, how = "left").query('transcript_id==possible_transcript')
    del(df_stat_all_transcript['possible_transcript'])
    return df_stat_all_transcript






def f_combine_pvalue_transcript(transcript_id, gene_name, num_exon_transcript, df_stat_all_region):
    sub_df_stat_all_region = df_stat_all_region[df_stat_all_region.apply(lambda x: (transcript_id in x.transcript_id.split(',')), axis = 1)].query('gene_name==@gene_name')
    list_exon_combined_pvalue_start = sub_df_stat_all_region.bbinom_pvalue_start.tolist()
    list_exon_combined_pvalue_end = sub_df_stat_all_region.bbinom_pvalue_end.tolist()
    if num_exon_transcript > len(list_exon_combined_pvalue_start):
        list_exon_combined_pvalue = list_exon_combined_pvalue_start + [1]*(num_exon_transcript - len(list_exon_combined_pvalue_start)) + list_exon_combined_pvalue_end + [1]*(num_exon_transcript - len(list_exon_combined_pvalue_end))
    else:
        list_exon_combined_pvalue = list_exon_combined_pvalue_start + list_exon_combined_pvalue_end 
    transcript_combined_pvalue = f_combin_p(list_exon_combined_pvalue)
    return(transcript_combined_pvalue)

def f_avg_minus_log10P(transcript_id, gene_name, num_exon_transcript, df_stat_all_region):
    sub_df_stat_all_region = df_stat_all_region[df_stat_all_region.apply(lambda x: (transcript_id in x.transcript_id.split(',')), axis = 1)].query('gene_name==@gene_name')
    list_exon_combined_pvalue_start = sub_df_stat_all_region.bbinom_pvalue_start.tolist()
    list_exon_combined_pvalue_end = sub_df_stat_all_region.bbinom_pvalue_end.tolist()
    if num_exon_transcript > len(list_exon_combined_pvalue_start):
        list_exon_combined_pvalue = list_exon_combined_pvalue_start + [1]*(num_exon_transcript - len(list_exon_combined_pvalue_start)) + list_exon_combined_pvalue_end + [1]*(num_exon_transcript - len(list_exon_combined_pvalue_end))
    else:
        list_exon_combined_pvalue = list_exon_combined_pvalue_start + list_exon_combined_pvalue_end 
    transcript_combined_pvalue = sum([math.log10(x+1E-10) for x in list_exon_combined_pvalue])/len(list_exon_combined_pvalue)
    return(transcript_combined_pvalue)




def f_avg_pvalue(transcript_id, gene_name, num_exon_transcript, df_stat_all_region):
    sub_df_stat_all_region = df_stat_all_region[df_stat_all_region.apply(lambda x: (transcript_id in x.transcript_id.split(',')), axis = 1)].query('gene_name==@gene_name')
    list_exon_combined_pvalue_start = sub_df_stat_all_region.bbinom_pvalue_start.tolist()
    list_exon_combined_pvalue_end = sub_df_stat_all_region.bbinom_pvalue_end.tolist()
    if num_exon_transcript > len(list_exon_combined_pvalue_start):
        list_exon_combined_pvalue = list_exon_combined_pvalue_start + [1]*(num_exon_transcript - len(list_exon_combined_pvalue_start)) + list_exon_combined_pvalue_end + [1]*(num_exon_transcript - len(list_exon_combined_pvalue_end))
    else:
        list_exon_combined_pvalue = list_exon_combined_pvalue_start + list_exon_combined_pvalue_end 
    transcript_combined_pvalue = sum(list_exon_combined_pvalue)/len(list_exon_combined_pvalue)
    return(transcript_combined_pvalue)





def f_avg_cDNA(transcript_id, gene_name, num_exon_transcript, df_stat_all_region):
    sub_df_stat_all_region = df_stat_all_region[df_stat_all_region.apply(lambda x: (transcript_id in x.transcript_id.split(',')), axis = 1)].query('gene_name==@gene_name')
    total_cDNA = sub_df_stat_all_region.num_clipped_start.tolist() + sub_df_stat_all_region.num_clipped_end.tolist()
    if len(total_cDNA)<(2*num_exon_transcript):
        total_cDNA = total_cDNA + [1]*(2*num_exon_transcript - len(total_cDNA))
    total_cDNA = np.array(total_cDNA)
    # x = np.isnan(total_cDNA)
    x = pd.isnull(total_cDNA)
    total_cDNA[x] = 0
    avg_cDNA_per_exon = np.nansum(total_cDNA)/(2*num_exon_transcript)
    return(avg_cDNA_per_exon)


def f_median_cDNA(transcript_id, gene_name, num_exon_transcript, df_stat_all_region):
    sub_df_stat_all_region = df_stat_all_region[df_stat_all_region.apply(lambda x: (transcript_id in x.transcript_id.split(',')), axis = 1)].query('gene_name==@gene_name')
    total_cDNA = sub_df_stat_all_region.num_clipped_start.tolist() + sub_df_stat_all_region.num_clipped_end.tolist()
    if len(total_cDNA)<(2*num_exon_transcript):
        total_cDNA = total_cDNA + [1]*(2*num_exon_transcript - len(total_cDNA))
    total_cDNA = np.array(total_cDNA)
    # x = np.isnan(total_cDNA)
    x = pd.isnull(total_cDNA)
    total_cDNA[x] = 0
    avg_cDNA_per_exon = np.median(total_cDNA)
    return(avg_cDNA_per_exon)





def f_recheck_exon(df_stat_single, df_stat_all_region):
    # df_stat_single = df_stat_all_transcript.iloc[2]
    if (df_stat_single.num_exon_detected) == 1:
        transcript_id = df_stat_single.transcript_id
        df_exon = df_stat_all_region[df_stat_all_region.apply(lambda x: transcript_id in x.transcript_id.split(',') , axis = 1)].iloc[0]
        if f_exon_filter_strict(df_exon, global_para) == False:
            num_exon_detected_adjust = 0
        else:
            num_exon_detected_adjust = 1
    else:
        num_exon_detected_adjust = df_stat_single.num_exon_detected
    return num_exon_detected_adjust



def f_exon_filter_strict(df_exon, global_para):
    if re.search("1", df_exon.is_exon_boundary_start):
        is_exon_match_region_start = True;
        is_pvalue_start = df_exon.bbinom_p_start <= (global_para.cutoff_pvalue/2)
    else:
        is_exon_match_region_start = df_exon.is_exon_match_start;
        is_pvalue_start = df_exon.bbinom_p_start <= (global_para.cutoff_pvalue/2)
    if re.search("1", df_exon.is_exon_boundary_end):
        is_exon_match_region_end = True;
        is_pvalue_end = df_exon.bbinom_p_end <= (global_para.cutoff_pvalue/2)
    else:
        is_exon_match_region_end = df_exon.is_exon_match_end;
        is_pvalue_end = df_exon.bbinom_p_end <= (global_para.cutoff_pvalue/2)
    is_keep = all([is_pvalue_start,is_exon_match_region_start, is_pvalue_end,is_exon_match_region_end])
    return(is_keep)


