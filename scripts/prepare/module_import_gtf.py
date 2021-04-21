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


def filename_gff2saf(gtf_file,output_dir):
    # used to show saf file path
    input_full, input_dir, input_id = infile_path(gtf_file)
    out_suffix = '.saf'
    out_gtf_file = outfile_path(output_dir, input_id, out_suffix)
    return(out_gtf_file)


def infile_path(file):
    # used to decode input file
    input_full = os.path.abspath(file)
    input_dir = os.path.dirname(file)
    input_id = os.path.splitext(os.path.basename(file))[0]
    return input_full, input_dir, input_id

def set_logger(out_log):
    # used to record messages from scripts.
    # out_log: log file path
    logging.basicConfig(format= '%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(out_log,'w')
    logger.addHandler(handler)
    handler.setLevel(logging.DEBUG)
    logFormatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(logFormatter)
    return logger

def outfile_path(output_dir, output_id,out_suffix):
    # used to show full output file
    output_dir = os.path.abspath(output_dir)
    out_filename = output_id + out_suffix
    output_full = os.path.join(output_dir,out_filename)
    return output_full



def f_import_transform_bed(global_para):
    with open(global_para.gtf_file) as f:
        first_line = f.readline()
    f.close()
    if first_line.startswith('track') or first_line.startswith('#'):
        df_gff = pd.read_csv(global_para.gtf_file,sep = '\t', skiprows = 1, header = None)
    else:
        df_gff = pd.read_csv(global_para.gtf_file,sep = '\t' , header = None)
    ## add more column names for bed format
    if df_gff.shape[1] >=3:
        df_gff_cds_expand = df_gff.iloc[:,0:4]
        df_gff_cds_expand.columns = ['seq_id','start','end','gene_name']
    elif df_gff.shape[1] == 3:
        df_gff_cds_expand = df_gff.iloc[:,0:3]
        df_gff_cds_expand.columns = ['seq_id','start','end']
        df_gff_cds_expand['gene_name'] = df_gff_cds_expand.apply(lambda x:"region|" + x.seq_id + ":" + str(x.start) + "-" + str(x.end), axis = 1)
    else:
        logger.error('wrong bed file')
        exit(1)
    df_gff_cds_expand['gene_id'] = df_gff_cds_expand['gene_name']
    df_gff_cds_expand['transcript_id'] = df_gff_cds_expand['gene_name']
    return df_gff_cds_expand



def f_import_transform_gtf(global_para):
    df_gff = gffpd.read_gff3(global_para.gtf_file)
    df_gff_all = df_gff.df
    df_gff_all = df_gff_all.astype({"seq_id":"str"})
    df_gff_cds = df_gff_all.query('type == @global_para.feature_type_cds')
    if global_para.source  != "all":
        tmp_list_source = global_para.source.split(',')
        df_gff_cds = df_gff_cds.query('source in @tmp_list_source')
    df_gff_cds_expand = f_gtf_attributes_expand(df_gff_cds)
    df_gff_cds_expand.rename(columns = {global_para.featureid_gene_id :'gene_id', global_para.featureid_gene_name :'gene_name', global_para.featureid_transcript_id :'transcript_id'},inplace = True)
    return df_gff_cds_expand


def f_read_genome_seq(genome_fasta):
    # get genome chromosome sequences into one dictionary
    genome_dict = dict()
    for seq_record in SeqIO.parse(genome_fasta, "fasta"):
        genome_dict[seq_record.id] = seq_record.seq
    return genome_dict


def get_chromosome_keep(list_chr_all, chr_exclude, chr_include):
    # used to get chromosome list which should be remained.
    if chr_include != "":
        list_chr_include = chr_include.split(',')
        list_chr_keep = list_chr_include
    elif chr_exclude == "":
        list_chr_keep = list_chr_all;
    else:
        list_chr_remove = chr_exclude.split(',')
        list_chr_keep = [x for x in list_chr_all if x not in list_chr_remove]
    return(list_chr_keep)

def f_is_exon_boundary_start(id_order, min):
    # define if exon boundary is the most left of transcript
    if id_order == min:
        return(1)
    else:
        return(0)


def f_is_exon_boundary_end(id_order, max):
    # define if exon boundary is the most right of transcript
    if id_order == max:
        return(1)
    else:
        return(0)

### filter chromosome list
def filter_chromosome(df_gff_exon_rename, list_chr_keep):
    # 5. map other results with ncbi seqname
    df_gff_exon_rename = df_gff_exon_rename[df_gff_exon_rename['seq_id'].isin(list_chr_keep)]
    # df_gff_exon_seqname.to_csv('ncbi_hg38_cds_exon_bestref_anno.tsv',sep = '\t',index = False)
    return(df_gff_exon_rename)


def f_gtf_attributes_expand(df_gff_cds):
    # expand attributes
    df_gff_exon_cds = df_gff_cds.copy()
    list_temp = ['']*len(df_gff_exon_cds)
    d_ref = dict()
    list_key = list();
    list_value = list();
    list_num = list();
    i = 0;    
    for attribute in df_gff_exon_cds['attributes']:
        # print(attribute)
        tmp_list_key = [re.sub('^ ','',x).split(' ',maxsplit = 1)[0] for x in re.sub(';$','',attribute).split(';')]
        tmp_list_value =  [re.sub('^ ','',x).split(' ',maxsplit = 1)[1].strip('"') for x in re.sub(';$','',attribute).split(';')]
        tmp_list_num = [i]*len(tmp_list_key)
        list_key.extend(tmp_list_key)
        list_value.extend(tmp_list_value)
        list_num.extend(tmp_list_num)
        i = i+1
    list_key_unique = list(set(list_key))
    d_temp = {i:'' for i in range(0,len(df_gff_exon_cds))}
    for key in list_key_unique:
        d_ref[key] = d_temp.copy()
    for i in range(0,len(list_key)):
        d_ref[list_key[i]][list_num[i]] = list_value[i]
    d_ref_list = dict.fromkeys(d_ref.keys(),[])
    for name in d_ref_list.keys():
        df_gff_exon_cds[name] = list(d_ref[name].values())
    return(df_gff_exon_cds)




def f_get_seq(seqname, start,end, genome_dict):
    # get genome sequences within one region
    region_seq = str(genome_dict[seqname][start:end])
    return region_seq



def f_exon_boundary_info(df_gff_cds_expand_filter, genome_dict):
    # used to add if boundary the most left/right of a transcript, and nearby exon sequences.
    if global_para.featureid_cds_ranking in df_gff_cds_expand_filter.columns.tolist():
        list_df_id_order = []
        df_gff_cds_expand_filter.rename(columns = {global_para.featureid_cds_ranking :'exon_number'}, inplace = True)
        df_gff_cds_expand_filter = df_gff_cds_expand_filter.astype({"exon_number":"int"})
        for transcript_id,sub_df in df_gff_cds_expand_filter.sort_values('exon_number').groupby(['transcript_id']):
            # print(transcript_id)
            if len(sub_df) == 1:
                sub_df_new = sub_df.copy()
                sub_df_new['id_order'] = 1
            else:
                sub_df_sort = sub_df.sort_values('exon_number')
                if sub_df_sort.iloc[0]['start'] >sub_df_sort.iloc[1]['start']:
                    sub_df_new = sub_df_sort.iloc[::-1]
                else:
                    sub_df_new = sub_df_sort.copy()
                sub_df_new['id_order'] = range(len(sub_df_new))
            list_df_id_order.append(sub_df_new)
        df_gff_cds_expand_filter = pd.concat(list_df_id_order, axis = 0)
    else:
        rank = df_gff_cds_expand_filter.groupby(['transcript_id'])['start'].rank(ascending = True, method = 'first')
        rank.name = "id_order"
        df_gff_cds_expand_filter = pd.concat([df_gff_cds_expand_filter,rank],axis = 1)
    df_gff_cds_expand_filter = df_gff_cds_expand_filter.astype({"id_order":"int"})
    # b. add boundary information: is start/end of one transcripts
    tmp_df_start_min_max = df_gff_cds_expand_filter.groupby('transcript_id')['id_order'].agg([('id_order_min','min'),('id_order_max','max')]).reset_index
    df_gff_cds_expand_filter = df_gff_cds_expand_filter.merge(tmp_df_start_min_max())
    df_gff_cds_expand_filter['is_exon_boundary_start'] = df_gff_cds_expand_filter.apply(lambda x:f_is_exon_boundary_start(x.id_order, x.id_order_min),axis = 1)
    df_gff_cds_expand_filter['is_exon_boundary_end'] = df_gff_cds_expand_filter.apply(lambda x:f_is_exon_boundary_end(x.id_order, x.id_order_max),axis = 1)
    # c. add boundary information: neary seq information from neiborhood exon
    df_gff_cds_expand_filter['exon_seq'] =  df_gff_cds_expand_filter.apply(lambda x:f_get_seq(x.seq_id,x.start -1, x.end, genome_dict),axis = 1)
    df_gff_cds_expand_filter['exon_seq_start20'] =  df_gff_cds_expand_filter.apply(lambda x:x.exon_seq[0:20],axis = 1)
    df_gff_cds_expand_filter['exon_seq_end20'] =  df_gff_cds_expand_filter.apply(lambda x:x.exon_seq[-20:],axis = 1)
    df_gff_cds_expand_filter['exon_id_order_shift_minus1'] = df_gff_cds_expand_filter.apply(lambda x:x.id_order -1, axis = 1)
    tmp_start_id_order_shift_seq = df_gff_cds_expand_filter.filter(['transcript_id','id_order','exon_seq_end20'])
    x = df_gff_cds_expand_filter.merge(tmp_start_id_order_shift_seq, how = 'left',left_on = ['transcript_id','exon_id_order_shift_minus1'], right_on = ['transcript_id','id_order'],suffixes = ('',"_1"),sort = False)
    x = x.filter(['exon_seq_end20_1']).fillna('')
    df_gff_cds_expand_filter['exon_boundary_start_nearseq20']  = x['exon_seq_end20_1']
    df_gff_cds_expand_filter['exon_id_order_shift_add1'] = df_gff_cds_expand_filter.apply(lambda x:x.id_order +1, axis = 1)
    tmp_end_id_order_shift_seq = df_gff_cds_expand_filter.filter(['transcript_id','id_order','exon_seq_start20'])
    x = df_gff_cds_expand_filter.merge(tmp_end_id_order_shift_seq, how = 'left', right_on = ['transcript_id','id_order'],left_on = ['transcript_id','exon_id_order_shift_add1'],suffixes = ('',"_1"),sort = False)
    x = x.filter(['exon_seq_start20_1']).fillna('')
    df_gff_cds_expand_filter['exon_boundary_end_nearseq20']  = x['exon_seq_start20_1']
    # return df with seq information
    return df_gff_cds_expand_filter


def f_exon_flank_seq(df_gff_cds_expand_filter, genome_dict):
    df_gff_cds_expand_filter['exon_flank_start20'] =  df_gff_cds_expand_filter.apply(lambda x:f_get_seq(x.seq_id,x.start -21, x.start-1, genome_dict),axis = 1)
    df_gff_cds_expand_filter['exon_flank_end20'] =  df_gff_cds_expand_filter.apply(lambda x:f_get_seq(x.seq_id,x.end , x.end + 20, genome_dict),axis = 1)
    return df_gff_cds_expand_filter




def f_gff_exon_to_saf(df_gff_cds_expand_filter):
    # convert gff dataframe into selected saf format
    # remove additional columns
    # remove duplicated rows
    # columns need to keep
    df_gff_cds_expand_filter.rename(columns = {"seq_id":"seqname"},inplace = True)
    list_attribute_keep = ['seqname','start','end','gene_id','gene_name','transcript_id','exon_flank_start20','exon_flank_end20' , 'is_exon_boundary_start', 'is_exon_boundary_end','exon_boundary_start_nearseq20','exon_boundary_end_nearseq20']
    list_attribute_sort = ['seqname','start','end','gene_id','transcript_id']
    df_gtf_cds_expand_filter_redup = df_gff_cds_expand_filter.filter(list_attribute_keep).drop_duplicates().sort_values(list_attribute_sort)
    return df_gtf_cds_expand_filter_redup




def f_sequence_upper_column(df_gtf_unique, list_column_upper):
    for column in list_column_upper:
        df_gtf_unique[column] = df_gtf_unique[column].str.upper()

    return df_gtf_unique