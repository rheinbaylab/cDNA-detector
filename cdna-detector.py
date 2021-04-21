#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import os
import sys
import re
from argparse import RawTextHelpFormatter

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def prepare(args):
	from scripts.prepare import step1
	from scripts.prepare import global_para
	global_para.gtf_file = args.annotation
	global_para.output_dir = args.output_dir
	global_para.chr_exclude = args.chr_exclude
	global_para.chr_include = args.chr_include
	global_para.genome_fasta = args.genome
	global_para.format = args.format
	global_para.source = args.source
	global_para.feature_type_cds = args.feature_type_cds
	global_para.featureid_gene_id = args.featureid_gene_id
	global_para.featureid_gene_name = args.featureid_gene_name
	global_para.featureid_cds_ranking = args.featureid_cds_ranking
	global_para.featureid_transcript_id = args.featureid_transcript_id

	step1.convert_gff2saf(global_para)


def detect(args):
	from scripts.detect import step2
	from scripts.detect import global_para
	global_para.script_path = get_script_path()
	global_para.gtf_gene_unique_file, global_para.blastn_database, global_para.file_source_known = global_para.f_get_input_files(args, global_para)	
	global_para.genome_bam_file = args.bam
	global_para.sampleid = args.sample_id
	global_para.min_quality = int(args.min_quality)
	global_para.output_dir = args.output_dir
	global_para.n_threads = int(args.n_threads)
	global_para.cutoff_num_exon_unaligned_reads = args.min_exon_cdna_reads
	global_para.cutoff_pvalue = args.pvalue
	global_para.cutoff_ratio_gene = args.min_ratio
	global_para.exclude_ehc = True if args.exclude_ehc == "True" else False
	global_para.ratio_ehc = args.ratio_ehc
	global_para.count_ehc = args.count_ehc
	global_para.source_inference_include = args.inferred_source_include
	global_para.source_inference_exclude = args.inferred_source_exclude
	global_para.source_known_databases_include = args.known_source_include
	global_para.source_known_databases_exclude = args.known_source_exclude
	# output files
	global_para.out_exon_stat = os.path.join(global_para.output_dir, (global_para.sampleid + '.exon_statistics.tsv'))
	global_para.out_gene_stat = os.path.join(global_para.output_dir, (global_para.sampleid + '.gene_statistics.tsv'))
	global_para.out_exon_stat_filter = os.path.join(global_para.output_dir, (global_para.sampleid + '.exon_statistics.filtered.tsv'))
	global_para.out_gene_stat_filter = os.path.join(global_para.output_dir, (global_para.sampleid + '.gene_statistics.filtered.tsv'))
	global_para.out_gene_stat_filter_source = os.path.join(global_para.output_dir, (global_para.sampleid + '.gene_statistics.filtered.source_filtered.tsv'))

	global_para.out_bed_merge = os.path.join(global_para.output_dir, (global_para.sampleid + '.merge_region.tsv'))
	global_para.out_log = os.path.join(global_para.output_dir, (global_para.sampleid + '.log'))
	global_para.out_blastn_seq = os.path.join(global_para.output_dir, (global_para.sampleid + '.clipped_seq'))
	global_para.out_blastn_seq_table = global_para.out_blastn_seq + ".blastn.tsv"
	global_para.out_blastn_seq_source = global_para.out_blastn_seq + ".source_inference.tsv"


	step2.detect_cdna(global_para)

def clean(args):
	from scripts.clean import step3
	from scripts.clean import global_para
	# global_para = global_para;
	global_para.file_region = args.region
	global_para.file_bam = args.bam
	# global_para.genome_bam_file = args.bam
	global_para.sampleid = args.sample_id
	global_para.bl_gDNA_remove = args.gDNA_remove
	global_para.clean_way = args.clean_way
	global_para.output_dir = args.output_dir
	global_para.file_bam_clean = os.path.join(global_para.output_dir, (global_para.sampleid + '.clean.bam'))
	global_para.file_region_clean = os.path.join(global_para.output_dir,(global_para.sampleid+'.clean_region.tsv'))
	global_para.out_log = os.path.join(global_para.output_dir, (global_para.sampleid + '.clean.log'))

	step3.cdna_remove_region(global_para)
	




def main():
	description_str = '''
cdna-detector.py is a tool to detect and clean cDNA contamination in DNA-Seq.

Version: 0.1.0
Code: https://github.com/rheinbaylab/cDNA_detector.git
Mail: mqi3@mgh.harvard.edu

Usage:
	cdna-detector.py <subcommand> [options]

Example:
	cdna-detector.py prepare -h
	cdna-detector.py detect -h
	cdna-detector.py clean -h
	'''
	parser = argparse.ArgumentParser(usage = argparse.SUPPRESS,add_help = True,description =  description_str,formatter_class=RawTextHelpFormatter)
	subparsers = parser.add_subparsers(title = "sub-commands include",metavar = "-------------------")
	# create the parser for the "step1: create gene model" command
	parser_a = subparsers.add_parser('prepare', help='prepare a gene model annotation file. Input file: NCBI gff3',usage = "cdna-detector.py prepare --annotation <gtf/bed> --genome <genome sequence>  [options]",add_help = True, formatter_class=argparse.RawTextHelpFormatter)
	parser_a.set_defaults(func=prepare)
	required_a = parser_a.add_argument_group('required arguments')
	required_a.add_argument('--annotation', metavar = "", help="gene annotation files\ninput format: gtf/bed\n",required = True)
	required_a.add_argument('--genome', type=str,  metavar = '', help='genome fasta file\ninput format: fa/fasta', required = True)
	parser_a.add_argument('--format', type=str,  metavar = '', help='input annotation format: gtf or bed\n If format is "bed", only 3 or 4 columns will be used.\n- default: gtf', default = "gtf")
	parser_a.add_argument('--output_dir', type=str,  metavar = '', help='output directory\n- default: .',default = '.')
	group = parser_a.add_mutually_exclusive_group()
	group.add_argument("--chr_exclude",  metavar = '', help = 'exclude chromosomes, multiple chromosomes are separated by ","\n- conflict with --chr_include\n- default: chrM\n- example: --chr_exclude chrX,chrY,chrM',default = "chrM")
	group.add_argument("--chr_include",  metavar = '', help = 'only include chromosomes, multiple chromosomes are separated by ","\n- conflict with --chr_exclude\n- example: --chr_include chr1,chr2',default = '')
	parser_a.add_argument('--source',metavar = '', default = "all", help = 'the program that generated this feature, it is located in 2nd column of gtf. multiple sources are separated by ","\n- default: all source\n- example: --source havana')
	parser_a.add_argument('--feature_type_cds',metavar = '', default = "CDS", help = 'feature name for identifying CDS regions. It is located in 3rd column of gtf. Multiple features are separated by ","\n- default: CDS')
	parser_a.add_argument('--featureid_gene_id',metavar = '', default = "gene_id", help = 'attribute to show gene id. \n- default: gene_id')
	parser_a.add_argument('--featureid_gene_name',metavar = '', default = "gene_name", help = 'attribute to show gene name. \n- default: gene_name')
	parser_a.add_argument('--featureid_cds_ranking',metavar = '', default = "exon_number", help = 'attribute to show exon number.\n- default exon_number')
	parser_a.add_argument('--featureid_transcript_id',metavar = '', default = "transcript_id", help = 'attribute to show transcript id. \n- default transcript_id')
	parser_a._action_groups.reverse()

	# create the parser for the "step2: detect cDNA in DNA sequencing datasets" command
	parser_b = subparsers.add_parser('detect', help='detect possible cDNA regions',usage = "cdna-detector.py detect --bam <bam> --sample_id <sample_id> --gene_model <gene_model> [options]",add_help = True, formatter_class=argparse.RawTextHelpFormatter)
	parser_b.set_defaults(func=detect)
	required_b = parser_b.add_argument_group('required arguments')
	required_b.add_argument('--bam',metavar = "",  help='The input file is in BAM format. Recommend results from software bwa',required = True)
	required_b.add_argument('--sample_id',metavar = "",  help='Identifier. Used as output prefix',required = True)
	required_b.add_argument('--gene_model',metavar = "", help='Link to gene model.\n - INPUT: hg19/hg38/mm10/(gene model generated from subcommand "prepare")')
	parser_b.add_argument('--min_quality',metavar = "", type = int, default = 0,  help='Minimum read mapping quality.\n- integer\n- default: 0\n')
	parser_b.add_argument('--pvalue',metavar = "", type = float, default = 0.05,  help='significant p values.\n- float\n- default: 0.05')
	parser_b.add_argument('--min_ratio_transcript',dest = "min_ratio",metavar = "", type = float, default = 0.3, help = "minimum ratio of detected exons of one transcript\n- float\n- default: 0.3")
	parser_b.add_argument('--output_dir',metavar = "", default = '.', help = 'output directory\n- default: "."')
	parser_b.add_argument('--n_threads',metavar = "", type = int, default = 1, help = 'number of threads\n- integer\n- default: 1')
	parser_b.add_argument('--min_exon_cdna_reads',metavar = "", type = int, default = 3,  help='minimum number of reads which may come from cDNA for each exon.\n- integer\n- default: 1')
	parser_b.add_argument('--exclude_ehc',default = "True", choices = ["True", "False"], help='to exclude extremely high coverage of clipped reads for single gene from background calculation.\n- string\n- default: True')
	parser_b.add_argument('--ratio_ehc',metavar = "", type = float, default = 0.05,  help='cutoff for ratio of extremely high coverage of clipped reads for single gene.\n- float\n- default: 0.05')
	parser_b.add_argument('--count_ehc',metavar = "", type = float, default = 10000,  help='cutoff for read count of extremely high coverage of clipped reads for single gene.\n- integer\n- default: 10000')
	parser_b.add_argument('--blastn_database', dest = "blastn_database", metavar = "", default = "", help = "databases for cDNA source inference.\n- default: corresponding databases for build-in gene models or default databases")
	parser_b.add_argument('--file_source_known', dest = "file_source_known", metavar = "", default = "", help = "file represent gene class.\n- default: for human, genes from retrocopy and blacklist are listed.")
	group_b_inference_source = parser_b.add_mutually_exclusive_group()
	group_b_inference_source.add_argument("--inferred_source_include", metavar = "", help = 'only include cDNAs with inferred sources, multiple input are separated by ","\n- conflict with --inferred_source_exclude\n- example: --inferred_source_include vector,unknown\n- default: output all cDNAs with inferred sources', default = '')
	group_b_inference_source.add_argument("--inferred_source_exclude", metavar = "", help = 'only exclude cDNAs with inference sources, multiple input are separated by ","\n- conflict with --inferred_source_include\n- example: --inferred_source_exclude retrocopy\n- default: null value', default = 'retrocopy')
	group_b_known_source = parser_b.add_mutually_exclusive_group()
	group_b_known_source.add_argument("--known_source_include", metavar = "", help = 'only include cDNAs with known sources, multiple input are separated by ","\n- conflict with --known_source_exclude\n- example: --known_source_include retrocopy\n- default: null value', default = '')
	group_b_known_source.add_argument("--known_source_exclude", metavar = "", help = 'only exclude cDNAs with inference sources, multiple input are separated by ","\n- conflict with --known_source_include\n- example: --known_source_exclude blacklist\n- default: blacklist', default = 'blacklist')
	parser_b._action_groups.reverse()


	# create the parser for the "step3: remove cDNA in DNA sequencing datasets" command
	parser_c = subparsers.add_parser('clean', help='remove possible cDNA from bam files',usage = "cdna-detector.py clean --bam <bam> --region <region_file> --sample_id <sample_id> [options]",add_help = True, formatter_class=argparse.RawTextHelpFormatter)
	parser_c.set_defaults(func=clean)
	required_c = parser_c.add_argument_group('required arguments')
	required_c.add_argument('--bam', metavar = "", help='The input file is bam format\n- Note: must be sorted',required = True)
	required_c.add_argument('--sample_id',metavar = "",  help='Identifier. Used as output prefix',required = True)
	required_c.add_argument('--region', metavar = "", help='files which show contaminated regions\n- which should be output of subcommand "detect"',required = True)
	group_c = parser_c.add_mutually_exclusive_group()
	group_c.add_argument('--method',dest = "clean_way",default = "automatic", choices = ["automatic","rpm","fraction"], help='method of removing reads from bam files\n- value: automatic, fraction, rpm\n- automatic: estimate cDNA based on cDNA in exon boundaries\n- fraction: fraction of reads kept after cleaning in exon regions\n- rpm: reads per million kept after cleaning in exon regions\n- default: automatic')
	parser_c.add_argument('--output_dir', default = os.getcwd(), metavar = "", help = 'output directory\n- default: "."')
	# parser_c.add_argument('--paired', default = 'False', action = 'store_true', help='if bam files are paired-end, if not set, bam files are recognized as single-end files')
	parser_c.add_argument('--gDNA_remove', default = 'False', action = 'store_true', help='if set, gDNA will be removed if satisfied cutoff')
	parser_c._action_groups.reverse()
	return(parser)

if __name__ == "__main__":
	parser = main()
	parser._action_groups.reverse()
	args = parser.parse_args()
	if len(sys.argv)==1:
		parser.print_help()
	else:
		args.func(args)

