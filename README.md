# cDNA-detector
Detect and remove cDNA contamination in DNA sequencing experiments (ATAC-Seq, ChIP-Seq, WES, etc.)
## General overview
`cdna-detector` is a Python tool for detecting and removing contaminating cDNA from vectors and other sources in DNA-Seq data (ATAC-Seq, ChIP-Seq, WES, etc). The tool analyzes BAM formatted alignment files (coordinate-sorted), detects candidate contaminant cDNAs and, if desired, removes suspected contaminant reads from BAM files.

## Prerequisites
`cdna-detector` requires Python3 and blastn. We have tested it in Python 3.7.6 and BLASTN (2.9.0).
`cdna-detector` also requires the Python modules listed below:
```
pandas >= 0.23.0
numpy >= 1.17.0
scipy >= 1.1.0
pysam >= 0.15.0
Statsmodels >= 0.11.1
Biopython >= 1.77
joblib >= 0.15.1
gffpandas >= 0.1.0
```
Modules can be installed using pip or conda. Note that gffpandas cannot be installed via conda, but is only necessary if generating custom gene model files. 
```
pip install --user pandas numpy scipy pysam statsmodels biopython gffpandas joblib
conda install -c bioconda blast==2.9.0
```
or
```
conda install  pandas numpy scipy statsmodels biopython joblib
conda install -c bioconda pysam
conda install -c bioconda blast==2.9.0
pip install gffpandas
```
Gene model files for the hg19 and hg38 human and mm10 mouse genome assemblies are provided with cDNA-detector. The subcommand `prepare` is only necessary when generating gene models for other species, assemblies or custom amplicons.

## Installation
```
git clone https://github.com/rheinbaylab/cDNA-detector.git
cd cDNA-detector
chmod a+x cdna-detector.py
```

Then add cdna-detector into your PATH.


## Usage
`cDNA-detector` consists of three steps (subcommands):

1. `prepare`: Generate gene model file [optional]
2. `detect`: Detect candidate contaminant cDNA in bam file
3. `clean`: Clean up cDNA in bam file [optional]

Details and examples are shown below.


### 1. prepare: Generate gene model file [optional]
`cdna-detector` relies on a species and assembly-specific gene model file to search for possible contaminant genes. We provide gene models for hg19/hg38/mm10 with this distribution. The 'prepare' step is thus only required for other species/assemblies or custom amplicon sets. The gene model file needs to be created only once. 
Built-in gene models are
```
hg19
hg38
hg19_nochr
hg38_nochr
mm10
mm39
```
Gene models `hg19_nochr` and `hg38_nochr` are used for reference genomes where chromosome names do not start with "chr".
#### Human hg19/hg38
Gene models ready for use are provided in diretory `data/gene_model`. 
And hg19/hg38 are equivalent to gene model files in the directory. Please see section 2 for usage examples. 

#### Mouse mm10/mm39
Gene models are generated in diretory `data/gene_model`. 
Users can use these files directly. And mm10/mm39 is equivalent to gene model file in the directory.

#### Other species/assemblies
**usage**
```
cdna-detector.py prepare --annotation <gtf/bed> --genome <genome sequenc filee>  [options]
```
**example**
```
cd cDNA-detector/
cdna-detector  prepare --annotation  example_data/gtf_example/test.gtf --output_dir example_data/output_prepare --genome  example_data/test_genome/test_genome.fa
```

**Input parameters**
```
usage: cdna-detector.py prepare --annotation <gtf/bed> --genome <genome sequence file>  [options]

required arguments:
  --annotation          gene annotation files
                        input format: gtf/bed
  --genome              genome fasta file
                        input format: fa/fasta

optional arguments:
  -h, --help            show this help message and exit
  --format              input annotation format: gtf or bed
                         If format is "bed", only 3 or 4 columns will be used.
                        - default: gtf
  --output_dir          output directory
                        - default: .
  --chr_exclude         exclude chromosomes, multiple chromosomes are separated by ","
                        - conflict with --chr_include
                        - default: chrM
                        - example: --chr_exclude chrX,chrY,chrM
  --chr_include         only include chromosomes, multiple chromosomes are separated by ","
                        - conflict with --chr_exclude
                        - example: --chr_include chr1,chr2
  --source              the program that generated this feature, it is located in 2nd column of gtf. multiple sources are separated by ","
                        - default: all source
                        - example: --source havana
  --feature_type_cds    feature name for identifying CDS regions. It is located in 3rd column of gtf. Multiple features are separated by ","
                        - default: CDS
  --featureid_gene_id   attribute to show gene id. 
                        - default: gene_id
  --featureid_gene_name 
                        attribute to show gene name. 
                        - default: gene_name
  --featureid_cds_ranking 
                        attribute to show exon number.
                        - default exon_number
  --featureid_transcript_id 
                        attribute to show transcript id. 
                        - default transcript_id
```
`annotation` input file can be in GTF or BED format, `genome` is the genome sequence. The chromosome names in the genome sequences must be consistent with chromosome names in GTF/BED files.
If the input file is in BED format, only the first 3 or 4 columns will be required and used. Details can be found in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1.com).

By default, chrM is excluded due to high read accumulation in sequencing experiments. It can be added back to the analysis workflow with `--chr_include chrM`.


**Output files**

1. `GtfFileName.saf`: Gene model file. This file contains gene exon regions and boundary positions. This files is required for step `detect`.
2. `GtfFileName.prepare.log`: Log file.

Here are the first several lines of the file `GtfFileName.saf`.

| seqname | start  | end    | gene_id | gene_name | transcript_id | exon_flank_start20   | exon_flank_end20     | is_exon_boundary_start | is_exon_boundary_end | exon_boundary_start_nearseq20 | exon_boundary_end_nearseq20 |
|---------|--------|--------|---------|-----------|---------------|----------------------|----------------------|------------------------|----------------------|-------------------------------|-----------------------------|
| chr1    | 901912 | 901994 | PLEKHN1 | PLEKHN1   | NM_001160184  | GACCTGTACGACTCTGGCCA | TGAGCGGGGCGTGGGTGCGG | 1                      | 0                    | GAGGACAGCGCGCGGATGTC          |                             |
| chr1    | 901912 | 901994 | PLEKHN1 | PLEKHN1   | NM_001367552  | GACCTGTACGACTCTGGCCA | TGAGCGGGGCGTGGGTGCGG | 1                      | 0                    | GAGGACAGCGCGCGGATGTC          |                             |
| chr1    | 901912 | 901994 | PLEKHN1 | PLEKHN1   | NM_032129     | GACCTGTACGACTCTGGCCA | TGAGCGGGGCGTGGGTGCGG | 1                      | 0                    | GAGGACAGCGCGCGGATGTC          |                             |
| chr1    | 902084 | 902183 | PLEKHN1 | PLEKHN1   | NM_001160184  | CCCCTTGCCTTGTCCCCAGA | TGAGCGCGGCGTGCACGGTG | 0                      | 0                    | CCTCGCTGAAGGGAAACAGG          | ACATCCTGGACCTGGAGAAC        |
| chr1    | 902084 | 902183 | PLEKHN1 | PLEKHN1   | NM_001367552  | CCCCTTGCCTTGTCCCCAGA | TGAGCGCGGCGTGCACGGTG | 0                      | 0                    | CCTCGCTGAAGGGAAACAGG          | ACATCCTGGACCTGGAGAAC        |

Users can provide several required columns for exon regions, such as seqname, start, end, gene_id, gene_name, then fill columns `is_exon_boundary_start` and `is_exon_boundary_end` with 0. Fill transcript with transcript id or gene id (If no transcript id.). Leave other columns blank.


### 2. detect: detecting cDNA contamination in BAM files
`cdna-detector detect` searches for possible cDNA contamination in BAM-formatted sequence alignments using the predefined gene model file. `cdna-detector detect` outputs several files, including original statistics about candidate contaminant exons and genes, high-confidence filtered exons and genes, and coordinate files containing the positions of genes with suspected contaminations for use with the `clean` step.

**Usage**
```
cdna-detector detect --bam <bam> --sample_id <sample_id> --gene_model <gene_model> [options]
```
**Example**
```
cd cDNA-detector/
cdna-detector detect --bam  example_data/bam_file/tmp.sample.bam --sample_id tmp_sample  --gene_model hg19 --output_dir example_data/output_detect
```

**Input options**

```
usage: cdna-detector.py detect --bam <bam> --sample_id <sample_id> --gene_model <gene_model> [options]

required arguments:
  --bam                 The input file is in BAM format. Recommend results from software bwa
  --sample_id           Identifier. Used as output prefix
  --gene_model          Link to gene model.
                         - INPUT: hg19/hg38/mm10/(gene model generated from subcommand "prepare")

optional arguments:
  -h, --help            show this help message and exit
  --min_quality         Minimum read mapping quality.
                        - integer
                        - default: 0
  --pvalue              significant p values.
                        - float
                        - default: 0.05
  --min_ratio_transcript 
                        minimum ratio of detected exons of one transcript
                        - float
                        - default: 0.3
  --output_dir          output directory
                        - default: "."
  --n_threads           number of threads
                        - integer
                        - default: 1
  --min_exon_cdna_reads 
                        minimum number of reads which may come from cDNA for each exon.
                        - integer
                        - default: 1
  --exclude_ehc {True,False}
                        to exclude extremely high coverage of clipped reads for single gene from background calculation.
                        - string
                        - default: True
  --ratio_ehc           cutoff for ratio of extremely high coverage of clipped reads for single gene.
                        - float
                        - default: 0.05
  --count_ehc           cutoff for read count of extremely high coverage of clipped reads for single gene.
                        - integer
                        - default: 10000
  --blastn_database     databases for cDNA source inference.
                        - default: corresponding databases for build-in gene models or default databases
  --file_source_known   file represent gene class.
                        - default: for human, genes from retrocopy and blacklist are listed.
  --inferred_source_include 
                        only include cDNAs with inferred sources, multiple input are separated by ","
                        - conflict with --inferred_source_exclude
                        - example: --inferred_source_include vector,unknown
                        - default: output all cDNAs with inferred sources
  --inferred_source_exclude 
                        only exclude cDNAs with inference sources, multiple input are separated by ","
                        - conflict with --inferred_source_include
                        - example: --inferred_source_exclude retrocopy
                        - default: null value
  --known_source_include 
                        only include cDNAs with known sources, multiple input are separated by ","
                        - conflict with --known_source_exclude
                        - example: --known_source_include retrocopy
                        - default: null value
  --known_source_exclude 
                        only exclude cDNAs with inference sources, multiple input are separated by ","
                        - conflict with --known_source_include
                        - example: --known_source_exclude blacklist
                        - default: blacklist
```

`gene_model` can be a **hg19**/**hg38**/**mm10**/**hg19_nochr**/**hg38_nochr**/**mm39** gene model provided with this tool or gene model file generated by subcommand `cdna-detector prepare` (see above). Note: Adaptor sequences can affect the sensitivity of cDNA-detector. It's recommend to remove adaptor sequences in BAM files.


**Output files**

If cDNAs are detected, `cdna-detector detect` outputs several files. If no cDNA is detected in the input, only the `SAMPLE_ID.log` log file will be generated.
Important files are explained below:

1. `SAMPLE_ID.log`: log file.
2. `SAMPLE_ID.gene_statistics.tsv` original statistics for all tested genes.
3. `SAMPLE_ID.exon_statistics.tsv` original statistics for all tested exons.
4. `SAMPLE_ID.exon_statistics.filter.tsv` filtered results containing only "high-confidence" candidate contaminant exons.
5. `SAMPLE_ID.gene_statistics.filter.tsv` filtered results containing only "high-confidence" candidate contaminant genes.
6. `SAMPLE_ID.gene_statistics.filtered.source_filtered.tsv`: filtered results containing only "high-confidence" candidate contaminant genes from user-defined sources.
7. `SAMPLE_ID.merge_region.tsv`:  merged exon regions which can be used in the "clean" step (see below).
8. `SAMPLE_ID.merge_region.bed`:  merged exon regions which can be used in IGV for manually checking results.
9. `SAMPLE_ID.clipped_seq.source_inference.tsv`: Inferred source (vector or retrogene) for each candidate cDNA, if enough evidence is available.


For most cases, users can only check `SAMPLE_ID.gene_statistics.filtered.source_filtered.tsv`, which describes high-confidence cDNAs in the BAM file. It is recommended to review cDNA calls in IGV with `SAMPLE_ID.merge_region.bed`. 

**Important Note**: cDNAs where only one exon scored as significant can represent a false positive and should be reviewed.

Description of columns in `SAMPLE_ID.gene_statistics.tsv`, `SAMPLE_ID.gene_statistics.filter.tsv` and `SAMPLE_ID.gene_statistics.filtered.source_filtered.tsv`:

| column                | description                                                          |
|---------------------  |----------------------------------------------------------------------|
| gene_id               | gene id                                                              |
| gene_name             | gene name                                                            |
| transcript_id         | transcript id                                                        |
| num_exon_detected     | number of exons detected                                             |
| num_exon_transcript   | number of exons in the transcript                                    |
| ratio                 | ratio of detected exons to all exons in transcript                   |
| source_inference      | possible source of detected cDNAs inferred by cDNA-detector          |
| source_known_databases| known source of detected cDNAs defined by users or built-in databases|
| combined_pvalue       | combined p-value for each transcript                                 |
| avg_log10pvalue       | average log10 pvalue for each transcript                             |
| avg_pvalue            | average pvalue for each transcript                                   |
| avg_cDNA              | average value of possible cDNA reads for each transcript             |
| median_cDNA           | median value of possible cDNA reads for each transcript              |




Description of columns in `SAMPLE_ID.exon_statistics.tsv` and `SAMPLE_ID.exon_statistics.filter.tsv`.

| column              | description                                                         |
|---------------------|---------------------------------------------------------------------|
| seqname             | chromosome name                                                     |
| exon_start          | start position of exon                                              |
| exon_end            | end position of exon                                                |
| pos_start           | start position of candidate region                                  |
| pos_end             | end position of candidate region                                    |
| gene_id             | gene id                                                             |
| gene_name           | gene name                                                           |
| transcript          | transcript name                                                     |
| num_clipped_start   | number of clipped reads overlapping start of candidate region       |
| num_total_start     | number of total reads at start position of candidate region         |
| bbinom_pvalue_start | binomial p-value at start position of candidate region             |
| num_clipped_end     | number of clipped reads at end position of candidate region         |
| num_total_end       | number of total reads at end position of candidate region           |
| bbinom_pvalue_end   | binomial p-value at end position of candidate region                |
| num_bg_clipped      | number of clipped reads in total exon regions (background)          |
| num_bg_total        | number of total reads in total exon regions (background)            |
| combine_pvalue      | combined p-value of start and end of candidate regions              |
| combine_pvalue      | combined q-value of start and end of candidate regions              |






### 3. clean up cDNA in BAM file
`cdna-detector clean` utilizes `SAMPLE_ID.merge_region.tsv` generated in the previous step `cdna-detector detect` and the source BAM file, then generates a "clean" BAM file by removing identified contaminant reads from candidate regions.

**Usage**
```
cdna-detector clean --bam <bam> --region <region_file> --sample_id <sample_id> [options]
```
**Example**
```
cdna-detector clean --bam example_data/bam_file/tmp.sample.bam  --sample_id tmp_sample --region example_data/output_detect/tmp_sample.merge_region.tsv   --output_dir example_data/output_clean
```
**Input options**
```
required arguments:
  --bam                 The input file in bam format
                        - Note: must be sorted (use samtools sort if necessary)
  --region              files which show contaminated regions
                        - output "SAMPLE_ID.merge_region.tsv" from the "detect" step
  --sample_id           Identifier. Used as output prefix

optional arguments:
  -h, --help            show this help message and exit
  --method {automatic,rpm,fraction}
                        method of removing reads from bam files
                        - value: automatic (recommended), fraction, rpm
                        - automatic: estimate cDNA based on cDNA in exon boundaries
                        - fraction: fraction of reads to keep after cleaning in exon regions
                        - rpm: reads per million to keep after cleaning in exon regions
                        - default: automatic
  --output_dir          output directory
                        - default: "."
  --gDNA_remove         by default, reads identified as genomic DNA (gDNA) will be kept in exon regions. If the user suspects all reads in a given exon are contaminant reads, add this flag
```

By default, `cdna-detector clean` estimates the fraction of cDNA in a candidate region based on cDNA detected in exon boundaries. Sometimes, this is not accurate. In this case, the user has the option to set `method` to rpm or fraction to estimate quantity of cDNA manually (i.e after alignment review in IGV). This can be useful when all reads in a candidate exon regions stem from contaminating cDNA.

```
cdna-detector clean --bam example_data/bam_file/tmp.sample.bam  --sample_id tmp_sample --region example_data/output_detect/tmp_sample.merge_region.tsv   --output_dir example_data/output_clean
```

**Output**

`cdna-detector clean` outputs 4 files:

1. `SAMPLEID.clean.log`: log file.
2. `SAMPLEID.clean.bam`: clean BAM file.
3. `SAMPLEID.clean.bam`: index for clean BAM file.
4. `SAMPLEID.clean_region.tsv`: cleaned regions and detailed information about how much cDNA was removed in each region. Useful for manual cDNA removal. [ If `method` is rpm or fraction, this file will not be generated ]

Example `SAMPLEID.clean_region.tsv`.

| seqname | start    | end      | gene_name | rpm                | fraction           |
|---------|----------|----------|-----------|--------------------|--------------------|
| chr1    | 11167436 | 11167558 | MTOR      | 3515.5354398240534 | 0.6                |
| chr1    | 11168234 | 11168344 | MTOR      | 7608.5018214549555 | 0.7173738991192954 |
| chr1    | 11169344 | 11169427 | MTOR      | 5655.426577108259  | 0.7458006718924972 |
| chr1    | 11169705 | 11169787 | MTOR      | 1316.20288205973   | 0.6007751937984496 |

When the automatic method setting appears incorrect, the user can manually edit fraction or rpm. For example, these values can be set to 0 if there is reason to believe all reads in a given exon stem from contaminating cDNA.

Description of columns in `SAMPLEID.clean_region.tsv`

| column              | description                                                         |
|---------------------|---------------------------------------------------------------------|
| seqname             | chromosome name                                                     |
| start               | start position of candidate/cleaned region                          |
| end                 | end position of candidate/cleaned region                            |
| gene_name           | gene name                                                           |
| rpm                 | reads per million kept in the clean bam file                        |
| fraction            | fraction of reads kept in the clean bam file                        |

## Citation
If you use the tool in your work, please cite:
> Qi M, Nayar U, Ludwig LS, Wagle N, Rheinbay E. 2021. [cDNA-detector: Detection and removal of cDNA contamination in DNA sequencing libraries](http://dx.doi.org/10.1101/2021.08.11.455962). bioRxiv doi: 10.1101/2021.08.11.455962

## Getting Help

If you have any questions, please submit issues via GitHub or send email to mqi3@mgh.harvard.edu.

## License Agreement

By using the software, you acknowledge that you agree to the terms below:

For academic and non-profit use, you are free to fork, download, modify, distribute and use the software without restrictions.

For commercial use, you are required to contact the authors at Massachusetts General Hospital to discuss licensing options.

