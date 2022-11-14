# sample_params.sh: parameters to descibe a single bisulphite
# sequence sample, including the raw fastq file(s) and any
# trimming that might be justified by the read qualities.

# identify if this was RRBS or WGBS: this determines a number
# of downstream operations, including adaptor trimming and
# differential methylation calculations.  Values are RRBS or WGBS.

dmap_run_type="WGBS";

# name of adaptor-trimmed output dir:

adtrimmed_out_dir="./";

# sample file name(s): if only maping read 1 then
# leave read 2 name blank/empty.

dmap_sample_files=("SRR18283145_1.fastq.gz"
		   "SRR18283145_2.fastq.gz");

# hard trim length, leave 0 for no trimming.

dmap_sample_trim_length=0;

# minimum read length for retention: usually set this to 20

dmap_min_read_length=20;

# mapping output directory

mapping_out_dir="./";
