# srr_145_146_diffmeth_params.sh: parameters for a pairwise differential
# methylation run for samples SRR18283145 & SRR18283146

# are samples RRBS or WGBS?

dmap_run_type="WGBS";

# type of diffmeth run: values are one of "pairwise" "chisquare" "Anova" "CpGlist" or "binlist"

diffmeth_run_type="pairwise";

# "pairwise" requires two mapping output files

file_list1=("SRR18283145_1_at.fastq_bsmapz.bam"
            "SRR18283146_1_at.fastq_bsmapz.bam");

# diffmeth requires a fragment size specification, even for WGBS (ignored then)

diffmeth_fragment_min=40;
diffmeth_fragment_max=220;

# "WGBS" requires a window (tile) length

wgbs_window=1000;

# output file name for diffmeth

diffmeth_out_name="srr_145_146_pairwise_diffmeth_bsmap.txt";

# threshold criteria for including fragment/regions

# No. of hits per CpG for fragment/region to be valid
# default is 1

min_cpg_hits=1;

# No. of CpG of fragment/region that must meet min_cpg_hits criterion
# 0 indicates no check for this

min_valid_cpgs=0;

# Minimum No. samples reaching criteria for a fragment/region to be included
# defaults to 2

min_sample_count=2;

# Vertebrate C methylation is largely confined to CpG dinucleotides,
# invertebrates and plants less so, so an option is needed for non_CpG methylation
# to be enabled.

non_CpG_methylation="";

