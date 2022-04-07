# diffmeth_run_params.sh: parameters to define differential
# methylation run type and mapped files for the DMAP diffmeth
# program.

# identify if this was RRBS or WGBS: this determines a number
# of downstream operations, including adaptor trimming and
# differential methylation calculations.  Values are "RRBS" or "WGBS"
# where:
# "RRBS" compares methylation on fragments (usually MspI) in size range
# "WGBS" compares methylation on defined windows or regions
#
# this is previously specified for the adaptor trimming and mapping step
# and obviously the values here should be consistent with those.

dmap_run_type="RRBS";
# dmap_run_type="WGBS";

# type of diffmeth run: values are one of "pairwise" "chisquare" "Anova" "CpGlist" or "binlist"
# where:
# "pairwise"  runs Fisher's Exact statistic on two samples
# "chisquare" runs a Chi Square test on multiple samples to identify fragments/regions with greatest variation
# "Anova"    run Analysis of Variance on two sets of multiple samples, producing F statistic
# "CpGlist"   lists methylation counts for each CpG
# "binlist"   lists methylation information for each fragment/region

diffmeth_run_type="Anova";

Anova)
# Anova requires a group number to define how many groups are involved
# and a list of at least 2 mapping output files for each group
# a further parameter defines the amount of information returned, being
# one of: "anova_simple" "anova_medium" or "anova_full"
# where
# "anova_simple" returns F statistic and its probability
# "anova_medium" adds sample counts and greater methylated group
# "anova_full"   further adds methylation figures for each CpG in fragment/region

anova_group_number=3;

file_list1=(
 "./rundata/comb_D1.bam"
 "./rundata/comb_D2.bam"
 "./rundata/comb_D3.bam"
 "./rundata/comb_D4.bam"
);

file_list2=(
 "./rundata/comb_D5.bam"
 "./rundata/comb_D6.bam"
 "./rundata/comb_D7.bam"
);

file_list3=(
 "./rundata/comb_D8.bam"
 "./rundata/comb_D9.bam"
 "./rundata/comb_D10.bam"
)

anova_detail="anova_full";

# fold difference is not relevant for > 2 anova groups

anova_fold_difference="no";

# diffmeth requires a fragment size specification, even for WGBS (ignored then)

diffmeth_fragment_min=40;
diffmeth_fragment_max=220;

# "WGBS" requires a window (tile) length

wgbs_window=1000;

# output file name for diffmeth

diffmeth_out_name="diffmeth_output_anova3.txt";

# threshold criteria for including fragment/regions

# No. of hits per CpG for fragment/region to be valid
# default is 1

min_cpg_hits=10;

# No. of CpG of fragment/region that must meet min_cpg_hits criterion
# 0 indicates no check for this

min_valid_cpgs=0;

# Minimum No. samples reaching criteria for a fragment/region to be included
# defaults to 2

min_sample_count=6;

# For RRBS 3' mapping reads, count initial CpG to previous fragment.

map_init_CpG_to_prev="yes";

# Probability threshold: omit fragments/regions with Pr greater.
# default (1.0) implies no filtering

#pr_threshold=1.0;

# Max read length for .SAM files (not needed for .BAM, though)
# defaults to 150bp, for 0

#max_sam_read_length=150;

# Vertebrate C methylation is largely confined to CpG dinucleotides,
# invertebrates and plants less so, so an option is needed for non_CpG methylation
# to be enabled.

non_CpG_methylation="";

