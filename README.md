#				DMAP2

A package of scripts to simplify use of the DMAP differential
methylation analysis package (https://github.com/peterstockwell/DMAP).

DMAP operates as a series of steps:
```
a. Genome preparation for bismark mapping

b. Adaptor trimming of bisulphite sequence reads for each sample

c. Mapping reads for each sample to the genome with bismark producing
.bam or .sam files

d. Comparing mapped data with the DMAP program diffmeth

e. Relating genomic fragments or regions of interest to genomic
features with the DMAP program identgeneloc
```
DMAP2 organises these steps with a series of bash scripts using bash
variables to control the analysis.

If prior work has already performed some steps (genome preparation,
mapping) then those can be omitted and appropriate information put in
the basic parameter file.  While the mapping scripts are based on
using the bismark bisulphite aligner, other aligners can be used
(e.g. bsmap) in which case they can be processed from the diffmeth
point at step (d).

Platforms: both DMAP and DMAP2 have been developed and tested in MacOS
X and RH Linux environments.  DMAP is written in C and should compile
under any typical Linux/Unix systems, while DMAP2 should operate with
any typical BASH shell.  RAM requirements will vary depending on
genome sizes and numbers of sequence reads involved in a project,
5-10Gb of available RAM should suffice for moderately complex work.

This distribution can by obtained either by downloading from
https://github.com/peterstockwell/DMAP2 and unzipping or by command
line download with:

git clone https://github.com/peterstockwell/DMAP2

In either case the resulting directory will contain:
```
README.md - this file
DMAP2_UserGuide.pdf - documentation
Scripts - directory containing the main DMAP2 scripts
Example_params - directory containing example parameter scripts for
                   the various steps
```
The Scripts directory contains:
```
dmap_index_build.sh - to generate bismark indices and organise the
                        genomic annotation
map_bisulphite_reads.sh - to adaptor trim reads and map with bismark
dmap_run_diffmeth.sh - to calculate differential methylation with the
                        DMAP program diffmeth
dmap_run_genloc.sh - to relate sequence regions or fragments to
                        genomic features
```
The Example_params directory contains:
```
dmap_basic_params.conf - to provide basic project details to all the
                           scripts in Scripts
diffmeth_run_params.conf - diffmeth config file in DMAP2_UserGuide.pdf
                             Appendix III
sample_params.conf - typical sample parameters for the mapping step for
                     map_bisulphite_reads.sh
dmap_anova2_params.conf - parameters for diffmeth with 2 group Anova
                     statistic
dmap_anova3_params.conf - for 3 group Anova statistic
dmap_chisq_params.conf - for differential methylation analysis by Chi
                     Square statistic
dmap_prwise_params.conf - for Fisher's Exact pairwise analysis
dmap_cpglist_params.conf - list CpG counts and positions for a single
                     sample for each fragment or region
dmap_binlist_params.conf - list methylation counts for a single sample
                     for each fragment or bin
genloc_params.conf - parameters for identgeneloc to relate
                     fragments/regions to genomic features

And parameter files for SRA examples:

srr18283145_sample_params.conf		- mapping parameters for SRR18283145.
srr_145_146_diffmeth_pairwise.conf	- parameters for pairwise diffmeth run.
srr_145_146_genloc_params.conf		- parameters for identgeneloc on diffmeth output
```
Documentation: in the file DMAP2_UserGuide.pdf and in comments in the
script files.
