# dmap_basic_params.sh: a series of definitions that are fundamental to running
# the DMAP package.  These assume that DMAP has already been downloaded
# (git clone https://github.com/peterstockwell/DMAP) and compiled according to
# enclosed instructions.

# verbose: set to "yes" for informational messages to stdout

verbose="yes";

# location of DMAP executables, leave empty if these are installed in your
# PATH and will run by typing their names.  The command 'which diffmeth' or
# 'diffmeth -h' will indicate if this is the case.

path_to_dmap="";

# location of source genome files and files giving chromosome seq
# information.  Terminate directory specs with '/'

dmap_genome_dir="/Users/stope71p/Documents/peter_sw/cpg/dmap_scripts/hs_genome_GRCh37/";
#dmap_genome_dir="/Users/stope71p/Documents/peter_sw/cpg/dmap_scripts/hs_genome_GRCh37/whole_genome/";

# list of genome fasta files present in above location.  DMAP will use each chromosome in a separate file
# or the complete genome in a single file.  Should be a .fa or .fasta extension.

dmap_genome_fasta_files=(
"Homo_sapiens.GRCh37.65.dna.chromosome.1.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.10.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.11.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.12.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.13.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.14.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.15.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.16.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.17.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.18.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.19.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.2.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.20.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.21.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.22.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.3.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.4.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.5.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.6.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.7.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.8.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.9.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.MT.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.X.fa"
"Homo_sapiens.GRCh37.65.dna.chromosome.Y.fa"
);

# As single file:

#dmap_genome_fasta_files=("whole_genome/Homo_sapiens_GRCh37.65.dna.allchromosomes.fa");

# Bisulphite mapper to be used: choices are bismark or bsmapz

bisulphite_mapper="bismark"
#bisulphite_mapper="bsmapz"

case $bisulphite_mapper in

  bismark)

# set stuff for bismark operation:
# location of bismark index directory, to be used for building or using the index. Is the same as the
# genome file location to which you need read and write access.

dmap_bismark_index_location="${dmap_genome_dir}";

# path to bismark executables: leave empty if these are already on your path.

dmap_path_to_bismark_exes="";

# bismark index building parameters: recent versions of bismark default to bowtie2.  bowtie1 is
# no longer an option.  Let's spell it out anyhow

bismark_build_options="--bowtie2";

# bismark run parameters: recent versions of bismark default to bowtie2.  bowtie1 is
# no longer an option.  Let's spell it out anyhow.  The -N 1 option is what we
# have routinely used to minimise uncertainty in the 'seed' part of the alignment.
# It also performs significantly faster than the default -n 2 in our trials.

bismark_run_options-"-N 1 --bowtie2 ";

;;

  bsmapz)

# set bsmap stuff

# path to run bsmapz

bsmapz_exe="bsmapz"

# number of processor cores for bsmapz run

bsmapz_cores=4;

# Mismatches for bsmapz: if this value is between 0 and 1, it's interpreted as the mismatch rate w.r.t to the read length.
# otherwise it's interpreted as the maximum number of mismatches allowed on a read, <=15.
# example: -v 5 (max #mismatches = 5), -v 0.1 (max #mismatches = read_length * 10%): Note, defaults to 0.08

bsmapz_mismatch_val=2;

bsmapz_run_options=" -p ${bsmapz_cores} -v ${bsmapz_mismatch_val}";

;;

  *)

 printf "Error: invalid bisulphite_mapper value '%s' - expected 'bismark' or 'bsmapz'\n" "${bisulphite_mapper}";
 exit 1;

;;

esac

# name of chromosome information file: this will be created during the dmap_index_build.sh
# step

dmap_chr_info_file="dmap_chr_info.txt";

# Annotation information for identgeneloc operation

# Type of annotation data: one of
# EMBL
# Genbank
# SeqMonk
# GFF3
# GTF
# none

#feature_annotation_type="GTF";

# Location of feature annotation files:

#annotation_file_location="/Volumes/dsm-pathology-ecclesRNA/Erin_Macaulay_placental/hs_gencode_GRCh37/";

# An array of annotation files.  GTF & GFF3 only need one entry, others
# need one file/chromosome

#annotation_files=(
#"gencode.v32lift37.annotation.gtf"
#);

dmap_annot_info_file="dmap_annot_info.txt";

## Example of parameters for SeqMonk feature information
#
feature_annotation_type="SeqMonk";
#
annotation_file_location="/Users/stope71p/seqmonk_genomes/Homo sapiens/GRCh37"/
#
annotation_files=(
"1.dat"
"10.dat"
"11.dat"
"12.dat"
"13.dat"
"14.dat"
"15.dat"
"16.dat"
"17.dat"
"18.dat"
"19.dat"
"2.dat"
"20.dat"
"21.dat"
"22.dat"
"3.dat"
"4.dat"
"5.dat"
"6.dat"
"7.dat"
"8.dat"
"9.dat"
"MT.dat"
"X.dat"
"Y.dat"
)

