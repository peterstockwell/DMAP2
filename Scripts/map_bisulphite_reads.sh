# map_bisulphite_reads.sh: script to adaptor trim one bisulphite
# sample and run the bismark aligner on it.
#
# Peter Stockwell: 6-Apr-2022


# needs a basic parameter file to specify positions of relevant
# executables and genome files. Also needs a file specific to
# each individual mapping run giving raw fastq file name(s) and
# any parameters for trimming

# usage: ./map_bisulphite_reads.sh <dmap_basic_params.sh> <sample_parameters>

# check for parameter files:

if [[ -z $1 || -z $2 ]]; then

echo "This script needs two parameter files:";
echo " 1. Name of basic parameter file";
printf " 2. Name of sample parameter file\n";

exit 1;

fi

# can we read this file?

if [[ ! -f "$1" ]]; then

  printf "Can't read basic parameter file '%s'\n" "$1";

  exit 1;

else

# pick up definitions from parameter file

. "$1";

fi

# check for sample parameter file ($2)

if [[ ! -f "$2" ]]; then

  printf "Can't read sample parameter file '%s'\n" "$2";

  exit 1;

else

# pick up sample definitions

. "$2";

fi

if [[ $verbose == "yes" ]]; then

printf "                DMAP\n";
printf " Differential Methylation Analysis Package\n";
printf "Adaptor trimmming and Mapping sequence reads\n\n";
printf "  Basic parameters from '%s'\n" "$1";
printf "  Sample parameters from '%s'\n" "$2";

fi

# create the contaminant file for adaptor trimming, if it doesn't exist already:

if [[ ! -f "contam.fa" ]]; then

cat << CONTAM.FA > contam.fa
>contam|SingleEndAdapter1
ACACTCTTTCCCTACACGACGCTGTTCCATCT
>contam|SingleEndAdapter2
CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT
>contam|SingleEndPCRPrimer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|SingleEndPCRPrimer2
CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT
>contam|SingleEndSequencingPrimer
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|PairedEndAdapter1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|PairedEndAdapter2
CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>contam|PariedEndPCRPrimer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|PairedEndPCRPrimer2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>contam|PariedEndSequencingPrimer1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|PairedEndSequencingPrimer2
CGGTCTCGGCATTCCTACTGAACCGCTCTTCCGATCT
>contam|DpnIIexpAdapter1
ACAGGTTCAGAGTTCTACAGTCCGAC
>contam|DpnIIexpAdapter2
CAAGCAGAAGACGGCATACGA
>contam|DpnIIexpPCRPrimer1
CAAGCAGAAGACGGCATACGA
>contam|DpnIIexpPCRPrimer2
AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA
>contam|DpnIIexpSequencingPrimer
CGACAGGTTCAGAGTTCTACAGTCCGACGATC
>contam|NlaIIIexpAdapter1
ACAGGTTCAGAGTTCTACAGTCCGACATG
>contam|NlaIIIexpAdapter2
CAAGCAGAAGACGGCATACGA
>contam|NlaIIIexpPCRPrimer1
CAAGCAGAAGACGGCATACGA
>contam|NlaIIIexpPCRPrimer2
AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA
>contam|NlaIIIexpSequencingPrimer
CCGACAGGTTCAGAGTTCTACAGTCCGACATG
>contam|SmallRNAAdapter1
GTTCAGAGTTCTACAGTCCGACGATC
>contam|SmallRNAAdapter2
TCGTATGCCGTCTTCTGCTTGT
>contam|SmallRNARTPrimer
CAAGCAGAAGACGGCATACGA
>contam|SmallRNAPCRPrimer1
CAAGCAGAAGACGGCATACGA
>contam|SmallRNAPCRPrimer2
AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA
>contam|SmallRNASequencingPrimer
CGACAGGTTCAGAGTTCTACAGTCCGACGATC
>contam|MultiplexingAdapter1
GATCGGAAGAGCACACGTCT
>contam|MultiplexingAdapter2
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|MultiplexingPCRPrimer1.01
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|MultiplexingPCRPrimer2.01
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>contam|MultiplexingRead1SequencingPrimer
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|MultiplexingIndexSequencingPrimer
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>contam|MultiplexingRead2SequencingPrimer
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>contam|PCRPrimerIndex1
CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC
>contam|PCRPrimerIndex2
CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC
>contam|PCRPrimerIndex3
CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC
>contam|PCRPrimerIndex4
CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC
>contam|PCRPrimerIndex5
CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC
>contam|PCRPrimerIndex6
CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC
>contam|PCRPrimerIndex7
CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC
>contam|PCRPrimerIndex8
CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC
>contam|PCRPrimerIndex9
CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC
>contam|PCRPrimerIndex10
CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC
>contam|PCRPrimerIndex11
CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC
>contam|PCRPrimerIndex12
CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC
>contam|DpnIIGexAdapter1
GATCGTCGGACTGTAGAACTCTGAAC
>contam|DpnIIGexAdapter1.01
ACAGGTTCAGAGTTCTACAGTCCGAC
>contam|DpnIIGexAdapter2
CAAGCAGAAGACGGCATACGA
>contam|DpnIIGexAdapter2.01
TCGTATGCCGTCTTCTGCTTG
>contam|DpnIIGexPCRPrimer1
CAAGCAGAAGACGGCATACGA
>contam|DpnIIGexPCRPrimer2
AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA
>contam|DpnIIGexSequencingPrimer
CGACAGGTTCAGAGTTCTACAGTCCGACGATC
>contam|NlaIIIGexAdapter1.01
TCGGACTGTAGAACTCTGAAC
>contam|NlaIIIGexAdapter1.02
ACAGGTTCAGAGTTCTACAGTCCGACATG
>contam|NlaIIIGexAdapter2.01
CAAGCAGAAGACGGCATACGA
>contam|NlaIIIGexAdapter2.02
TCGTATGCCGTCTTCTGCTTG
>contam|NlaIIIGexPCRPrimer1
CAAGCAGAAGACGGCATACGA
>contam|NlaIIIGexPCRPrimer2
AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA
>contam|NlaIIIGexSequencingPrimer
CCGACAGGTTCAGAGTTCTACAGTCCGACATG
>contam|SmallRNA_RT_Primer
CAAGCAGAAGACGGCATACGA
>contam|5pRNAAdapter
GTTCAGAGTTCTACAGTCCGACGATC
>contam|RNAAdapter1
TCGTATGCCGTCTTCTGCTTGT
>contam|SmallRNA3pAdapter1
ATCTCGTATGCCGTCTTCTGCTTG
>contam|SmallRNAPCRPrimer1
CAAGCAGAAGACGGCATACGA
>contam|SmallRNAPCRPrimer2
AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA
>contam|SmallRNASequencingPrimer
CGACAGGTTCAGAGTTCTACAGTCCGACGATC
>contam|TruSeqUniversalAdapter
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>contam|TruSeq_AdapterIndex1
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex2
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex3
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex4
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex5
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex6
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex7
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex8
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex9
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex10
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex11
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex12
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex13
GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex14
GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex15
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex16
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex18
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex19
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex20
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex21
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex22
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex23
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex25
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG
>contam|TruSeq_AdapterIndex27
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG
>contam|RNA_RTPrimer
GCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCRPrimer
AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA
>contam|RNA_PCR_PrimerIndex1
CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex2
CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex3
CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex4
CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex5
CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex6
CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex7
CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex8
CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex9
CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex10
CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex11
CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex12
CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex13
CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex14
CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex15
CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex16
CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex17
CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex18
CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex19
CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex20
CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex21
CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex22
CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex23
CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex24
CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex25
CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex26
CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex27
CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex28
CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex29
CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex30
CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex31
CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex32
CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex33
CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex34
CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex35
CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex36
CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex37
CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex38
CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex39
CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex40
CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex41
CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex42
CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex43
CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex44
CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex45
CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex46
CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex47
CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|RNA_PCR_PrimerIndex48
CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
>contam|Nextera_transposase_R1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>contam|Nextera_transposase_R2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>NxtTransposonEndSequence
AGATGTGTATAAGAGACAG
>NxtPrimer_1*
AATGATACGGCGACCACCGA
>NxtAdaptor_1
AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAG
>NxtPrimer_2
CAAGCAGAAGACGGCATACGA
>NxtAdaptor_2 (minus bar code)*:
CAAGCAGAAGACGGCATACGAGATCGGTCTGCCTTGCCAGCCCGCTCAG
>Nextera_R1_Primer
GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG
>Nextera_R2_Primer
GCCTTGCCAGCCCGCTCAGAGATGTGTATAAGAGACAG
>NexteraIndexRead_Primer
CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGCAGACCG
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
CONTAM.FA

fi

# adaptor trimming with cleanadaptors:  Set hard trim if required

if [[ dmap_sample_trim_length -gt 0 ]]; then

sample_trim="-l ""${dmap_sample_trim_length}"" ";

else

sample_trim="";

fi

if [[ $dmap_run_type == "RRBS" ]]; then

trim_back="-t 3 ";

else

trim_back="";

fi

if [[ dmap_min_read_length -gt 0 ]]; then

min_length="-x ""${dmap_min_read_length}"" ";

else

min_length="";

fi

# is read 1 file compressed?

printf "R1 sample is: '%s'\n" "${dmap_sample_files[0]}";

if [[ "${dmap_sample_files[0]}" == *fastq.gz ]]; then

comp_option=("-z ");
read_trailer=(".fastq.gz");

else

comp_option=("-Z ");
read_trailer=(".fastq");

fi

# generate trimmed name for read1

read_out_name=(
"${adtrimmed_out_dir}""$(basename "${dmap_sample_files[0]}" "${read_trailer[0]}")""_at.fastq"
)

printf "Trimming R1 to '%s'\n" "${read_out_name[0]}";

# do we have paired-end data??

if [[ ${#dmap_sample_files[@]} -gt 1 ]]; then

if [[ -f "${dmap_sample_files[1]}" ]]; then

# paired end trim:

# is read 2 file compressed?

if [[ "${dmap_sample_files[1]}" == *fastq.gz ]]; then

comp_option+=("-z ");
read_trailer+=(".fastq.gz");

else

comp_option+=("-Z ");
read_trailer+=(".fastq");

fi

read_out_name+=("${adtrimmed_out_dir}""$(basename "${dmap_sample_files[1]}" "${read_trailer[1]}")""_at.fastq");

if [[ $verbose == "yes" ]]; then

printf "Paired end adaptor trimming to '%s' and '%s'\n" "${read_out_name[0]}" "${read_out_name[1]}";

fi

trim_command="${path_to_dmap}""cleanadaptors -I contam.fa ""${sample_trim}""${trim_back}""${min_length}""${comp_option[0]}""-F ""${dmap_sample_files[0]}"" ${comp_option[1]}""-G ""${dmap_sample_files[1]}"" -Z -o ""${read_out_name[0]}"" -O ""${read_out_name[1]}";

bismark_options=(
"-1 ""${read_out_name[0]}"
"-2 ""${read_out_name[1]}");

else

printf "Can't open read 2 file '%s'\n" "${dmap_sample_files[1]}";

exit 1;

fi

else

# single ended trim

if [[ $verbose == "yes" ]]; then

printf "Single ended adaptor trimming to '%s'\n" "${read_out_name[0]}";

fi

trim_command="${path_to_dmap}""cleanadaptors -I contam.fa ""${sample_trim}""${trim_back}""${min_length}""${comp_option[0]}""-F ""${dmap_sample_files[0]}"" -Z -o ""${read_out_name[0]}";

bismark_options=("${read_out_name[0]}" "");

fi

if [[ $verbose == "yes" ]]; then

printf "Executing: %s\n" "${trim_command}";

fi

eval "${trim_command}";

# we need to run bismark on the adaptor trimmed file(s)

# check the bismark output dir exists

mkdir -p "${bismark_out_dir}";

# create 'unique' name for bismark command file:

bismark_cmd_file="${adtrimmed_out_dir}""$(basename "${dmap_sample_files[0]}" "${read_trailer[0]}")""_bismark_cmd.sh";

printf "%sbismark %s -o %s %s %s %s\n" "${dmap_path_to_bismark_exes}" "${dmap_bismark_index_location}" "${bismark_out_dir}" "${bismark_run_options}" "${bismark_options[0]}" "${bismark_options[1]}" > "${bismark_cmd_file}";

if [[ $verbose == "yes" ]]; then

    printf "Executing: %sbismark %s -o %s %s %s %s\n" "${dmap_path_to_bismark_exes}" "${dmap_bismark_index_location}" "${bismark_out_dir}" "${bismark_run_options}" "${bismark_options[0]}" "${bismark_options[1]}";

fi

/bin/bash "${bismark_cmd_file}"

# if verbose, then run check on mapping stats:

if [[ $verbose == "yes" ]]; then

# make the awk script

cat << 'BISMARKSTATS.AWK'  > bismark_stats.awk
# bismark_stats.awk: to retrieve information from bismark mapping reports
#

BEGIN{printf("File(s)\tTotalReads\tUniqReads\tMapping\tCpGmeth\tCHGmeth\tCHHmeth\n");}

$0~/Bismark report for:/&&$0!~/ and /{printf("%s\t",$4);}

$0~/Bismark report for:/&&$0~/ and /{printf("%s&%s\t",$4,$6);}

$0~/Sequences analysed in total:/||$0~/Sequence pairs analysed in total:/{printf("%s\t",$NF);}

$0~/Number of alignments with a unique/||$0~/Number of paired-end alignments with a unique/{printf("%s\t",$NF);}

$0~/Mapping efficiency:/{printf("%s\t",$NF);}

$0~/C methylated in CpG context:/{printf("%s\t",$NF);}

$0~/C methylated in CHG context:/{printf("%s\t",$NF);}

$0~/C methylated in CHH context:/{printf("%s\n",$NF);}
BISMARKSTATS.AWK

if [[ ${#dmap_sample_files[@]} -gt 1 ]]; then

report_file_name="${adtrimmed_out_dir}""$(basename "${dmap_sample_files[0]}" "${read_trailer[0]}")""_at_bismark_bt2_PE_report.txt";

else

report_file_name="${adtrimmed_out_dir}""$(basename "${dmap_sample_files[0]}" "${read_trailer[0]}")""_at_bismark_bt2_SE_report.txt";

fi

# generate the report

awk -f bismark_stats.awk "${report_file_name}"

fi

exit 0;
