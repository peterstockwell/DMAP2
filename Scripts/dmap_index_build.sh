# dmap_index_build.sh: file to generate bismark indices for bismark
# mapping, or appropriate genome files for bsmapz mapping.
# Plus other information files required for DMAP running.
# Also generate annotation file information.
# Peter Stockwell: Aug-2022


# needs a basic parameter file to specify positions of relevant
# executables and genome files.

# usage: ./dmap_index_build <dmap_basic_params.sh>

# check for parameter file:

if [[ -z $1 ]]; then

echo "This script needs a parameter:";
echo " 1. Name of basic parameter file";

exit 1;

fi

# can we read this file?

if [[ ! -f "$1" ]]; then

  printf "Can't read parameter file '%s'\n" "$1";

  exit 1;

else

# pick up definitions from parameter file

. "$1";

fi

printf "                DMAP\n";
printf " Differential Methylation Analysis Package\n";
printf "\nWorking in dir '%s'\nBasic parameters from '%s'\n\n" "$PWD" "$1";

case $bisulphite_mapper in

  bismark)
    printf "      Building bismark index\n";
    ;;
  bsmapz)
    printf "      Prepping genome for bsmapz\n";
    ;;
esac

if [[ $verbose == "yes" ]]; then

printf "Writing chromosome information to '%s'\n" "${dmap_chr_info_file}";

fi

# generate chromosome information file from of genome sequence files
# start by deleting any existing file

rm -f "${dmap_chr_info_file}";

# are genome file(s) compressed?

genome_file_ext="${dmap_genome_fasta_files[0]##*.}"

if [[ $genome_file_ext == "gz" ]]; then

comp_ext=".gz"

uncomp_file_name="${dmap_genome_fasta_files[0]%.*}";

else

comp_ext="";

uncomp_file_name="${dmap_genome_fasta_files[0]}";

fi

# find genome file extension ('.fa', '.fasta', '.fa.gz' or '.fasta.gz')

genome_file_ext="${uncomp_file_name##*.}"

if [[ $genome_file_ext != "fa" && $genome_file_ext != "fasta" ]]; then
  printf "Genome fasta file extension '.%s' invalid, should be '.fa' or '.fasta'\n" "${genome_file_ext}";
  exit 1;
fi

# add the leading '.'

genome_file_ext=".""${genome_file_ext}";

# is genome in multiple files or single file?
# DMAP:diffmeth requires each chromosome in separate file,
# bismark is happy either way but bsmapz needs a single
# whole genome fasta file, so create both forms if
# single file or if multiple and bsmapz is mapper.

if [[ ${#dmap_genome_fasta_files[@]} -gt 1 ]]; then

# create name for single genome fasta file

  base_genome_name="$(basename "${dmap_genome_fasta_files[0]}")";
  whole_genome_filepath="${dmap_genome_dir}""genome_single/""${base_genome_name%%.*}"".wholegenome""${genome_file_ext}";

# create this dir if necessary:

  mkdir -p "${dmap_genome_dir}""genome_single/"

  if [[ $bisulphite_mapper == "bismark" ]]; then

    if [[ $verbose == "yes" ]]; then

    printf "Executing bismark genome indexing: this may take some hours\n";

    fi

  printf "%sbismark_genome_preparation %s %s\n" "${dmap_path_to_bismark_exes}" "${bismark_build_options}" "${dmap_bismark_index_location}" | /bin/bash

#printf "%sbismark_genome_preparation %s %s\n" "${dmap_path_to_bismark_exes}" "${bismark_build_options}" "${dmap_bismark_index_location}" > bismark_index_build.sh

  fi

  if [[ $verbose == "yes" && $bisulphite_mapper == "bsmapz" ]]; then

  printf "Combining chromosome fasta files into '%s'\n" "${whole_genome_filepath}";

# make sure we don't already have such a file

    if [[ -f "${whole_genome_filepath}" ]]; then

      printf "Deleting existing whole genome file '%s' and recreating it\n" "${whole_genome_path}";

    fi

  rm -f "${whole_genome_filepath}";

  fi


for genomefile in "${dmap_genome_fasta_files[@]}"; do

# compressed?? if so need to generate uncompressed versions:

genome_path_and_file="${dmap_genome_dir}""$(basename "${genomefile}")";

# need to append "/multi_genome/" to path for multi destination

dest_path_and_file="$(dirname "${genome_path_and_file}")""/multi_genome/""$(basename "${genome_path_and_file}")";

# check the destination dir exists:

mkdir -p "$(dirname "${genome_path_and_file}")""/multi_genome/";

  if [[ -f "${genome_path_and_file}" ]]; then

    if [[ $comp_ext != "" ]]; then

      gzip -dc "${genome_path_and_file}" > "${dest_path_and_file%.*}";

      gzip -dc "${genome_path_and_file}" | head -1 | awk '{printf("\"%s\"\t\"%s\"\n",substr($1,2),fullname);}' fullname="${dest_path_and_file%.*}" >> "${dmap_chr_info_file}";

      if [[ $bisulphite_mapper == "bsmapz" ]]; then

        gzip -dc "${genome_path_and_file}" >> "${whole_genome_filepath}";

      fi

    else

      head -1 "${genome_path_and_file}" | awk '{printf("\"%s\"\t\"%s\"\n",substr($1,2),fullname);}' fullname="${genome_path_and_file}" >> "${dmap_chr_info_file}";

      if [[ $bisulphite_mapper == "bsmapz" ]]; then

        cat "${genome_path_and_file}" >> "${whole_genome_filepath}";

      fi

    fi

  else

    printf "Can't open genome file: '%s'\n" "${genome_path_and_file}";

    exit 1;

  fi

done

else   # single genome fasta file

whole_genome_filepath="${dmap_genome_dir}""$(basename "${dmap_genome_fasta_files[0]}")";

printf "Single genome file: '%s'\n" "${whole_genome_filepath}";

# check we can open this

if [[ -f "${whole_genome_filepath}" ]]; then

  if [[ $bisulphite_mapper == "bismark" ]]; then

    if [[ $verbose == "yes" ]]; then

    printf "Executing bismark genome indexing: this may take some hours\n";

    fi

  printf "%sbismark_genome_preparation %s %s\n" "${dmap_path_to_bismark_exes}" "${bismark_build_options}" "${dmap_bismark_index_location}" | /bin/bash

#  printf "%sbismark_genome_preparation %s %s\n" "${dmap_path_to_bismark_exes}" "${bismark_build_options}" "${dmap_bismark_index_location}" > bismark_index_build.sh

  fi

# want to split the files now, after bismark genome prep (if required).

# Single fasta genome file: start by making awk script to split it

cat << 'SPLIT_FASTA.AWK' > split_fasta.awk
# split_fasta.awk: script to break multi-entry fasta
# files into separate files, each containing 1 entry.
# files named by the fasta IDs, with optional header
#
# Peter Stockwell: 21-Feb-2022
#
# usage:
# awk -f split_fasta.awk <multi_fasta_file>
#
# or piped, e.g.
#
# cat <multi_fasta_file> | awk -f split_fasta.awk
#
# Can define optional useful parameters in command line:
# namehdr - prefixed to each output file name, can be a directory or string or both (e.g. namehdr="singles_fasta/")
# extension - use as file extension, default is ".fa"
# infofile - will write table of Chr IDs and full file paths to this file
#
# use as: e.g.
#
# awk -f split_fasta.awk namehdr="genome_singles/" extension=".fasta" infofile="chr_info.txt" <multi_fasta_file>
#
# to produce chr_info.txt which is suitable for the diffmeth -G option
# If namehdr is a directory, it must already exist.
#

BEGIN {namehdr="";
extension=".fa";
outfile="";
infofile=""
}

$0~/>/{if (outfile != "")
  close(outfile);
outfile = sprintf("%s.%s%s",namehdr,substr($1,2),extension);
printf("%s\n",$0) > outfile;
if (infofile != "")
  printf("\"%s\"\t\"%s\"\n",substr($1,2),outfile) > infofile;
}

$0!~/>/{
printf("%s\n",$0) > outfile;
}

END{
if (outfile != "")
  close(outfile);
if (infofile != "")
  close(infofile);
}
SPLIT_FASTA.AWK

  genome_base_file="$(basename "${dmap_genome_fasta_files[0]}")";
  genome_base_head="${genome_base_file%%.*}";

  genome_multi_dir="${dmap_genome_dir}""multi_genome/";
  whole_genome_path_hdr="${genome_multi_dir}""/""${genome_base_head}";

  if [[ $verbose == "yes" ]]; then

  printf "Splitting whole genome fasta file to separate chromosomes in '%s'\n" "${genome_multi_dir}";

  fi

# create the multi_genome_dir if necessary

  mkdir "${genome_multi_dir}";

# if genome files are compressed, need to uncompress them

if [[ $comp_ext != "" ]]; then

  gzip -dc "${whole_genome_filepath}" | awk -f split_fasta.awk namehdr="${whole_genome_path_hdr}" infofile="${dmap_chr_info_file}" extension="${genome_file_ext}"

else

  awk -f split_fasta.awk namehdr="${whole_genome_path_hdr}" infofile="${dmap_chr_info_file}" extension="${genome_file_ext}" "${whole_genome_filepath}";

fi

  else

  printf "Can't open genome file '%s'\n" "${whole_genome_file}";
  exit 1;

  fi

fi

# Generate feature annotation information if required

case $feature_annotation_type in

  Genbank | EMBL)
    if [[ $verbose == "yes" ]]; then
      printf "Writing annotation information to '%s'\n" "${dmap_annot_info_file}";
    fi
# remove any existing annotation info file
    rm -f "${dmap_annot_info_file}";

    for annot_file in "${annotation_files[@]}"; do
      annot_path_and_file="${annotation_file_location}""${annot_file}";
      if [[ -f "${annot_path_and_file}" ]]; then
# check for .gz
        if [[ "${annot_path_and_file##*.}" == "gz" ]]; then
          if [[ $verbose == "yes" ]]; then
	    printf "Writing decompressed file to %s\n" "${annot_path_and_file%.*}";
	  fi
	  gzip -dc "${annot_path_and_file}" > "${annot_path_and_file%.*}";
          head -1 "${annot_path_and_file%.*}" | awk '{printf("\"%s\"\t\"%s\"\n",$2,fullname);}' fullname="${annot_path_and_file%.*}" >> "${dmap_annot_info_file}";
	else
          head -1 "${annot_path_and_file}" | awk '{printf("\"%s\"\t\"%s\"\n",$2,fullname);}' fullname="${annot_path_and_file}" >> "${dmap_annot_info_file}";
	fi
      else
        printf "Can't open annotation file '%s'\n" "${annot_path_and_file}";
	exit 1;
      fi
    done
;;
  SeqMonk)
    if [[ $verbose == "yes" ]]; then
      printf "Writing annotation information to '%s'\n" "${dmap_annot_info_file}";
    fi
# remove any existing annotation info file
    rm -f "${dmap_annot_info_file}";

    for annot_file in "${annotation_files[@]}"; do
      annot_path_and_file="${annotation_file_location}""${annot_file}";
      if [[ -f "${annot_path_and_file}" ]]; then
# make awk script to do this

cat << 'GET_SEQMONK_INFO' > get_seqmonk_info.awk
{split($2,s2,":");printf("\"%s\"\t\"%s\"\n",s2[2],fullname);}
GET_SEQMONK_INFO

        head -1 "${annot_path_and_file}" | awk -f get_seqmonk_info.awk fullname="${annot_path_and_file}" >> "${dmap_annot_info_file}";
      else
        printf "Can't open annotation file '%s'\n" "${annot_path_and_file}";
	exit 1;
      fi
    done
;;
 GFF3 | GTF)
# don't need to do anything for these, but check for .gz
  if [[ "${annotation_files[0]##*.}" == "gz" ]]; then
    annotation_path_and_file="${annotation_file_location}""${annotation_files[0]}";
    if [[ $verbose == "yes" ]]; then
      printf "Writing decompressed file to %s\n" "${annotation_path_and_file%.*}";
    fi
    gzip -dc "${annotation_path_and_file}" > "${annotation_path_and_file%.*}";
  fi
;;

 none)
# do nothing
;;

*)
   printf "%s feature annotation type not yet implemented\n" ${feature_annotation_type};
;;

esac

if [[ $verbose == "yes" ]]; then

printf "Genome preparation complete\n";

fi

exit 0;
