# dmap_index_build.sh: file to generate bismark indices for
# mapping, plus other information files required for DMAP running.
# Also generate annotation file information.
# Peter Stockwell: Aug-2021, Oct-2021


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
printf "      Building bismark index\n";

# generate chromosome information file from of genome sequence files
# start by deleting any existing file

if [[ $verbose == "yes" ]]; then

printf "Writing chromosome information to '%s'\n" "${dmap_chr_info_file}";

fi

rm -f "${dmap_chr_info_file}";

for genomefile in "${dmap_genome_fasta_files[@]}"; do

genome_path_and_file="${dmap_genome_file_location}""${genomefile}";

if [[ -f "${genome_path_and_file}" ]]; then

head -1 "${genome_path_and_file}" | awk '{printf("\"%s\"\t\"%s\"\n",substr($1,2),fullname);}' fullname="${genome_path_and_file}" >> "${dmap_chr_info_file}";

else

printf "Can't open genome file: '%s'\n" "${genome_path_and_file}";

exit 1;

fi

done

if [[ $verbose == "yes" ]]; then

printf "Executing bismark genome indexing: this may take appreciable time\n";

fi

printf "%sbismark_genome_preparation %s %s\n" "${dmap_path_to_bismark_exes}" "${bismark_build_options}" "${dmap_bismark_index_location}" | /bin/bash

# printf "%sbismark_genome_preparation %s %s\n" "${dmap_path_to_bismark_exes}" "${bismark_build_options}" "${dmap_bismark_index_location}" > bismark_index_build.sh

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
        head -1 "${annot_path_and_file}" | awk '{printf("\"%s\"\t\"%s\"\n",$2,fullname);}' fullname="${annot_path_and_file}" >> "${dmap_annot_info_file}";
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
        head -1 "${annot_path_and_file}" | awk '{split($2,s2,":");printf("\"%s\"\t\"%s\"\n",s2[2],fullname);}' fullname="${annot_path_and_file}" >> "${dmap_annot_info_file}";
      else
        printf "Can't open annotation file '%s'\n" "${annot_path_and_file}";
	exit 1;
      fi
    done
;;
 GFF3 | GTF)
# don't need to do anything for these
;;
 *)
   printf "%s feature annotation type not yet implemented\n" ${feature_annotation_type};
;;

esac

exit 0;
