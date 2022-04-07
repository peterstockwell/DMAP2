# dmap_run_diffmeth.sh: script to run the differential
# methylation application diffmeth in various modes
# on previously-mapped sample reads.
#
# Peter Stockwell: Aug/Sep/Oct/Nov/Dec-2021,Jan-2022


# needs a basic parameter file to specify positions of relevant
# executables and genome files. Also needs a parameter file specific to
# each diffmeth run to specify run type and mapping files and the
# criteria for accepting/rejecting fragments or regions.

# usage: ./dmap_run_diffmeth.sh <dmap_basic_params.sh> <diffmeth_run_parameters>

# check for parameter files:

if [[ -z $1 || -z $2 ]]; then

echo "This script needs two parameter files:";
echo " 1. Name of basic parameter file";
printf " 2. Name of diffmeth parameter file\n";

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

# check for diffmeth parameter file ($2)

if [[ ! -f "$2" ]]; then

  printf "Can't read diffmeth parameter file '%s'\n" "$2";

  exit 1;

else

# pick up sample definitions

. "$2";

fi

if [[ $verbose == "yes" ]]; then

printf "                    DMAP\n";
printf "       Differential Methylation Analysis Package\n";
printf "Running differential methylation analysis on mapped reads\n\n";
printf "  Basic parameters from '%s'\n" "$1";
printf "  Sample parameters from '%s'\n" "$2";

fi

if [[ $verbose == "yes" ]]; then

   printf "Run mode is %s with %s data \n" "$diffmeth_run_type" "$dmap_run_type";

fi

# check for valid run type and appropriate file definitions

# define useful function(s)

function append_sam_or_bam () {

# check for sam or bam extension and that the file is readable

if [[ -f "$2" ]]; then
  if [[ $2 == *.bam ]]; then
    dmeth_cmd="${dmeth_cmd}"" -z ""$1 $2";
  else
    dmeth_cmd="${dmeth_cmd}"" -Z ""$1 $2";
  fi
else
  printf "Can't read file '%s'\n" $2;
  exit 1;
fi
}
  
# set some basic parameters for diffmeth

diffmeth_run_options="";

if [[ $dmap_run_type == "RRBS" ]]; then

  if [[ $map_init_CpG_to_prev == "yes" ]]; then

    diffmeth_run_options=" -N ";

  else

    diffmeth_run_options="";

  fi

else

  diffmeth_run_options=" -W ""$wgbs_window";

fi

if [[ $min_valid_cpgs -gt 0 ]]; then

diffmeth_run_options="${diffmeth_run_options}"" -F ""${min_valid_cpgs}";

fi

if [[ $min_cpg_hits -gt 0 ]]; then

diffmeth_run_options="${diffmeth_run_options}"" -t ""${min_cpg_hits}";

fi

# nonCpG methylaton??

if [[ $non_CpG_methylation == "yes" ]]; then

diffmeth_run_options=${diffmeth_run_options}" -J ";

fi

# min_sample_count is not wanted for pairwise or listing options

case $diffmeth_run_type in

  Anova | chisquare)
    if [[ $min_sample_count -gt 0 ]]; then
      diffmeth_run_options="${diffmeth_run_options}"" -I ""${min_sample_count}";
    fi

;;

  *)
;;

esac

if [[ -n $pr_threshold ]]; then

diffmeth_run_options="${diffmeth_run_options}"" -U ""${pr_threshold}";

fi

if [[ $max_sam_read_length -gt 0 ]]; then

diffmeth_run_options=-"${diffmeth_run_options}"" -b ""${max_sam_read_length}";

fi

case $diffmeth_run_type in

 pairwise)
   if (( ${#file_list1[@]} < 2 )); then
     printf "file_list1 does not contain 2 sam/bam file names\n";
     exit 1;
   else
     dmeth_cmd="${path_to_dmap}""diffmeth -G ""${dmap_chr_info_file}""${diffmeth_run_options}"" -P ""$diffmeth_fragment_min"",""$diffmeth_fragment_max";
     if [[ $verbose == "yes" ]]; then
       printf " with files";
     fi
     fileoption="-R";
     for file in "${file_list1[@]}"; do
       if [[ $verbose == "yes" ]]; then
         printf " '%s'" "$file";
       fi
       append_sam_or_bam "$fileoption" "$file";
       fileoption="-S";
     done

    dmeth_cmd="${dmeth_cmd}"" > ${diffmeth_out_name}";

  fi
;;

 chisquare)
# check for at least 2 input files:

  if (( ${#file_list1[@]} < 2 )); then
    printf "ChiSquare requires at least 2 file names in file_list1 array\n";
    exit 1;
  else
    dmeth_cmd="${path_to_dmap}""diffmeth -G ""${dmap_chr_info_file}""${diffmeth_run_options}"" -X ""$diffmeth_fragment_min"",""$diffmeth_fragment_max";
      if [[ $verbose == "yes" ]]; then
        printf "with files:";
      fi
      for file in "${file_list1[@]}"
        do
        if [[ $verbose == "yes" ]]; then
          printf " '%s'" "$file";
        fi
        append_sam_or_bam "-R" "$file";
        done
      if [[ $verbose == "yes" ]]; then
        printf "\n";
      fi

    dmeth_cmd="${dmeth_cmd}"" > ${diffmeth_out_name}";
  fi

;;

  Anova)
# check for two or more files in each list
    for ((agroup=1; agroup<=anova_group_number; agroup++)); do
      newname="file_list""$agroup";
      namarray=$(eval echo \${${newname}[@]});
      IFS=', ' read -a namarray <<< $(eval echo \${${newname}[@]})
      if (( ${#namarray[@]}<2 )); then
        printf "Need 2 or more bam/sam files in %s\n" "${newname}";
        exit 1
        fi
      done
      dmeth_cmd="${path_to_dmap}""diffmeth -G ""${dmap_chr_info_file}""${diffmeth_run_options}";
      case $anova_detail in
        anova_simple)
	  dmeth_cmd="${dmeth_cmd}"" -a ";
          anova_fold_difference="";
	  ;;
        anova_medium)
	  dmeth_cmd="${dmeth_cmd}"" -A ";
	  ;;
	anova_full)
	  dmeth_cmd="${dmeth_cmd}"" -B ";
	  ;;
      esac
# do we want fold difference?
    if [[ $anova_fold_difference == "yes" ]]; then

cat << 'METHPROPS2FOLD' > methprops2fold.awk
# methprops2fold.awk: to take diffmeth output with the final column
# headed by R=rprop,S=sprop (e.g. R=0.1250,S=0.2093;Ra=0.1000,Re=0.1429;Sc=0.0833,Se=0.2581)
# and return fold difference.
#
# Note: 0.0 values will return -99.0
#
# Usage:
#
# awk -f methprops2fold.awk <diffmeth_output_file>
#
# or pipe diffmeth output directly to 'awk -f methprops2fold.awk'
#
# Peter Stockwell: 17-Feb-2022
#

$0~/#/{printf("%s\tFoldDiff\n",$0);}

$0!~/#/{split($NF,snf,";");
  fdiff = split_fold_diff(snf[1]);
  if (fdiff < 0.0)
    printf("%s\t-\n",$0);
  else
    printf("%s\t%.4f\n",$0,fdiff);
  }

function fold_diff(x,y)
# return fold difference between x & y (real values), -99.0 for
# invalid.
{
if ((x == 0.0) || (y == 0.0))
  return(-99.0);
else
  if (x > y)
    return(x/y);
  else
    return(y/x);
}

function split_fold_diff(propstring)
# split propstring into R & S real values and
# return fold difference
{
ns = split(propstring,pstr,/[=,]/);
return(fold_diff(pstr[2]+0.0,pstr[ns]+0.0));
}
METHPROPS2FOLD
      fi
    dmeth_cmd="${dmeth_cmd}""$diffmeth_fragment_min"",""$diffmeth_fragment_max";
    for ((agroup=1; agroup<=anova_group_number; agroup++)); do
      if [[ $verbose == "yes" ]]; then
        printf "with group%d files:\n" ${agroup};
        fi
      newname="file_list""$agroup";
      IFS=', ' read -a namarray <<< $(eval echo \${${newname}[@]})
      for file in "${namarray[@]}"; do
        if [[ $verbose == "yes" ]]; then
          printf " '%s'\n" "$file";
          fi
        if (( ${anova_group_number}==2 )); then
          if (( ${agroup}==1 )); then
	    fileoption="-R";
	  else
	    fileoption="-S";
	  fi
	else
	  fileoption="-R""${agroup}";
	fi
        append_sam_or_bam "${fileoption}" "$file";
        done
      done

    if [[ $anova_fold_difference == "yes" ]]; then
        dmeth_cmd="$dmeth_cmd"" | awk -f methprops2fold.awk ";
      fi
    dmeth_cmd="${dmeth_cmd}"" > ${diffmeth_out_name}";
;;

  CpGlist)
# need a valid input bam/sam file
    if [[ -z "${file_list1[0]}" ]]; then
      printf "%s operation needs a 'file_list1[0]'\n" ${diffmeth_run_type};
      exit 1;
    else
      if [[ $verbose == "yes" ]]; then
        printf " listing CpGs for '%s'\n" "${file_list1[0]}";
      fi
      dmeth_cmd="${path_to_dmap}""diffmeth -G ""${dmap_chr_info_file}""${diffmeth_run_options}"" -E ""$diffmeth_fragment_min"",""$diffmeth_fragment_max";
      append_sam_or_bam "-R" "${file_list1[0]}";
      dmeth_cmd="${dmeth_cmd}"" > ${diffmeth_out_name}";
    fi
;;

  binlist)
# need a valid input bam/sam file
    if [[ -z "${file_list1[0]}" ]]; then
      printf "%s operation needs a 'file_list1[0]'\n" ${diffmeth_run_type};
      exit 1;
    else
      if [[ $verbose == "yes" ]]; then
        printf " listing CpGs for '%s'\n" "${file_list1[0]}";
      fi
      dmeth_cmd="${path_to_dmap}""diffmeth -G ""${dmap_chr_info_file}""${diffmeth_run_options}"" -L ""$diffmeth_fragment_min"",""$diffmeth_fragment_max";
      append_sam_or_bam "-R" "${file_list1[0]}";
      dmeth_cmd="${dmeth_cmd}"" > ${diffmeth_out_name}";
    fi
;;

  *)
  printf "Invalid diffmeth run mode: '%s'\n" "$diffmeth_run_type";
  exit 1;
  ;;
esac

if [[ $verbose == "yes" ]]; then
  printf "Executing: %s\n" "${dmeth_cmd}";
fi

eval "${dmeth_cmd}";
  
exit 0;
