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
cwd=ENVIRON["PWD"];
}

$0~/>/{if (outfile != "")
  close(outfile);
outfile = sprintf("%s%s%s",namehdr,substr($1,2),extension);
printf("%s\n",$0) > outfile;
if (infofile != "")
  printf("\"%s\"\t\"%s/%s\"\n",substr($1,2),cwd,outfile) > infofile;
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
