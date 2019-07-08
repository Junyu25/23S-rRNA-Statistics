#!/bin/bash

#   Copyright {2015} Yuxiang Tan
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

#This script will run the ATCG freq database and the conserved region searching automatically

#This script starts from full aligned subclass fa
#Library file requirement:
#NOTE1: Rscript must be activated, generally, can activate QIIME

if [ $# -ne 4 ] #if the parameter is not equal to 4
then
  echo ""
    echo "Usage: conserved_region_summary.sh input_full_aligned_fasta Ecoli_MG1655_make_order step_width script_folder"
    echo "Example:sh /thinker/net/biosoft/ytan_tools/conserved_region_summary.sh Archaea.fa Ecoli.K-12substr.MG1655.fa_mark_order.txt 7 /thinker/net/biosoft/ytan_tools/"
    echo ""
    echo "input_full_aligned_fasta - The full-aligned fa."
    echo "Ecoli_MG1655_make_order - The path to Ecoli_MG1655_make_order file."
    echo "step_width - step_width used for summary, default is 7."
    echo "script_folder - the folder of scripts."
       
    exit 1
fi

#name the parameters
fa_f=$1 #Assign fa_f to the first parameter 
Ecoli_f=$2
step_width=$3
script_folder=$4

#check files
if [ ! -s $fa_f ] 
then
  echo ""
  echo "Warning: The file $fa_f does not exist in conserved_region_summary.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $Ecoli_f ] 
then
  echo ""
  echo "Warning: The file $Ecoli_f does not exist in conserved_region_summary.sh, exit."
  echo ""
  exit 1
fi

#Check folders
if [ ! -d $script_folder ] 
then
  echo ""
  echo "Warning: The directory $script_folder does not exist in conserved_region_summary.sh, exit."
  echo ""
  exit 1
fi

################
echo "step 1, check the fa file length."
echo `date`
###############
#call the command to check the length of aligned_ref
python $script_folder"/summary_aligned_ref_length.py" $fa_f 
length_sum_f="database_len_summary.txt" 
length_count=`cut -f2 -d" " $length_sum_f | sort -u | wc -l` 

if [ $length_count -ne 1 ] 
then
  echo ""
  echo "Warning: The file $fa_f has unequal length in conserved_region_summary.sh, exit."
  echo ""
  exit 1
fi

################
echo "step 2, count the ATCG frequency by R."
echo `date`
###############
#call the Rscript to calculate the ATCG frequency of every position
Rscript $script_folder"/count_freq_micro_database.R" file_name=$fa_f 

if [ ! -s $fa_f"_freq_summary.txt" ] 
then
  echo ""
  echo "Warning: The file $fa_f _freq_summary.txt did not generated in conserved_region_summary.sh, exit."
  echo ""
  exit 1
fi

################
echo "step 3, summary the ATCG frequency by R and find the highly conserved region."
echo `date`
###############
#Summary the ATCG frequency by R and find the highly conserved region
Rscript $script_folder"/summary_freq.R" count_freq_in=$fa_f"_freq_summary.txt" ref_matrix_file=$Ecoli_f step_width=$step_width
