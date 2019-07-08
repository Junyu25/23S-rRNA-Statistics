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

#This script will run the primer coverage estimator based on a list of primers in a for loop 
#a primer file for each primer will be generated

#self contain command
#Library file requirement:
#NOTE1: Rscript must be activated, generally, can activate QIIME

if [ $# -ne 7 ]
then
  echo ""
    echo "Usage: primer_coverage_estimator_loop_on_primerlist.sh input_unaligned_fasta primer_list script_folder input_strain_tax input_unaligned_fasta_full_ID all_ID_folder fa_f_sub "
    echo "Example:sh /thinker/net/biosoft/ytan_tools/primer_coverage_estimator_loop_on_primerlist.sh /thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/SILVA132-EZMay2018-GG13_8_merged_grouped.fa Probebase-201811-key-probes.txt /thinker/net/biosoft/ytan_tools/ /thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/SILVA132-EZMay2018-GG13_8_merged_7LV_uniq_strain.tax /thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/SILVA132-EZMay2018-GG13_8_merged_grouped_full_ID.txt /thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/intermediate/ SILVA132-EZMay2018-GG13_8_merged_grouped" 
    echo ""
    echo "input_unaligned_fasta - The unaligned fa, all Us must be convered to Ts already."
    echo "primer_list - The primer-list info from probase."
    echo "script_folder - the folder of scripts."
    echo "input_strain_tax - the strain tax file (such as SILVA132-EZMay2018-GG13_8_merged_7LV_uniq_strain.tax) of the input_unaligned_fasta."
    echo "input_unaligned_fasta_full_ID - the full ID file (such as SILVA132-EZMay2018-GG13_8_merged_grouped_full_ID.txt) of the input_unaligned_fasta."
    echo "all_ID_folder - the folder path of all_Dx_id_txt of all level"
    echo "fa_f_sub - the sub-str of fa file name：(such as SILVA132-EZMay2018-GG13_8_merged_grouped)"
    exit 1
fi

#name the parameters
fa_f=$1
primer_list=$2
script_folder=$3
tax_f=$4
fa_full_ID=$5
all_ID_folder=$6
fa_f_sub=$7

#check files
if [ ! -s $fa_f ] 
then
  echo ""
  echo "Warning: The file $fa_f does not exist in primer_coverage_estimator_loop_on_primerlist.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $tax_f ] 
then
  echo ""
  echo "Warning: The file $tax_f does not exist in primer_coverage_estimator_loop_on_primerlist.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $Ecoli_f ] 
then
  echo ""
  echo "Warning: The file $Ecoli_f does not exist in primer_coverage_estimator_loop_on_primerlist.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $primer_list ] 
then
  echo ""
  echo "Warning: The file $primer_list does not exist in primer_coverage_estimator_loop_on_primerlist.sh, exit."
  echo ""
  exit 1
fi

#Check folders
if [ ! -d $script_folder ] 
then
  echo ""
  echo "Warning: The directory $script_folder does not exist in primer_coverage_estimator_loop_on_primerlist.sh, exit."
  echo ""
  exit 1
fi

################
echo "step 1, get only the shortname and the seq from the primer_list file."
echo `date`
###############
cut -f1 $primer_list | sort -u > $primer_list"_ID.txt"
cut -f1,4 $primer_list | sort -u > $primer_list"_ID_seq.txt"

for i in $(cat $primer_list"_ID.txt")
do 
    #for samtools, the region file must be the chr:from-to format or only the chr, as a result, here, only the readname without > is used. 
    grep $i $primer_list"_ID_seq.txt" > "Primers_"$i".txt"

    if [ ! -s "Primers_"$i".txt" ] 
    then
      echo ""
      echo "Warning: The file Primers_ $i .txt did not generated in primer_coverage_estimator_loop_on_primerlist.sh, exit."
      echo ""
      exit 1
    fi
   
    ################
    echo "step 2, run primer_coverage_estimator.sh on the Primers_ $i .txt."
    echo `date`
    ###############
    sh $script_folder"/primer_coverage_estimator.sh" $fa_f "Primers_"$i".txt" $i $script_folder $tax_f $fa_full_ID $all_ID_folder $fa_f_sub

done