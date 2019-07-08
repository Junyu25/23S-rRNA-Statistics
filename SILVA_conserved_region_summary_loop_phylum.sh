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

#This script will run the ATCG freq database and the conserved region searching automatically for each phylum in a for loop 

#self contain command
#Library file requirement:
#NOTE1: Rscript must be activated, generally, can activate QIIME

if [ $# -ne 5 ]
then
  echo ""
    echo "Usage: SILVA_conserved_region_summary_loop_phylum.sh input_full_aligned_fasta Ecoli_MG1655_make_order step_width script_folder input_tax"
    echo "Example:sh /thinker/net/biosoft/ytan_tools/SILVA_conserved_region_summary_loop_phylum.sh Bacteria.fa Ecoli.K-12substr.MG1655.fa_mark_order.txt 7 /thinker/net/biosoft/ytan_tools/ Bacteria.tax"
    echo ""
    echo "input_full_aligned_fasta - The full-aligned fa."
    echo "Ecoli_MG1655_make_order - The path to Ecoli_MG1655_make_order file."
    echo "step_width - step_width used for summary, default is 7."
    echo "script_folder - the folder of scripts."
    echo "input_tax - The tax file of the full-aligned fa."
    
       
    exit 1
fi

#name the parameters
fa_f=$1
Ecoli_f=$2
step_width=$3
script_folder=$4
tax_f=$5

#check files
if [ ! -s $fa_f ] 
then
  echo ""
  echo "Warning: The file $fa_f does not exist in SILVA_conserved_region_summary_loop_phylum.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $tax_f ] 
then
  echo ""
  echo "Warning: The file $tax_f does not exist in SILVA_conserved_region_summary_loop_phylum.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $Ecoli_f ] 
then
  echo ""
  echo "Warning: The file $Ecoli_f does not exist in SILVA_conserved_region_summary_loop_phylum.sh, exit."
  echo ""
  exit 1
fi

#Check folders
if [ ! -d $script_folder ] 
then
  echo ""
  echo "Warning: The directory $script_folder does not exist in SILVA_conserved_region_summary_loop_phylum.sh, exit."
  echo ""
  exit 1
fi

################
echo "step 1, get the phylum name from the tax file."
echo `date`
###############
cut -f2 $tax_f | cut -f2 -d ";" | sort -u > $tax_f"_phylum.txt"

for i in $(cat $tax_f"_phylum.txt")
do 
    #for samtools, the region file must be the chr:from-to format or only the chr, as a result, here, only the readname without > is used. 
    grep $i";" $tax_f | cut -f1 | cut -f2 -d">" > $i"_ID.txt"

    if [ ! -s $i"_ID.txt" ] 
    then
      echo ""
      echo "Warning: The file $i _ID.txt did not generated in SILVA_conserved_region_summary_loop_phylum.sh, exit."
      echo ""
      break
    fi
    
    ################
    echo "step 2, extract corresponding fa file on $i ."
    echo `date`
    ###############
    #for samtools, the region file must be the chr:from-to format or only the chr, as a result, here, only the readname without > is used.
    #in order to have full line per read, must use -n option
    #in order to recognize the name of kindom and phylum,use two IDs in a file name
    samtools faidx $fa_f -r $i"_ID.txt" -n 100000 -o $fa_f"_"$i".fa"

    if [ ! -s $fa_f"_"$i".fa" ] 
    then
      echo ""
      echo "Warning: The file $fa_f"_"$i .fa did not generated in SILVA_conserved_region_summary_loop_phylum.sh, exit."
      echo ""
      break
    fi
    
    ################
    echo "step 3, run conserved_region_summary.sh on the phylum $i fa."
    echo `date`
    ###############
    sh $script_folder"/conserved_region_summary.sh" $fa_f"_"$i".fa" $Ecoli_f $step_width $script_folder

done