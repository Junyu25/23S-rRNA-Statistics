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

#This script will run the in silico primer matching and statistics automatically

#This script starts from unaligned full class fa
#Library file requirement:
#NOTE1:  QIIME must be activated
#recommended workding directory is the primer_coverage folder

if [ $# -ne 8 ]
then
  echo ""
    echo "Usage: primer_coverage_estimator.sh input_unaligned_fasta primer_f primer_r script_folder input_strain_tax input_unaligned_fasta_full_ID all_ID_folder fa_f_sub "
    echo "Example:sh /thinker/net/biosoft/ytan_tools/primer_coverage_estimator.sh /thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/SILVA132-EZMay2018-GG13_8_merged_grouped.fa /thinker/net/ytan/microbiome_proj/Primer_prospector_test/primers/Primers_515Yf.txt 515Yf /thinker/net/biosoft/ytan_tools/ /thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/SILVA132-EZMay2018-GG13_8_merged_7LV_uniq_strain.tax /thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/SILVA132-EZMay2018-GG13_8_merged_grouped_full_ID.txt /thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/intermediate/ SILVA132-EZMay2018-GG13_8_merged_grouped"
    echo ""
    echo "input_unaligned_fasta - The unaligned fa, all Us must be convered to Ts already."
    echo "primer - The path to the forward primer."
    echo "primer_name - The name of the forward primer."
    echo "script_folder - the folder of scripts."
    echo "input_strain_tax - the strain tax file (such as SILVA132-EZMay2018-GG13_8_merged_7LV_uniq_strain.tax) of the input_unaligned_fasta."
    echo "input_unaligned_fasta_full_ID - the full ID file (such as SILVA132-EZMay2018-GG13_8_merged_grouped_full_ID.txt) of the input_unaligned_fasta."
    echo "all_ID_folder - the folder path of all_Dx_id_txt of all level"
    echo "fa_f_sub - the sub-str of fa file name：(such as SILVA132-EZMay2018-GG13_8_merged_grouped)"
    exit 1
fi

#name the parameters
fa_f=$1
primer_f=$2
primer_f_name=$3
script_folder=$4
tax_f=$5
fa_full_ID=$6
all_ID_folder=$7
fa_f_sub=$8

#check files
if [ ! -s $fa_f ] 
then
  echo ""
  echo "Warning: The file $fa_f does not exist in primer_coverage_estimator.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $primer_f ] 
then
  echo ""
  echo "Warning: The file $primer_f does not exist in primer_coverage_estimator.sh, exit."
  echo ""
  exit 1
fi

#Check folders
if [ ! -d $script_folder ] 
then
  echo ""
  echo "Warning: The directory $script_folder does not exist in primer_coverage_estimator.sh, exit."
  echo ""
  exit 1
fi

################
echo "step 1, analyze_primers on target fa using Primer seqs."
echo `date`
###############
if [ ! -s $primer_f_name"_"$fa_f_sub"_hits.txt" ] 
then
  analyze_primers.py -f $fa_f -P $primer_f
fi

if [ ! -s $primer_f_name"_"$fa_f_sub"_hits.txt" ] 
then
  echo ""
  echo "Warning: The file $primer_f_name _ $fa_f_sub _hits.txt was not generated in primer_coverage_estimator.sh step2 failed, exit."
  echo ""
  exit 1
fi

################
echo "step 2, get amplicons using primer_hits."
echo `date`
###############
#get the output name format
echo $primer_f_name"_"$fa_f_sub"_hits.txt"

#get fa
out_primer_f=$primer_f_name"/"$primer_f_name

if [ ! -s $out_primer_f"_amplicons.fasta" ] 
then
  get_amplicons_and_reads.py -f $fa_f -i $primer_f_name"_"$fa_f_sub"_hits.txt" -o $primer_f_name -R 250
fi

if [ ! -s $out_primer_f"_amplicons.fasta" ] 
then
  echo ""
  echo "Warning: The file $out_primer_f _amplicons.fasta was not generated in primer_coverage_estimator.sh step2 failed, exit."
  echo ""
  exit 1
fi

################
echo "step 2.1, do coverage statistics on primers."
echo `date`
###############
function coverage_stat()
{
  #extract ID from amplicons.fasta
  grep "^>" $1"_amplicons.fasta" > $1"_ID.txt"
  #get taxa from faID
  Rscript $4"/extract_taxa_by_faID.R" file.list1=$1"_ID.txt" file.list2=$2 file.list3=$3
  out_taxa=$1"_ID.txt_merge_taxa.txt"
  
  #do stat on all levels
  allID_count $6 $5"/all_D0_id_txt" $out_taxa "D0" 
  allID_count $6 $5"/all_D1_id_txt" $out_taxa "D1"
  allID_count $6 $5"/all_D2_id_txt" $out_taxa "D2"
  allID_count $6 $5"/all_D3_id_txt" $out_taxa "D3"
  allID_count $6 $5"/all_D4_id_txt" $out_taxa "D4"
  allID_count $6 $5"/all_D5_id_txt" $out_taxa "D5"
  #in order to find the missed species
  allID_count $6 $5"/all_D6_id_txt" $out_taxa "D6"
  
  Rscript $4"/coverage_calculate.R" file.list1=$6"/D0_stat_all.txt" file.list2=$5"/D0/D0_stat_all.txt"
  Rscript $4"/coverage_calculate.R" file.list1=$6"/D1_stat_all.txt" file.list2=$5"/D1/D1_stat_all.txt"
  Rscript $4"/coverage_calculate.R" file.list1=$6"/D2_stat_all.txt" file.list2=$5"/D2/D2_stat_all.txt"
  Rscript $4"/coverage_calculate.R" file.list1=$6"/D3_stat_all.txt" file.list2=$5"/D3/D3_stat_all.txt"
  Rscript $4"/coverage_calculate.R" file.list1=$6"/D4_stat_all.txt" file.list2=$5"/D4/D4_stat_all.txt"
  Rscript $4"/coverage_calculate.R" file.list1=$6"/D5_stat_all.txt" file.list2=$5"/D5/D5_stat_all.txt"
}

function allID_count()
{
    #check folder
    if [ ! -d $1 ] 
    then
      echo ""
      echo "Warning: The directory $1 does not exist in primer_coverage_estimator.sh, generate it."
      echo ""
      mkdir -p $1
    fi
    cat $2 | while read line
    do
        uniq_species=`cut -f2 -d" " $3 | grep -F ${line}| cut -f7 -d ";" |cut -f1 -d ":" | sort -u | wc -l `
        uniq_strain=`cut -f2 -d" " $3 | grep -F ${line}| cut -f7 -d ";" | sort -u | wc -l `
        echo -e ${line}"\t"${uniq_species}"\t"${uniq_strain}
    done > $1"/"$4"_stat_all.txt"
}


#on primers, need to know their hit file structure and extract ID
coverage_stat $out_primer_f $fa_full_ID $tax_f $script_folder $all_ID_folder $primer_f_name"/"

#The following steps are for a primer-pair
#################
#echo "step 3, use vsearch to collapse redundance."
#echo `date`
################
#/thinker/net/biosoft/miniconda2/pkgs/vsearch-2.4.3-0/bin/vsearch --derep_fulllength $out_primer_f"_amplicons.fasta" --output $out_primer_f"_amplicons_uniq.fna" --uc $out_primer_f"_amplicons.fasta_uniq.uc" --threads 10
#
#if [ ! -s $out_primer_f"_amplicons.fasta_uniq.uc" ] 
#then
#  echo ""
#  echo "Warning: The file $out_primer_f _amplicons.fasta_uniq.uc was not generated in primer_coverage_estimator.sh step3 failed, exit."
#  echo ""
#  exit 1
#fi
#
#################
#echo "step 4, get the ID information of each IU."
#echo `date`
################
#python /thinker/net/biosoft/gist-github-qiime/parse_otu_mapping_from_uc.py $out_primer_f"_amplicons.fasta_uniq.uc" $out_primer_f"_IU_mapping.txt"
#
#if [ ! -s $out_primer_f"_IU_mapping.txt" ] 
#then
#  echo ""
#  echo "Warning: The file $out_primer_f _IU_mapping.txt was not generated in primer_coverage_estimator.sh step4 failed, exit."
#  echo ""
#  exit 1
#fi
#
#################
#echo "step 5, get the representative fa of each IU."
#echo `date`
################
#pick_rep_set.py -f $out_primer_f"_amplicons_uniq.fna" -i $out_primer_f"_IU_mapping.txt" -o $out_primer_f"_IU_grouped.fa"
#
#if [ ! -s $out_header"_IU_grouped.fa" ] 
#then
#  echo ""
#  echo "Warning: The file $out_header _IU_grouped.fa was not generated in primer_coverage_estimator.sh step5 failed, exit."
#  echo ""
#  exit 1
#fi

#########
#echo "step 5.1, do coverage statistics on IUs."
#echo `date`
###############
#this step is not necessary for a primer
#这里需要_IU_mapping.txt 515Yf_907r_IU_grouped.fa的ID $fa_full_ID $tax_f,似乎有点混乱了。。。要整理一下思路

