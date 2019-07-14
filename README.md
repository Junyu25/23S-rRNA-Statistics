# 23S-rRNA-Statistics

* SILVA_EZ_GG_database_merge.sh 数据库合并的流程  
  
* mark_ref_order.R 对unaligned的一条序列进行位置注释，是保守区统计的预处理的一个步骤，并不在流程中出现
  

* class_SILVA_aligned_transform.py 用于把SILVA的aligned initial拆分成3个文件再用于统计，注意：每个文件的结尾缺少一个换行符，所以如果直接文件与文件拼接的话，会少了一行。这个程序还会把.都改成-

* conserved_region_summary.sh 计算好ATCG出现频率以后，统计保守区的整个流程

* summary_aligned_ref_length.py 存粹就是检查aligned_ref的长度的，用于做整个保守区统计的流程conserved_region_summary.sh中的第一步
  
* count_freq_micro_database.R 计算aligned_silva的每个位置的ATCG出现频率，用于做整个保守区统计的流程conserved_region_summary.sh中的第二步

* summary_freq.R 对具体的ATCG比例的连贯性进行，整个保守区统计的流程conserved_region_summary.sh中的第三步，也是核心。


---

* primer_coverage_estimator.sh 统计单一引物的覆盖度情况的整体流程

* coverage_calculate.R 计算引物覆盖数和整体背景的比例，用于引物覆盖度计算primer_coverage_estimator.sh的流程中的2.1步骤

* extract_taxa_by_faID.R 利用fa文件里面的ID去获取全长ID，最后获取实际的taxonomy列表，是文件间比对的任务，用于引物覆盖度计算primer_coverage_estimator.sh的流程中的2.1步骤


* primer_pair_coverage_estimator.sh 统计单一引物对的覆盖度情况的整体流程

* intersect_uniq_taxaID.R uniqID和taxaID找交集，是一个比较特殊情况下用的功能，SILVA_EZ_GG_database_merge.sh数据库合并的流程中的3.2步骤


* SILVA_conserved_region_summary_loop_phylum.sh 对SILVAalign的每个门下面的每个界进行保守度统计的命令

