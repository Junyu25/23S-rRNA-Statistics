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


#usage: Rscript mark_ref_order.R aligned_ref=Ecoli_MS107-1.fa 
#This script is to annotate the unaligned base order for the input ref

#check arguments
for (e in commandArgs()) {
        ta = strsplit(e,"=",fixed=TRUE)
        if(! is.na(ta[[1]][2])) {
                temp = ta[[1]][2]
                if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
                temp = as.integer(temp)
                }
        if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
                temp = as.numeric(temp)
                }
        assign(ta[[1]][1],temp)
        } else {
        assign(ta[[1]][1],TRUE)
        }
}

#check whether the file in is exist
if (!exists("aligned_ref")) {
	stop("\nRscript mark_ref_order.R aligned_ref=Ecoli_MS107-1.fa \nWarning: Usage: aligned_ref file is not exist, please check the path. \n\n")
}

print(aligned_ref)

ref_str <- read.delim2(file=aligned_ref,header=TRUE,stringsAsFactors=FALSE)

ref_str_list = strsplit(ref_str[1,1], split="")[[1]]

ref_num_list = ref_str_list

start_num = 0
for (n in 1:length(ref_str_list)){
    if (ref_str_list[n]=="-"){
        ref_num_list[n]=start_num
    } else {
        start_num = start_num+1
        ref_num_list[n]=start_num
    }
}

ref_matrix <- rbind(ref_str_list, ref_num_list)

write.table(ref_matrix, file=paste(aligned_ref,"_mark_order.txt",sep=""),sep="\t", row.names= TRUE )
