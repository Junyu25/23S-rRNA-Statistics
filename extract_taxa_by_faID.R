##Generally use this R version: /share/apps/R/R-2.12.2_gnu412_x86_64_vanilla/bin/Rscript --vanilla 

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


#this script is to use the faID to get full_ID and then get taxa
#first, read in the three files
#second, merge file1 and 2 by the first col
#third, merge the step2 out and the file 3
#example: Rscript $script_folder"/extract_taxa_by_faID.R" file.list1="515Yf_926r/515Yf_926r_ID.txt" file.list2="/thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/SILVA132-EZMay2018-GG13_8_merged_grouped_full_ID.txt" file.list3="/thinker/net/biosoft/16S_ref_databases/SILVA132_EZ18_GG138/SILVA132-EZMay2018-GG13_8_merged_7LV_uniq_strain.tax"

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

#check whether file in is exist
if (!exists("file.list1")) {
	stop("\n\nWarning: Usage: file.list1( _ID.txt) is not exist, please check the path. \n\n")
}

#check whether file in is exist
if (!exists("file.list2")) {
	stop("\n\nWarning: Usage: file.list2(_grouped_full_ID.txt) is not exist, please check the path. \n\n")
}

#check whether file in is exist
if (!exists("file.list3")) {
	stop("\n\nWarning: Usage: file.list3(_7LV_uniq_strain.tax) is not exist, please check the path. \n\n")
}


#read two list of ID
list1<-read.delim2(file=file.list1,header=FALSE,stringsAsFactors=FALSE)
list2<-read.delim2(file=file.list2,header=FALSE,stringsAsFactors=FALSE,sep = " ")
list3<-read.delim2(file=file.list3,header=FALSE,stringsAsFactors=FALSE,sep = " ")

#list1<-as.matrix(list1)
#list2<-as.matrix(list2)


#merge by the first column (denovo ID)
m1 <- merge(list1, list2, by.x = "V1", by.y = "V1")
#colnames(m1)<- c("denovoID","TaxaName")
print(dim(m1))
#merge the taxa by ID
m2 <- merge(m1, list3, by.x = "V2", by.y = "V1")

#output the file
write.table(m2[,c(1,3)], file=paste(file.list1,"_merge_taxa.txt",sep = ""),  row.names = FALSE , col.names = FALSE, quote=FALSE)