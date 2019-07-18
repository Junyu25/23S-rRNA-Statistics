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

#usage: Rscript count_freq_micro_database.R file_name=Bacteria.fa 
#use aligned_SILVA file as input and calculate frequency of each base and output into a table.
#the output will be the input of summray_freq.R

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
if (!exists("file_name")) {
	stop("\nRscript count_freq_micro_database.R file_name=Bacteria.fa\nWarning: Usage: file_name file is not exist, please check the path. \n\n")
}

con <- file(file_name, "r")

#generate a matrix which have 149999 col and 6 row
col_num=149999
out_matrix <- matrix(0, nrow=6,ncol=col_num)
rownames(out_matrix)=c("A","T","C","G","-",".")
colnames(out_matrix)=1:col_num

line=readLines(con,n=1)

while( length(line) != 0 ) {
     #timestart<-Sys.time()
     if (substring(line, 1, 1)!=">")
     {
          #deal with the string one time only, in order to save time (saved 200 fold time)
          str_list=strsplit(line, split="")[[1]]
          
          #check each base and add into the matrix
          for(i in 1:col_num){
               if(str_list[i]=="A")
               {
                    out_matrix[1,i]=out_matrix[1,i]+1
               } else if(str_list[i]=="T")
               {
                    out_matrix[2,i]=out_matrix[2,i]+1
               } else if(str_list[i]=="C")
               {
                    out_matrix[3,i]=out_matrix[3,i]+1
               } else if(str_list[i]=="G")
               {
                    out_matrix[4,i]=out_matrix[4,i]+1
               } else if(str_list[i]=="-")
               {
                    out_matrix[5,i]=out_matrix[5,i]+1
               } else if(str_list[i]==".")
               {
                    out_matrix[6,i]=out_matrix[6,i]+1
               }
          }
     }
     
     #timeend<-Sys.time()
     #runningtime<-timeend-timestart
     #print(runningtime)
     
     #go to next line
     line=readLines(con,n=1)
}

write.table(out_matrix, file=paste(file_name,"_freq_summary.txt",sep = ""),sep="\t", row.names= TRUE )

close(con)