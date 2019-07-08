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


#usage: Rscript summary_freq.R count_freq_in=Bacteria.fa_freq_summary.txt ref_matrix_file=Ecoli_MS107-1.fa_mark_order.txt step_width=7(default)
#This script is to analyze the frequency summary table
#Three steps in all:
#Step 1: calculate the frequency percentage and transfer it into a prob and a max matrix, the prob matrix can be used to calculate the ATCG plot.
#Step 2: In the max matrix, find the high conserved regions by stepwidth, and mark it in the summary table when it meet the cutoff.
#Step 3: Using the summary table ,output the sorted result by region length, and annotate the location and conresponding sequence on Ecoli.

step_width=7

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
if (!exists("count_freq_in")) {
	stop("\nRscript summray_freq.R count_freq_in=Bacteria.fa_freq_summary.txt ref_matrix_file=Ecoli_MS107-1.fa_mark_order.txt step_width=7(default)\nWarning: Usage: count_freq_in file is not exist, please check the path. \n\n")
}

print(step_width)
step_width=as.numeric(step_width)

#step 1
#read in freq table
freq_matrix<-read.delim2(file=count_freq_in,header=TRUE,stringsAsFactors=FALSE)
col_num=dim(freq_matrix)[2]

prob_matrix <- matrix(0, nrow=5,ncol=col_num)
rownames(prob_matrix) <- c("A","T","C","G","-")
colnames(prob_matrix) <- 1:col_num

max_matrix <- as.data.frame(matrix(0, ncol=12,nrow=col_num))
colnames(max_matrix) <- c("Max_prob","100","99","97","95","90","25","40","50","60","70","Max_Char")
rownames(max_matrix) <- 1:col_num

#calculate frequency
for(i in 1:col_num){
    sum_num <- sum(freq_matrix[,i])
    for (j in 1:5){
        prob_matrix[j,i]<-freq_matrix[j,i]/sum_num
    }
    
    #get the max prob and the coresponding char
    max_i <- max(prob_matrix[1:5,i])
    max_char <- names(which(prob_matrix[,i]==max(prob_matrix[1:5,i])))
    
    max_matrix[i,1] <- max_i
    #in order to aviod the situation of equal max_prob and report error,just take the first charater
    max_matrix[i,12] <- max_char[1]
    if(max_i==1) {
        max_matrix[i,2:6]=1
    } else if(max_i>=0.99) {
        max_matrix[i,3:6]=1
    } else if(max_i>=0.97) {
        max_matrix[i,4:6]=1
    } else if(max_i>=0.95) {
        max_matrix[i,5:6]=1
    } else if(max_i>=0.90) {
        max_matrix[i,6]=1
    } 
}

#check the density of key characters
den_matrix <- matrix(0, nrow=10,ncol=2)
rownames(den_matrix) <- c("0.3","0.4","0.5","0.6","0.7","0.8","0.9","0.95","0.99","1")
colnames(den_matrix) <- c("num of smaller and equal than","num of in between")
for(rn in 1:10){
	den_matrix[rn,1] = length(which(max_matrix[which(max_matrix[,12]!="-"),1]<=as.numeric(rownames(den_matrix)[rn])))
}

#get the in between number
den_matrix[1,2]=den_matrix[1,1]
for(rn in 2:10){
	den_matrix[rn,2] = den_matrix[rn,1]-den_matrix[rn-1,1]
}
write.table(den_matrix, file=paste(count_freq_in,"_density.txt",sep=""),sep="\t", col.names= TRUE, row.names= TRUE, quote = FALSE)


#get the various region
max_matrix[which(max_matrix[,1]<0.25),7]=1
max_matrix[which(max_matrix[,1]<0.4),8]=1
max_matrix[which(max_matrix[,1]<0.5),9]=1
max_matrix[which(max_matrix[,1]<0.6),10]=1
max_matrix[which(max_matrix[,1]<0.7),11]=1

#fill in the "-" the gaps as various
max_matrix[which(max_matrix[,12]=="-"),7:11]=1

ref_matrix<-read.delim2(file=ref_matrix_file,header=TRUE,stringsAsFactors=FALSE)
prob_matrix_out <- matrix(0, nrow=6,ncol=col_num)
prob_matrix_out[1:5,] <- prob_matrix
prob_matrix_out[6,]  <- unlist(ref_matrix[2,])
rownames(prob_matrix_out) <- c("A","T","C","G","-","E.coli_loci")

write.table(t(prob_matrix_out), file=paste(count_freq_in,"_prob.txt",sep=""),sep="\t", col.names= TRUE )
write.table(max_matrix, file=paste(count_freq_in,"_max_matrix.txt",sep=""),sep="\t", row.names= TRUE, quote = FALSE )

#Step 2 In the max matrix, find the high conserved regions by stepwidth, and mark it in the summary table when it meet the cutoff.
summary_matrix <- matrix(0, nrow=5,ncol=col_num)
rownames(summary_matrix) <- c("100","99","97","95","90")
colnames(summary_matrix) <- 1:col_num

#build summary matix to record the length of short conserved regions
for (j in 2:6){
	i=1
	min_prob <- as.numeric(rownames(summary_matrix)[j-1])/100
	while(i <=(col_num-step_width)){
		#check whether this base is over the meaningful persentage 
		#if(max_matrix[j,i]>=min_prob){
		if(max_matrix[i,j]==1 & max_matrix[i,12]!="-"){
			seq_begin=i
			seq_end=seq_begin+step_width-1
			sum_count_ori = sum(max_matrix[seq_begin:seq_end,j])
			sum_count = sum(max_matrix[seq_begin:seq_end,j])
			#cut_off <-(step_width-2)*as.numeric(rownames(summary_matrix)[j-1])/100
			#cut_off <-(step_width-1) since most the fails are truely caused by variables, 1 exception quota is not necessary.
			cut_off <-step_width
			#check whether it can be enlongated
			while( sum_count >= cut_off & seq_end<=col_num) {
				seq_begin=seq_begin+1
				seq_end=seq_begin+step_width-1
				sum_count = sum(max_matrix[seq_begin:seq_end,j])
			}
			#the last base do not need to be counted, so +1 is not neccessary
			seq_len=seq_end-i
			
			#check whether it can be recorded			
			if(sum_count_ori >= cut_off){
				summary_matrix[j-1,i] <- seq_len
			
				#reset the i to the connective round.
				i=seq_begin
			} else {i=i+1}
		} else {i=i+1}
	}
}



#Step 3
#using the summary matrix, generate groups of conserved regions seqs(Ecoli as example), sort by length,then output base_num and short seq
#optional, i should add the position of base num on Ecoli in the output.
for (prob_i in 1:5){
	outfile=paste(count_freq_in,"_",rownames(summary_matrix)[prob_i],".txt",sep="")
	outfile1=paste(count_freq_in,"_",rownames(summary_matrix)[prob_i],"_details.txt",sep="")
	outfile2=paste(count_freq_in,"_",rownames(summary_matrix)[prob_i],"_stat.txt",sep="")
	len_list <- sort(unique(summary_matrix[prob_i,]), decreasing = TRUE)
	
	#if (file.exists(outfile)) {
	#	file.create(outfile)
	#}
	#if (file.exists(outfile1)) {
	#	file.create(outfile1)
	#}
	#if (file.exists(outfile2)) {
	#	file.create(outfile2)
	#}
	
	write("#header:len start_pos-end_pos Ecoli-seq max_prob Max_char_seq",file = outfile,append = FALSE)
	write("#header:len x-meanning y-meaning:Max_char_seq prob_matrix",file = outfile1,append = FALSE)
	write("#header:len start_pos-end_pos len_without_meanninglessgap;len_without_gap;unaligned_start-begin       meaningless gaps are - positions that are 100%",file = outfile2,append = FALSE)
	
	prob_group <- as.numeric(rownames(summary_matrix)[prob_i])/100
	for (len_i in len_list){
		if(len_i>=step_width){
			head_str <- c(paste(len_i,"bp ",sep=""))
			out_str <- head_str
			group_list <- which(summary_matrix[prob_i,]==len_i)
			for (seq_begin_i in group_list){
				seq_end <- seq_begin_i + summary_matrix[prob_i,seq_begin_i] -1
				seq_str <- paste(seq_begin_i,"-",seq_end,":", paste(ref_matrix[1, seq_begin_i:seq_end],collapse=""),"|",paste(max_matrix[seq_begin_i:seq_end,1],collapse=","),"|",paste(max_matrix[seq_begin_i:seq_end,12],collapse=""),";",sep="")
				out_str <- paste(out_str, seq_str, sep="")
				out_str1 <- paste(head_str, "x:ATCG-", paste("y:",seq_begin_i,"-",seq_end,":",paste(max_matrix[seq_begin_i:seq_end,12],collapse=""),sep=""), sep=" ")
				out_str2 <- paste(head_str, seq_begin_i,"-",seq_end,":", len_i-sum(prob_matrix[5,seq_begin_i:seq_end]==1),";", len_i-sum(prob_matrix[5,seq_begin_i:seq_end]>=prob_group), ";", paste(ref_matrix[2, seq_begin_i], ref_matrix[2,seq_end],sep="-"), sep="")
				write(out_str1,file = outfile1,append = TRUE)
				write(prob_matrix[,seq_begin_i:seq_end],file = outfile1,append = TRUE)
				write(out_str2,file = outfile2,append = TRUE)
				
			}
			write(out_str,file = outfile,append = TRUE)
		}
	}
	stat_matrix<-read.delim2(file=outfile2,header=TRUE,stringsAsFactors=FALSE,sep = ";")
	outfile3=paste(count_freq_in,"_",rownames(summary_matrix)[prob_i],"_highlight.txt",sep="")
	highlight_matirx <- stat_matrix[which(stat_matrix[,2]>=11),2:3]
	highlight_out <- highlight_matirx[order(highlight_matirx[,1], decreasing = TRUE),]
	write.table(highlight_out, file=outfile3,sep="\t", row.names = FALSE, col.names= FALSE, quote = FALSE )
}


#Step 4 In the max matrix, find the high various regions by stepwidth, and mark it in the summary table when it meet the cutoff.
summary_var_matrix <- matrix(0, nrow=5,ncol=col_num)
rownames(summary_var_matrix) <- c("25","40","50","60","70")
colnames(summary_var_matrix) <- 1:col_num

#build summary matix to record the length of short various regions
for (j in 7:11){
	i=1
	min_prob <- as.numeric(rownames(summary_var_matrix)[j-6])/100
	while(i <=(col_num-step_width)){
		#check whether this base is over the meaningful persentage 
		#if(max_matrix[j,i]>=min_prob){
		if(max_matrix[i,j]==1 & max_matrix[i,12]!="-"){
			seq_begin=i
			seq_end=seq_begin+step_width-1
			sum_count_ori = sum(max_matrix[seq_begin:seq_end,j])
			sum_count = sum(max_matrix[seq_begin:seq_end,j])
			#cut_off <-(step_width-2)*as.numeric(rownames(summary_var_matrix)[j-1])/100
			#cut_off <-(step_width-1) since most the fails are truely caused by variables, 1 exception quota is not necessary.
			cut_off <-step_width
			##add a check point of "-" percentage 
			#permax_matrix[i,12]!="-"
			#check whether it can be enlongated
			while( sum_count >= cut_off & seq_end<=col_num) {
				seq_begin=seq_begin+1
				seq_end=seq_begin+step_width-1
				sum_count = sum(max_matrix[seq_begin:seq_end,j])
			}
			#the last base do not need to be counted, so +1 is not neccessary
			seq_len=seq_end-i
			
			#check whether it can be recorded			
			if(sum_count_ori >= cut_off){
				summary_var_matrix[j-6,i] <- seq_len
			
				#reset the i to the connective round.
				i=seq_begin
			} else {i=i+1}
		} else {i=i+1}
	}
}

#Step 5
#using the summary_var_matrix, generate groups of highly various regions seqs(Ecoli as example), sort by length,then output base_num and short seq
#optional, i should add the position of base num on Ecoli in the output.


for (prob_i in 1:5){
	outfile=paste(count_freq_in,"_",rownames(summary_var_matrix)[prob_i],".txt",sep="")
	outfile1=paste(count_freq_in,"_",rownames(summary_var_matrix)[prob_i],"_details.txt",sep="")
	outfile2=paste(count_freq_in,"_",rownames(summary_var_matrix)[prob_i],"_stat.txt",sep="")
	len_list <- sort(unique(summary_var_matrix[prob_i,]), decreasing = TRUE)
	
	#if (file.exists(outfile)) {
	#	file.create(outfile)
	#}
	#if (file.exists(outfile1)) {
	#	file.create(outfile1)
	#}
	#if (file.exists(outfile2)) {
	#	file.create(outfile2)
	#}
	
	write("#header:len start_pos-end_pos Ecoli-seq max_prob Max_char_seq",file = outfile,append = FALSE)
	write("#header:len x-meanning y-meaning:Max_char_seq prob_matrix",file = outfile1,append = FALSE)
	write("#header:len start_pos-end_pos len_without_meanninglessgap;len_without_gap;unaligned_start-begin       meaningless gaps are - positions that are 100%",file = outfile2,append = FALSE)
	
	prob_group <- as.numeric(rownames(summary_var_matrix)[prob_i])/100
	for (len_i in len_list){
		if(len_i>=step_width){
			head_str <- c(paste(len_i,"bp ",sep=""))
			out_str <- head_str
			group_list <- which(summary_var_matrix[prob_i,]==len_i)
			for (seq_begin_i in group_list){
				seq_end <- seq_begin_i + summary_var_matrix[prob_i,seq_begin_i] -1
				seq_str <- paste(seq_begin_i,"-",seq_end,":", paste(ref_matrix[1, seq_begin_i:seq_end],collapse=""),"|",paste(max_matrix[seq_begin_i:seq_end,1],collapse=","),"|",paste(max_matrix[seq_begin_i:seq_end,12],collapse=""),";",sep="")
				out_str <- paste(out_str, seq_str, sep="")
				out_str1 <- paste(head_str, "x:ATCG-", paste("y:",seq_begin_i,"-",seq_end,":",paste(max_matrix[seq_begin_i:seq_end,12],collapse=""),sep=""), sep=" ")
				#out_str2 <- paste(head_str, seq_begin_i,"-",seq_end,":", len_i-sum(prob_matrix[5,seq_begin_i:seq_end]==1),";", len_i-sum(prob_matrix[5,seq_begin_i:seq_end]>=prob_group), ";", paste(ref_matrix[2, seq_begin_i], ref_matrix[2,seq_end],sep="-"), sep="")
				#because the logic of calculation of various region is diff from the conserved region,the way of annotation must be adjusted.
				out_str2 <- paste(head_str, seq_begin_i,"-",seq_end,":", len_i-sum(max_matrix[seq_begin_i:seq_end,12]=="-"),";", len_i-sum(ref_matrix[1,seq_begin_i:seq_end]=="-"), ";", paste(ref_matrix[2, seq_begin_i], ref_matrix[2,seq_end],sep="-"), sep="")
				write(out_str1,file = outfile1,append = TRUE)
				write(prob_matrix[,seq_begin_i:seq_end],file = outfile1,append = TRUE)
				write(out_str2,file = outfile2,append = TRUE)
				
			}
			write(out_str,file = outfile,append = TRUE)
		}
	}
	stat_matrix<-read.delim2(file=outfile2,header=TRUE,stringsAsFactors=FALSE,sep = ";")
	outfile3=paste(count_freq_in,"_",rownames(summary_var_matrix)[prob_i],"_highlight.txt",sep="")
	highlight_matirx <- stat_matrix[which(stat_matrix[,2]>=11),2:3]
	highlight_out <- highlight_matirx[order(highlight_matirx[,1], decreasing = TRUE),]
	write.table(highlight_out, file=outfile3,sep="\t", row.names = FALSE, col.names= FALSE, quote = FALSE )
}