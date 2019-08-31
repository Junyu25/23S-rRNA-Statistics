inf_name="Eukaryota.fasta_freq_summary.txt_prob.txt"
data <- read.table(file=inf_name,header=T, sep="\t",quote = "", stringsAsFactors=F)
data1 <- apply(data,c(1,2),as.numeric)
data2 <-data1[,1:4]
data_max <- apply(data2,1,max)

data_out <- data.frame(data_max,data1[,6])

#count insertion freq
freq_mtx <- table(data1[,6])
pdf(file=paste(inf_name,"_freq.pdf",sep = ""))
plot(freq_mtx)
dev.off()
write.table(freq_mtx,file=paste(inf_name,"_freq_matrix.txt",sep = ""), sep='\t')

#filter out the positions with rare inserts
data_out_filter <- data_out[which(data_out[,1]>0.2),]
freq_f_mtx <- table(data_out_filter[,2])
pdf(file=paste(inf_name,"_filtered_freq.pdf",sep = ""))
plot(freq_f_mtx)
dev.off()
write.table(freq_f_mtx,file=paste(inf_name,"_filtered_freq_matrix.txt",sep = ""), sep='\t')

pdf(file=paste(inf_name,"__max_filtered_boxplot.pdf",sep = ""))
boxplot(data_out_filter[,1]~data_out_filter[,2],data_out_filter)
dev.off()
write.table(data_out_filter,file=paste(inf_name,"_max_filtered.txt",sep = ""), sep='\t')


pdf(file=paste(inf_name,"__max_filtered_plot.pdf",sep = ""),width=35)
plot(data_out_filter[,1]~data_out_filter[,2],data_out_filter,type="l")
#axis(1, at = data_out_filter[,2],las=2)
dev.off()

#separate each 500 bases
data_out_filter_500 <- data_out_filter[which(data_out_filter[,2]<=500),]
data_out_filter_1000 <-  data_out_filter[intersect(which(data_out_filter[,2]>500) , which(data_out_filter[,2]<=1000)),]
data_out_filter_1500 <- data_out_filter[which(data_out_filter[,2]>1000),]

pdf(file=paste(inf_name,"__max_filtered_plot_500.pdf",sep = ""),width=14)
plot(data_out_filter_500[,1]~data_out_filter_500[,2],data_out_filter_500,type="l")
dev.off()
pdf(file=paste(inf_name,"__max_filtered_plot_1000.pdf",sep = ""),width=14)
plot(data_out_filter_1000[,1]~data_out_filter_1000[,2],data_out_filter_1000,type="l")
dev.off()
pdf(file=paste(inf_name,"__max_filtered_plot_1500.pdf",sep = ""),width=14)
plot(data_out_filter_1500[,1]~data_out_filter_1500[,2],data_out_filter_1500,type="l")
dev.off()


