
wd<-getwd()
d<-dir()
ref<-read.csv("reference.txt" ,sep=";", header=T)
a<-grep("copraRNA2_ooi_post_filtering.r", d)
b<-grep("reference.txt", d)
d<-d[-c(a,b)]
out<-matrix(,nrow(ref),5)
colnames(out)<-c("ooi_p_filtered","ooi_p_cons_fltered","ooi_ooi_cons_p_filtered","ooi_cons_p_filtered","IntaRNA_result")
for(i in 1:length(d)){
	print(i)
	setwd(paste(wd,"/",d[i],sep=""))
	
	nofil<-read.csv("CopraRNA2_final_all_ooi.csv")
	fil1<-read.csv("CopraRNA2_final_all_ooi_filtered1.csv")
	fil2<-read.csv("CopraRNA2_final_all_ooi_filtered2.csv")
	fil3<-read.csv("CopraRNA2_final_all_ooi_ooi_consensus_filtered.csv")
	fil4<-read.csv("CopraRNA2_final_all_ooi_consensus_filtered.csv")
	e<-grep("Additional.homologs", colnames(fil1))
	

	temp<-which(ref[,1]==d[i])
	
	for(j in 1:length(temp)){
		temp2<-grep(ref[temp[j],2],apply(nofil[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],5]<-as.character(nofil[temp2,3])
		}
		
		
		temp2<-grep(ref[temp[j],2],apply(fil1[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],1]<-temp2
		}
		temp2<-grep(ref[temp[j],2], apply(fil2[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],2]<-temp2
		}
		temp2<-grep(ref[temp[j],2], apply(fil3[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],3]<-temp2
		}
		temp2<-grep(ref[temp[j],2], apply(fil4[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],4]<-temp2
		}
	}
	setwd(wd)

}

out<-cbind(ref,out)

write.table(out, "benchmark_no_stm.csv", sep=";", quote=F,row.names=F)
