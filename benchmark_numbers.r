
wd<-getwd()
d<-list.dirs( recursive=FALSE)
ref<-read.csv("reference_04_05_2018.txt" ,sep="\t", header=T)
IntaRNA<-c()
ooi_org<-"NC_000913"



out<-matrix(,nrow(ref),10)
colnames(out)<-c("c2_ooi","c2_ooi_cons","c2_ooi_ooicons","c2_bal","c2_bal_cons","ooi_filtered","ooi_ooi_cons_p_filtered","bal_filtered","bal_cons_filtered","IntaRNA")
for(i in 1:length(d)){
	print(i)
	setwd(paste(wd,"/",d[i],sep=""))
	
	
	
	ooi<-read.csv("CopraRNA2_final_all_ooi.csv")
	ooi_cons<-read.csv("CopraRNA2_final_all_ooi_consensus.csv")
	ooi_ooi_cons<-read.csv("CopraRNA2_final_ooi_ooiconsensus.csv")
	bal<-read.csv("CopraRNA2_final_all_balanced.csv")
	bal_cons<-read.csv("CopraRNA2_final_all_balanced_consensus.csv")
	
	fil1<-read.csv("ooi_filtered_filtered.csv")
	fil2<-read.csv("ooi_cons_filtered_filtered.csv")
	fil3<-read.csv("all_filtered_filtered.csv")
	fil4<-read.csv("all_cons_filtered_filtered.csv")
	fil<-paste(ooi_org, "_upfromstartpos_200_down_100.fa.intarna.sorted.csv",sep="")
	datInt<-read.csv(fil)
	e<-grep("Additional.homologs", colnames(fil1))
	

	temp<-which(ref[,1]==gsub("\\./","",d[i]))
	
	for(j in 1:length(temp)){
		temp2<-grep(ref[temp[j],2],apply(ooi[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],1]<-temp2
		}
		temp2<-grep(ref[temp[j],2],apply(ooi_cons[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],2]<-temp2
		}
		temp2<-grep(ref[temp[j],2],apply(ooi_ooi_cons[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],3]<-temp2
		}
		temp2<-grep(ref[temp[j],2],apply(bal[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],4]<-temp2
		}
		temp2<-grep(ref[temp[j],2],apply(bal_cons[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],5]<-temp2
		}
		temp2<-grep(ref[temp[j],2],apply(fil1[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],6]<-temp2
		}
		temp2<-grep(ref[temp[j],2], apply(fil2[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],7]<-temp2
		}
		temp2<-grep(ref[temp[j],2], apply(fil3[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],8]<-temp2
		}
		temp2<-grep(ref[temp[j],2], apply(fil4[,3:e], 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],9]<-temp2
		}
		
		temp2<-grep(ref[temp[j],2], apply(datInt, 1, paste, collapse=""))[1]
		if(length(temp2)>0){
			out[temp[j],10]<-temp2
		}
		
	}
	setwd(wd)

}

out<-cbind(ref[,1:3],out[,1:10],ref[,"reference"])

write.table(out, "benchmark_no_stm.csv", sep=";", quote=F,row.names=F)

