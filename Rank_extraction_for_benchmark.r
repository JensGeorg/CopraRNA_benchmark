



# extracting ranks from the result files


wd<-getwd()
d<-dir()
ref<-read.csv("reference.txt" ,sep=";", header=T)
filelist<-c("CopraRNA2_final_all_ooi.csv","CopraRNA2_final_all_ooi_consensus.csv","CopraRNA2_final_all_ooi_ooiconsensus.csv","CopraRNA2_final_all_balanced.csv","CopraRNA2_final_all_balanced_consensus.csv")
filelist<-c(filelist,filelist)
out<-matrix(,nrow(ref),length(filelist))

sRNA<-unique(ref[,1])

for(i in 1:length(sRNA)){
	print(sRNA[i])
	if(sRNA[i]!="ChiX"){
	temp<-match(sRNA[i],d)
	npos<-which(as.character(ref[,1])==as.character(sRNA[i]))
	setwd(paste(wd,"/",d[temp],sep=""))
	#system("R --slave -f  ../conservation_test.r --args thres=0.32 thres2=0.51 thres_overlap=0.45")
	load("heatmap_data.Rdata")
	out_table2<-res[[1]]
	for(j in 1:(length(filelist)/2)){
		dat<-read.csv(filelist[j])
		neg<-which(out_table2[,1]>4)
		dat2<-dat
		
			if(length(neg)>0){
				dat2<-dat2[order(dat2[,"initial_sorting"]),]
				
				dat2<-dat2[-neg,]
				dat2<-dat2[order(dat2[,"p.value"]),]
				
			}
		
		
		for(ii in 1:length(npos)){
			temp2<-grep((ref[npos[ii],2]),(dat[,3]))[1]
			
			temp3<-grep(as.character(ref[npos[ii],2]),(dat2[,3]))[1]
			
			print(temp3)
			if(is.na(temp2)==F){
				out[npos[ii],j]<-temp2
			}
			if(is.na(temp3)==F){
				out[npos[ii],(j+length(filelist)/2)]<-temp3
			}
		}
	}
	print(out[npos,])
	setwd(wd)
	}
}

out<-cbind(ref,out)



write.table(out,file="test.txt", sep="\t")
