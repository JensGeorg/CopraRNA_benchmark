

# Benchmark Melamed et al. results
benchfile<-"reference.txt"
benchdat<-read.csv(benchfile, sep="\t")
rildata<-read.csv("melamed_results.txt",sep="\t")
nor<-cbind(rildata[,c(6,8,10)])

mnor<-c()
for(i in 1:nrow(nor)){
	mnor<-c(mnor,max(nor[i,],na.rm = TRUE))
}

rildata<-cbind(rildata,mnor)
rildata<-rildata[order(rildata[,"mnor"], decreasing=T),]



srna<-as.character(na.omit(unique(benchdat[,1])))

positive_list<-rep(NA,nrow(benchdat))
for(i in 1:length(srna)){
	print(i)
	targets<-which(benchdat[,1]==srna[i])
	if(srna[i]=="Spot42"){
		srna[i]<-"spf"
	}
	h<-paste("^",srna[i],sep="")
	temp<-grep(h,rildata[,1],ignore.case=T)
	temp2<-grep(h,rildata[,2],ignore.case=T)
	temp_table<-matrix(,length(temp)+length(temp2),4)
	temp_table[,1]<-c(as.character(rildata[temp,1]),as.character(rildata[temp2,2]))
	temp_table[,2]<-c(as.character(rildata[temp,2]),as.character(rildata[temp2,1]))
	temp_table[,3]<-c(as.character(rildata[temp,"mnor"]),as.character(rildata[temp2,"mnor"]))
	temp_table[,4]<-c(as.character(rildata[temp,4]),as.character(rildata[temp2,3]))
	temp_table<-temp_table[order(as.numeric(temp_table[,3]),decreasing=T),]
	if(length(temp_table)>0){
		fil<-which(temp_table[,4]=="AS" | temp_table[,4]=="other-ncRNA" | temp_table[,4]=="sRNA" | temp_table[,4]=="tRNA" | temp_table[,4]=="cis_AS_with_trans_t")
		if(length(fil)>0){
			temp_table<-temp_table[-fil,]
		}
		if(is(temp_table)[2]=="vector"){
			temp_table<-matrix(temp_table,1,4)
		}
		
		for(j in 1:length(targets)){
			temp_target<-strsplit(as.character(benchdat[targets[j],3]),"/")[[1]]
			for(ii in 1:length(temp_target)){
				temp_pos<-grep(temp_target[ii],temp_table[,2],ignore.case=T)
				print(temp_pos)
				if(length(temp_pos)>0){
					positive_list[targets[j]]<-temp_pos[1]
				}
			}
		}
	}
}
interactome_melamed<-positive_list

benchdat<-cbind(benchdat,interactome_melamed)

write.table(benchdat,file="benchmark_rilseq.csv",sep=";")
