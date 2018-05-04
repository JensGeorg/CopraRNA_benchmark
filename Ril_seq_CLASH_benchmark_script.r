
# Benchmark Melamed et al. results
benchfile<-"reference_03_05_2018.txt"
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
	
	if(srna[i]=="OmrAB"){
		#srna[i]<-"omrA"
		h<-paste("^","omrA",sep="")
		temp<-grep(h,rildata[,1],ignore.case=T)
		#h<-paste("^","omrA",sep="")
		temp2<-grep(h,rildata[,2],ignore.case=T)
		#srna[i]<-"omrB"
		h<-paste("^","omrB",sep="")
		temp<-c(temp,grep(h,rildata[,1],ignore.case=T))
		temp2<-c(temp2,grep(h,rildata[,2],ignore.case=T))
		temp_table<-matrix(,length(temp)+length(temp2),4)
		temp_table[,1]<-c(as.character(rildata[temp,1]),as.character(rildata[temp2,2]))
		temp_table[,2]<-c(as.character(rildata[temp,2]),as.character(rildata[temp2,1]))
		temp_table[,3]<-c(as.character(rildata[temp,"mnor"]),as.character(rildata[temp2,"mnor"]))
		temp_table[,4]<-c(as.character(rildata[temp,4]),as.character(rildata[temp2,3]))
		temp_table<-temp_table[order(as.numeric(temp_table[,3]),decreasing=T),]
	}
	
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


################

srna<-as.character(na.omit(unique(benchdat[,1])))
clash<-read.csv("clash_rnase_sRNA_mRNA.txt",sep="\t", comment.char="#")
clash<-as.matrix(clash)


s1<-grep("sRNA", clash[,2])
s2<-grep("sRNA", clash[,5])


clashtab<-matrix(,length(s1)+length(s2),12)

for(i in 1:length(s1)){
#print(i)
	clashtab[i,1]<-clash[s1[i],2]
	clashtab[i,2]<-clash[s1[i],5]
	clashtab[i,3:12]<-clash[s1[i],8:17]
}


for(i in 1:length(s2)){
	
	clashtab[i+length(s1),1]<-clash[s2[i],5]
	clashtab[i+length(s1),2]<-clash[s2[i],2]
	clashtab[i+length(s1),3:12]<-clash[s1[i],8:17]
}


for(i in 1:nrow(clashtab)){
	
	clashtab[i,1]<-strsplit(clashtab[i,1],"\\|",)[[1]][2]
	clashtab[i,2]<-strsplit(clashtab[i,2],"\\|",)[[1]][2]
}

ryeb<-which(clashtab[,1]=="ryeB")
clashtab[ryeb,1]<-"SdsR"

spf<-which(clashtab[,1]=="spf")
clashtab[spf,1]<-"Spot42"

sc<-as.numeric(clashtab[,ncol(clashtab)])
clashtab<-clashtab[order(sc),]

positive_list<-rep(NA,nrow(benchdat))
for(i in 1:length(srna)){
	print(i)
	targets<-which(benchdat[,1]==srna[i])
	
	
	temp<-grep(srna[i],clashtab[,1],ignore.case=T)
	if(srna[i]=="OmrAB"){
		
		temp<-grep("omrA",clashtab[,1],ignore.case=T)
		
		temp<-c(temp,grep("omrB",clashtab[,1],ignore.case=T))
	}
	temp_table<-clashtab[temp,]
	if(length(temp_table)>0){
		
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

interactome_clash<-positive_list
benchdat<-cbind(benchdat,interactome_clash)

write.table(benchdat,file="benchmark_rilseq_clash.csv",sep=";", quote=F, row.names=F)
