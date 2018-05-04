


ref<-read.csv("../../reference_04_05_2018.txt" ,sep="\t", header=T)
d<-dir()

sRNA<-gsub("\\.txt","",d)


targetRNA2<-c()
st<-function(x){
		temp<-strsplit(x,"\t")
		temp<-temp[[1]]
		temp
	}

for(i in 1:length(d)){
	npos<-which(as.character(ref[,1])==as.character(sRNA[i]))
	test<-read.delim(d[i], skip=20, sep="!")
	test<-test[,1]
	en<-grep("\\%", test)[1]
	test<-test[1:(en-1)]
	test<-as.character(test)
	test<-gsub("\t\t","\t",test)
	temp2<-lapply(test,st)
	temp2<-matrix(unlist(temp2),length(temp2),9,byrow=T)
	
	
	for(ii in 1:length(npos)){
		
		temp3<-grep(gsub("/.*","",as.character(ref[npos[ii],2])),temp2[,3])[1]
		targetRNA2<-c(targetRNA2, temp3)
		
	}
	
	
}
write.table(out,file="targetRNA2.txt", sep="\t")
