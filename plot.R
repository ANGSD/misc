

require(plotrix)
library(plotrix)


########### do not change ################3
l<-commandArgs(TRUE)
getArgs<-function(x,l)
  unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
Args<-function(l,args){
 if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
  cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument")
  q("no")
}
 arguments<-list()
 for(a in names(args))
   arguments[[a]]<-getArgs(a,l)

 if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
   cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
   q("no")
 }
 for(a in names(args))
   if(is.null(arguments[[a]]))
     arguments[[a]]<-args[[match(a,names(args))]]

   
 arguments
}
 
print.args<-function(args,des){
  if(missing(des)){
    des<-as.list(rep("",length(args)))
    names(des)<-names(args)
  }
  cat("->  needed arguments:\n")
  mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
  cat("->  optional arguments (defaults):\n")
  mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
  q("no")
}
###### ####### ###### ###### ###### #######
# choose your parameters and defaults
# NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments
args<-list(ifile=NULL,
           out="insertstat",
           max_insert=5e3,
           main="Error rate using an outgroup and a high quality genome",
           height=7,
           width=11,
           srt=90,
           cex=1,
           doPng=FALSE
           )
#if no argument aree given prints the need arguments and the optional ones with default
des<-list(ifile="the tmpfile",
          out="Name of the outfile",
          max_insert="Cutoff for max insertsize",
          width="width of the pdf",
          height="height of the pdf",
          srt="angle of ind names",
          cex="scale of names",
          doPng="Make png instead of pdf (not implemented)"
          )
######################################
#######get arguments and add to workspace
### do not change
if(length(l)==0) print.args(args,des)
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
  cat(" Arguments: output prefix\n")
  q("no")
}


con <- file(ifile)
d <-  read.table(file=con,nrows=-1,as.is=T)
misc <- d[1:6,]
d <- d[-c(1:6),]
m<-apply(d,2,as.numeric)


p <- which(diff(m[,1])<0)
insert <- m[1:p,]
rlen <- m[-c(1:p),]
above_thres <- sum(insert[,2][insert[,1]>max_insert])
insert <- insert[insert[,1]<max_insert,]

misc <- rbind(misc,c("nInserts>5k:\t",above_thres))

onam <-paste0(out,".pdf")
cat("output pdf is: ",onam,"\n")
pdf(onam,width=width,height=height)
par(mfrow=c(1,2))
barplot(rlen[,2],names=rlen[,1],xlab="Read lengths",ylab="Counts")

barplot(insert[,2],names=insert[,1],xlab="Insert Size",ylab="Counts")
addtable2plot(x=nrow(insert)/2,y=max(insert[,2])/2,misc)
dev.off()
