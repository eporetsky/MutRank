# Create an empty default_files.csv file if it doesn't exist
if(!file.exists("default_files.csv")){
  default_files <- data.frame("file_name"=c(rep(NA,8)),row.names=c("annotations","categories","data","domains","foldchange","GO_db","GO_genes","symbols"))
  write.table(default_files,"default_files.csv",sep=",", col.names=NA, quote=F)
}
# Read the default_files.csv that specifies which files to load when MutRank starts
defaults <- read.table("default_files.csv",row.names=1,header=T,sep=",",colClasses = "character")