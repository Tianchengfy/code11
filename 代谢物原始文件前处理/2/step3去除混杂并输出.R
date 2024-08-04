library(pacman)
library(TwoSampleMR)
setwd("E:/代谢物最新/去混杂")
p_load(data.table,purrr,readr,MRInstruments,dplyr,TwoSampleMR)

a<-getwd()
a<-paste0(getwd(),"/clumped/")
a
files<-list.files(path=a,pattern = "^GCST.*\\.txt",full.names = T)
files
if(!dir.exists("by")){
  dir.create("by")}
phdata <- fread("E:/代谢物最新/去混杂/clumped/out/bb.csv")
confounders<-c("Body mass index","Smoking status")
pattern<-paste0("^(",paste(confounders,collapse="|"),")$")
remove_snps<-phdata%>%
  filter(grepl(pattern,trait))
remove_snps

get_print_filepath<-function(files){
  filepath<-paste0("./by/",basename(files))
  cat("File path:",filepath,"\n")
  return(filepath)}
get_print_filepath
for(i in seq_along(files)){
  expo_data<-read_exposure_data(filename=files[i],sep = "\t",snp_col = "SNP",  beta_col = "beta.exposure",  
                                se_col = "se.exposure", 
                                phenotype_col="exposure",
                                id_col = "processed_name",
                                effect_allele_col = "effect_allele.exposure",  
                                other_allele_col = "other_allele.exposure",  
                                eaf_col = "eaf.exposure",  
                                pval_col = "pval.exposure", 
                                samplesize_col = "samplesize.exposure",   
  )
  if(nrow(remove_snps)>0){
    remove_snps$phenotype<-files[i]
    expo_data<-expo_data[!(expo_data$SNP%in%remove_snps$snp),]
    
  }
  
  filepath<-get_print_filepath(files[i])
  fwrite(expo_data,filepath,na="",quote=F,row.name=F,sep = "\t")
}
