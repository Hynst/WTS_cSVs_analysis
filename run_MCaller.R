#!/usr/bin/env Rscript

print("Running Metacaller..")
options(warn=-1)

# check libraries
libraries <- c("stringr", "reshape2", "reshape", "tools", "optparse")

for(j in libraries){
  #print(j)
  status<-j %in% rownames(installed.packages())
  if(status == TRUE){
    suppressMessages(library(j, character.only=TRUE))
  } else{
    install.packages(j)
    suppressMessages(library(j, character.only=TRUE))
  }    
}

# arguments
option_list <- list(
  make_option(c("-E", "--folder_eric"), type="character", default=NULL, help="EricScript folder",
              metavar = "character"),
  make_option(c("-J", "--folder_jaffa"), type="character", default=NULL,
              help="JAFFA folder", metavar="character"),
  make_option(c("-F", "--folder_fc"), type="character", default=NULL, help="FusionCatcher folder",
              metavar = "character"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL, 
              help="out folder", metavar = "character")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

if (length(opt) <= 4){
  print_help(opt_parser)
  stop("All arguments must be provided", call.=FALSE)
}

# results data
eric_dir <- opt$folder_eric
jaffa_dir <- opt$folder_jaffa
fc_dir <- opt$folder_fc
out_folder <- opt$out_folder

# create sample sheet
eric_samplesheet <- list.files(eric_dir)
jaffa_samplesheet <- list.files(jaffa_dir)
fc_samplesheet <- list.files(fc_dir)

setwd(out_folder)

for(i in 1:length(eric_samplesheet)){
  
  ### EricScrip filtering and table rearanging ###
  
  eric_results <- read.table(paste(eric_dir, "/", eric_samplesheet[i], sep=""), sep="\t", quote = "")
  eric_results <- eric_results[-1,]
  # filtering to only high confidence fusions ES > 0.90
  eric_results <- eric_results[as.numeric(as.character(eric_results$V25)) > 0.90,]
  eric_results <- eric_results[-grep("Unable to predict breakpoint position", eric_results[,4]),]
  eric_results <- eric_results[-grep("Unable to predict breakpoint position", eric_results[,7]),]
  
  #Extract sample name
  sample_name<-basename(paste0(eric_samplesheet[i]))
  sample_string<-gsub("_.*","", sample_name)
  
  ### JAFFA filtering and table rearanging ###
  jaffa_results <- read.table(paste(jaffa_dir, "/", jaffa_samplesheet[i], sep=""), sep="\t", quote = "" )
  # filtering to only high confidence fusions -
  jaffa_results <- jaffa_results[grep("HighConfidence",jaffa_results[,17]),]
  
  ### FusionCatcher filtering a table rearanging ###
  fusioncatcher_results <- read.table(paste(fc_dir, "/", fc_samplesheet[i], sep=""), sep="\t", quote = "")
  
  caller1 <- data.frame(sample=sample_string, gene1=eric_results[,1], position1 = eric_results[,4],
                        chromosome1=eric_results[,3],
                        gene2=eric_results[,2], position2 = eric_results[,7],
                        chromosome2=eric_results[,6],
                        tool="eric")
  caller1 <-  caller1[-1,]    
  caller1 <-caller1[!duplicated(caller1[,c(2,4,5,7)]),]
  
  caller2 <- data.frame(sample=sample_string, gene1=fusioncatcher_results[,1], position1 = gsub(":", "", str_sub(fusioncatcher_results[,9],3)),
                        chromosome1=gsub(":", "",str_sub(fusioncatcher_results[,9],1,2)),
                        gene2=fusioncatcher_results[,2], position2 = gsub(":", "", str_sub(fusioncatcher_results[,10],3)),
                        chromosome2=gsub(":", "",str_sub(fusioncatcher_results[,10],1,2)),
                        tool="fusioncatcher")
  #caller2<-caller2[gsub("+", "", gsub("-", "", caller2$position1)),]
  #caller2<-caller2[gsub("+", "", gsub("-", "", caller2$position1)),]
  caller2 <-  caller2[-1,]
  caller2 <-caller2[!duplicated(caller2[,c(2,4,5,7)]),]
  
  #split fusion genes from jaffa results to two diferent columns    
  genes<-str_split_fixed(jaffa_results[,2], ":", 2)
  gene1<-gsub("^.","",genes[,1])
  gene2<-gsub(".$","",genes[,2])
  
  caller3 <- data.frame(sample=sample_string, gene1=gene1, position1 = jaffa_results[,4],
                        chromosome1=gsub('"', "",gsub("chr","",jaffa_results[,3])),
                        gene2=gene2, position2 = jaffa_results[,7],
                        chromosome2=gsub('"', "",gsub("chr","",jaffa_results[,6])),
                        tool="jaffa")
  caller3 <-  caller3[-1,]
  caller3 <-caller3[!duplicated(caller3[,c(2,4,5,7)]),]
  
  # paste
  callers <-Reduce(function(x, y) merge(x, y, all=TRUE), list(caller1, caller2, caller3))
  header<-data.frame(sample="sample", gene1="gene1", gene1_position="gene1_pos", chromosome1="chr1", gene2="gene2", gene2_position="gene2_pos", chromosome2="chr", tool="tool")
  
  file <- paste0(out_folder, "/Sample",i,"_allmethods_results_pos.txt")
  file_out <- paste("Sample",i,"_method_overlap_pos.txt", sep="")
  
  write.table(header, file="header", sep="\t",row.names=F, quote=F)
  write.table(callers, file=file, sep="\t",row.names=F, quote=F)
  
  # metacaller
  command_awk <- paste("awk 'n=x[$2,$4,$5,$7]{print n", '"\\n\"', "$0;} {x[$2,$4,$5,$7]=$0;}'", sep = "")
  system(paste(command_awk, file, ">", file_out, sep=" "))

}

system(paste("head -n 1 header > header1"))
system(paste("cat header1 *method_overlap_pos.txt >> all_samples_results_metacaller.txt"))
system(paste("rm header header1 *_allmethods_results_pos.txt *method_overlap_pos.txt"))

# create final table
data_res<-read.table("all_samples_results_metacaller.txt", header = T, sep = "\t")
test<-melt(data_res, id_var=c("gene1", "gene2", "chromosome1", "chromosome2", "sample"), measure.vars = "tool")
final_df<-suppressMessages(dcast(test, sample+gene1+chromosome1+gene2+chromosome2~value, value.var= "variable"))
## change values from 2 to 1!!!

write.table(final_df, file="Metacaller_final_table.txt", sep="\t", quote=F, row.names = F)

print("..done")
