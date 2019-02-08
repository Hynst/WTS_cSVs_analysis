#!/usr/bin/env Rscript

# 1/function to perform t-test and get p-values for normalize expresion dataset
# 2/function for partional genes overlap across samples
# 3/function to identify outliers (diferentialy expresse gene) in datasets without distinguishable biological groups

######### check and install dependencies ##########
libraries <- c("data.table", "optparse")
for(x in libraries){
  status <- x %in% rownames(installed.packages())
    if(status == TRUE){
      library(x, character.only=TRUE)
  
    } else {
      install.packages(x)
      library(x, character.only=TRUE)
    }
}
######## arguments #########
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="expression table name",
              metavar = "character"),
  make_option(c("-o", "--output_file"), type="character", default="dif_expr_genes.txt",
              help="output file name [default = %default]", metavar="character"),
  make_option(c("-p", "--p_value"), type="numeric", default="0.05", help="p-value [default = %default]",
              metavar = "numeric"),
  make_option(c("-t", "--tolerance_ratio"), type="numeric", default="1", 
              help="tolerance ration [default = %default]; values: 0-1", metavar = "numeric")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)
out_dir <- dirname(opt$output_file)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (--file -> expression table name)", call.=FALSE)
}

################ Load functions #####################

# function 1 (getPvalues)
getPvalues<-function(x){
  # only values without row names
  
  expr_matrix <- x[,c(2:ncol(x))]
  # count p-values (t-test) for expression value for particular gene in particular sample against rest of samples
  p_values<-lapply(1:nrow(expr_matrix), function(ex){
    
    pval<-lapply(c(1:ncol(expr_matrix)), function(zet){
      j <- c(1:ncol(expr_matrix))
      remove <- c(zet)
      #print(zet)
      k <-as.integer(paste(j [! j %in% remove], sep=""))
      #print(k)
      t_test <-t.test(expr_matrix[ex,c(k)], mu= expr_matrix[ex,zet])
      return(t_test$p.value)
      #print(t_test$p.value, header=FALSE)
    })
    return(pval)
  })
  
  col_names<-colnames(expr_matrix)
  tmp <- data.frame(matrix(unlist(p_values), nrow=nrow(x), byrow=T))
  #asi ne nutny krok
  p_valuesDF <- tmp[,c(1:length(col_names))]
  names(p_valuesDF) <- col_names 
  return(p_valuesDF)
}

#function 2 (part_overlap) 
part_overlap <- function(x, num_s = y){
    list_len <- length(x)
    tab_it <- table(unlist(x))
    tab_it_perc <- tab_it / list_len
    names(tab_it_perc[tab_it_perc >= num_s])
  }

# function 3 (recognizeOutLiers)
recognizeOutLiers <- function(x,y,z,p){ 
  
  expr_matrix <- x[,c(2:ncol(x))]
  sampleID<-colnames(x[,-1])
  dif_gene_down <- list()
  dif_gene_up <- list()
  
  for(i in c(1:ncol(expr_matrix))){
    j <- c(1:ncol(expr_matrix))
    remove <- c(i)
    k <-paste(j [! j %in% remove], sep=" ")
    
    print(sampleID[i])
    print(paste("p_value",z, sep=" "))
    outliers_list_up <- list()
    outliers_list_down <- list()
    
    for(j in k){
      #print(i)
      #print(j)
      #set variables
      j<-as.integer(j)
      e<- as.character(paste(sampleID[i])) #y
      f<- as.character(paste(sampleID[j])) #x
      nazev3 <- paste("",e,"_",f,"_outliers_up", sep="")
      nazev4 <- paste("",e,"_",f,"_outliers_down", sep="")
      
      #statistics
      fit <- lm(x[,c(paste(e))] ~ x[,c(paste(f))], x)
      preds1 <- predict(fit, interval = 'confidence', level=0.99)
      df <- data.frame(outliers=cbind(preds1[,2] <= x[,(e)] & preds1[,3] >= x[,(e)]), genes_symbol=x$Gene_symbol,
                       interval_lwr = preds1[,2], interval_upr=preds1[,3], y_value=x[,(e)],
                       ttest_pvalue=y[,(e)], UP_regulation=cbind(preds1[,3] < x[,(e)]))
      
      df_up <- df[df[,5]> 0 & df[,1] == "FALSE" & df[,6] < z & df[,7] == "TRUE",]
      df_down <- df[df[,5]> 0 & df[,1] == "FALSE" & df[,6] < z & df[,7] == "FALSE",]
      
      outliers_list_up[[nazev3]]<-as.vector(df_up$genes_symbol)
      outliers_list_down[[nazev4]]<-as.vector(df_down$genes_symbol)
      
    }
    
    # merge genes to get genes from all pair combinnations
    dif_gene_up[[paste(sampleID[i])]] <- part_overlap(outliers_list_up, as.numeric(opt$tolerance_ratio))
    dif_gene_down[[paste(sampleID[i])]] <- part_overlap(outliers_list_down, as.numeric(opt$tolerance_ratio))
    
    #save outputs (mozna zakomentovat!!!)
    lapply(dif_gene_up, function(r){
      write.table(r, file=paste0(out_dir, "/", sampleID[i],"_upregulated_genes.txt"), sep="\t", quote = F, row.names = F)
    })
    
    lapply(dif_gene_down, function(q){
      write.table(q, file=paste0(out_dir, "/", sampleID[i],"_downregulated_genes.txt"), sep="\t", quote = F, row.names = F)
    })
    
  }
  
  tab_up <- rbindlist(lapply(dif_gene_up, function(s){
    
    a<-parent.frame()$i[]
    data.table(gene=s, sample= sampleID[a], reg = "up" )
    
  }))
  
  tab_down <- rbindlist(lapply(dif_gene_down, function(q){
    
    b<-parent.frame()$i[]
    data.table(gene=q, sample= sampleID[b], reg = "down" )
    
  }))
  
  final_table_tmp<-rbind(tab_up, tab_down)
  final_table_tmp<-final_table_tmp[!duplicated(final_table_tmp$gene),]
  
  # write final results
  final_table<-data.table::dcast.data.table(final_table_tmp, formula=gene~sample, value.var = "reg", fill = "0" )
  write.table(final_table, file=p, sep="\t", quote = F, row.names = F)
  
  
}

######## Run pipeline ##########

print("Loading expression table..")
eset<-read.table(opt$file, header=TRUE)
dimen<-dim(eset)
print(paste("Number of genes:", dimen[1], sep=" "))
print(paste("Number of samples:", dimen[2]-1, sep=" "))
print(paste("P-value:", opt$p_value, sep=" "))
print(paste("Tolerance ratio:", opt$tolerance_ratio, sep=" "))

#set directory
cd <- getwd()
setwd(cd)

#execute functions
print("Computing p-values (t-test)..")
p_valuesDF<-getPvalues(eset)
print("Outliers recognition in..")
recognizeOutLiers(eset, p_valuesDF, opt$p_value, opt$output_file)

