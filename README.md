# RNA-Seq data analysis in cases with Complex structural variants (cSVs)
These tools provide way to analyse differential gene expresion and allowing rational filtering and combinig results of fusion genes analysis with special address to cSVs cases.

  1) perform differential gene expresion in cohort without appropriate control group
  2) combine and filter fusion genes pipelines results using metacaller approach
  
## Analysis 1
Uncoventional statistical approach for differential gene expresion based on comparing individual cases in sample set without control group (for datasets, where in principle no adequate control group exist, for example cSVs as chromothripsis or chromoanasynthesis). More detail will come soon. (manuscript in preparation)

### Run analysis 1
No need for instalation
  1) Download bla.R RScript from repository, ensure that R library "data.table" is include in your R library path
  2) Copy inputs to bla.R RScript folder:
     <br /> Only neccesary input file is normalized expression table on log2 scale in tab separetated format          "expression_table.tsv")
  3) Run bla.R RScript with appropriate parametres:
     <br /> Rscript --vanilla bla.R {expression_table.tsv} "p_value" {outputfile_name.txt}
     <br /> where p_value is required value of sigificance for differentialy expressed genes (default: 0,05)         

Output file is in .tsv format with lines represented all differential expresed genes in dataset, and with columns for individual samples (up/down regulation tag is used to highlight)
 
## Analysis 2
This tool for fusion genes identification was develop with aim to rationaly filter and combine multiple results from different pipelines to maximize results reliability. Such an approach of combining results from different state of art fully automated pipelines (EricScript, JAFFA, FusionCatcher) prevent to call false positive fusion events.

### Run analysis 2
No need for instalation
   1) Download bla_met.R and metacall_wrapper.sh from repository, ensure that R libraries "stringr", "reshape", "reshape2" are include in your R library path
   2) Run metacall_wrapper.sh in this way (ensure that metacall_wrapper.sh and bla_met.R is in the same folder):
      ./metacall_wrapper.sh /PATH_TO_ERICSCRIPT_RESULTS /PATH_TO_JAFFA_RESULTS /PATH_TO_FUSIONCATCHER_RESULTS /OUT_FOLDER
   <br /> NOTE: Ensure that parameteres of script are in exactly same order as mentioned. If you will provide different numbers of parameteres script will not work properly   


