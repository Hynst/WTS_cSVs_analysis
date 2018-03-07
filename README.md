# RNA-Seq data analysis in cases with Complex structural variants (cSVs) - REPOSITORY IN PREPARATION!!!
These tools provide way to analyse differential gene expresion and allowing rational filtering and combinig results of fusion genes analysis with special address to cSVs cases. 

  1) perform differential gene expresion in cohort without appropriate control group
  2) combine and filter fusion genes pipelines results using metacaller approach

More detail about both methods will come soon. (manuscript in preparation)
  
## Analysis 1
Uncoventional statistical approach for differential gene expresion based on comparing individual cases in sample set without control group (for datasets, where in principle no adequate control group exist, for example cSVs as chromothripsis or chromoanasynthesis).

### Run analysis 1
No need for instalation
  1) Download run_DE.R RScript from repository, ensure that R library "data.table" is include in your R library path
  2) Copy inputs to bla.R RScript folder:
     <br /> Only neccesary input file is normalized expression table on log2 scale in tab separetated format          "expression_table.tsv")
  3) Run run_DE.R RScript with appropriate parametres:
     <br /> `Rscript --vanilla run_DE.R {expression_table.tsv} "p_value" {outputfile_name.txt}`
     <br /> where p_value is required value of sigificance for differentialy expressed genes (default: 0,05)         

Output file is in .tsv format with lines represented all differential expresed genes in dataset, and with columns for individual samples (up/down regulation tag is used to highlight)
 
## Analysis 2
This tool for fusion genes identification was develop with aim to rationaly filter and combine multiple results from different pipelines to maximize results reliability. Such an approach of combining results from different state of art fully automated pipelines (EricScript, JAFFA, FusionCatcher) prevent to call false positive fusion events.

### Run analysis 2
No need for instalation
   1) Download filter_results.R and run_MCaller_wrapper.sh from repository, make run_MCaller_wrapper.sh executable `chmod 744 run_MCaller_wrapper.sh`
   2) Copy all samples results from each pipelines to new folders ie. ~/eric , ~/fc , ~/jaffa
   3) Please name files prefix in each folder consistently ie. S1_eric.results.tsv, S1_jaffa.results.tsv,     S1_fc.results.tsv
   <br /> NOTE: here "S1" is recognize as sample name, _ is important for sample name recognition. Please provide _ inmediately after sample name
   <br /> NOTE: final data folder structure should look like:
<br />.
<br />├── eric
<br />│   └── S1_eric.results.tsv
<br />├── fc
<br />│   └── S1_fc.results.tsv
<br />└── jaffa
<br />    └── S1_jaffa.results.tsv
<br />
   4) Run run_MCaller_wrapper.sh in this way (ensure that metacall_wrapper.sh and filter_results.R is in the same folder):
```
./run_MCaller_wrapper.sh \
/PATH_TO_ERICSCRIPT_RESULTS \
/PATH_TO_JAFFA_RESULTS \
/PATH_TO_FUSIONCATCHER_RESULTS \
/OUT_FOLDER
```
      
   NOTE: Ensure that parameteres of script are in exactly same order as mentioned. If you will provide different numbers of parameteres script will not work properly   


