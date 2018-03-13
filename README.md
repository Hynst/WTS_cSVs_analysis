# RNA-Seq data analysis with Complex structural variants (cSVs)
These tools provide way to inspect differential gene expresion and allow rational filtering and combinig results of fusion gene events with special address to cSVs cases. 

  <br /> Analysis 1) is performing differential gene expresion in cohort without appropriate control group
  <br /> Analysis 2) is combining and filtering fusion genes pipelines results using metacaller approach

More detail about both methods will come soon. (manuscript in preparation)
  
## Analysis 1
Differential expression analysis for cSVs datasets, for example datasets with chromothripsis or chromoanasynthesis, where in principle no adequate control group exist.

### Run analysis 1
  1) Download run_DE.R RScript from repository
  2) Copy inputs to run_DE.R RScript folder:
     <br /> Input file is normalized expression table on log2 scale in tab separetated format (.tsv)         
  3) Run run_DE.R RScript with parametres:
     <br /> `Rscript --vanilla run_DE.R {input_expr_table.tsv} "p_value" {outputfile_name.txt}`
     <br />  p-value default: 0.05         

Output file is in tab separated format with lines representing all differential expresed genes in dataset, and with columns for individual samples (up/down regulation tag is used)
 
## Analysis 2
This tool for fusion genes identification was develop with aim to rationaly filter and combine multiple results from different pipelines (EricScript, JAFFA, FusionCatcher) to maximize results reliability and to prevent call of false positive fusion events.

### Run analysis 2
   1) Download filter_results.R and run_MCaller_wrapper.sh from repository, make run_MCaller_wrapper.sh executable `chmod 744 run_MCaller_wrapper.sh`
   2) Copy all samples results from each pipelines to new folders ie. ~/PATH/eric , ~/PATH/fc , ~/PATH/jaffa
   3) Please name files prefix in each folder consistently ie. S1_eric.results.tsv, S1_jaffa.results.tsv,     S1_fc.results.tsv
   <br /> NOTE1: here "S1" is recognized as sample name, _ is important for sample name identification in string. Please provide _ inmediately after sample name
   <br /> NOTE2: final data folder structure should look like:
   ```
.
├── eric
│   └── S1_eric.results.tsv
├── fc
│   └── S1_fc.results.tsv
└── jaffa
└── S1_jaffa.results.tsv
   ```
   4) Run run_MCaller_wrapper.sh in this way (ensure that metacall_wrapper.sh and filter_results.R are in the same folder):

   ```
   ./run_MCaller_wrapper.sh \
   /PATH_TO_ERICSCRIPT_RESULTS \
   /PATH_TO_JAFFA_RESULTS \
   /PATH_TO_FUSIONCATCHER_RESULTS \
   /OUT_FOLDER
   ```
      
   NOTE: Ensure that parameteres of script are in exactly same order as mentioned. If you will provide different numbers of parameteres script will not work properly.   


