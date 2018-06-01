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
     <br /> Input file is normalized expression table on log2 scale in tab separetated format (gene names in first column) (.tsv)         
  3) Run run_DE.R RScript with parameter -h :
     <br /> `Rscript --vanilla run_DE.R -h`
     <br /> ..to view list of parameters to run a script

Output file is in tab separated format with lines representing all differential expresed genes in dataset, and with columns for individual samples (up/down regulation tag is used)
 
## Analysis 2
This tool for fusion genes identification was develop with aim to rationaly filter and combine multiple results from different pipelines (EricScript, JAFFA, FusionCatcher) to maximize results reliability and to prevent call of false positive fusion events.

### Run analysis 2
   1) Download run_MCaller.R from repository
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
   4) Run run_MCaller.R RScript with parameter -h:
   <br /> `Rscript --vanilla run_MCaller.R -h`
   <br /> ..to view list of parameters to run a script and provide all of them to run analysis
     
  
