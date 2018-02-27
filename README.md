# RNA-Seq data analysis in cases with Complex structural variants (cSVs)
These tools provide way to analyse differential gene expresion and allowing rational filtering and combinig results of fusion genes analysis with special address to cSVs cases.

  1) perform differential gene expresion in cohort without appropriate control group
  2) combine and filter fusion genes pipeline results using metacaller approach
  
## Analysis 1
Uncoventional statistical approach for differential gene expresion based on comparing individual cases in sample set without control group (for datasets, where in principle no adequate control group exist, for example cSVs as chromothripsis or chromoanasynthesis). More detail will come soon. (manuscript in preparation)

### Run analysis 1
No need for instalation
  1) Download bla.R RScript from repository 
  2) Copy inputs to bla.R RScript folder:
     <br /> Only neccesary input file is normalized expression table on log2 scale in tab separetated format          expression_table.tsv)
  3) Run bla.R RScript with appropriate parametres:
     <br /> Rscript --vanilla bla.R {expression_table.tsv} "p_value" {outputfile_name.txt}
     <br />   where: p_value is required value of sigificance for differentialy expressed gene (default: 0,05)         

Output file is in .tsv format with lines represented all differential expresed genes in dataset, and with columns for individual samples (up/down tag is mentioned)
 
## Analysis 2
uvod
how to run - inut files a parameters
outputs description
