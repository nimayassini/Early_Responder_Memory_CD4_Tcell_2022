## Early responder memory CD4+ T cells modulate heterologous disease outcome

In this repository, the R scripts used to analyze the bulk-RNAseq data from the project __"Early responder memory CD4+ T cells modulate heterologous disease outcome"__ can be found.

The main code is in the script [bulkRNAseq_analysis_201021.R](https://github.com/nimayassini/Early_Responder_Memory_CD4_Tcell_2022/blob/main/scripts/bulkRNAseq_analysis_201021.R), which also calls the R scripts
- [0_upset_deg.R](https://github.com/nimayassini/Early_Responder_Memory_CD4_Tcell_2022/blob/main/scripts/0_upset_deg.R)
- [0_upset_deg_updown.R](https://github.com/nimayassini/Early_Responder_Memory_CD4_Tcell_2022/blob/main/scripts/0_upset_deg_updown.R)
- [0_upset_deg_updown_acttax.R](https://github.com/nimayassini/Early_Responder_Memory_CD4_Tcell_2022/blob/main/scripts/0_upset_deg_updown_acttax.R)
to generate UpSet plots.

The csv file "setup.csv" is used to retrieve information about the samples.
