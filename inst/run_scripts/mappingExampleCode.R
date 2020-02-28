library(tidyverse)
library(fgcz.gsea.ora)
library(httr)

# more information

# https://www.uniprot.org/help/uploadlists
# and here

fc_estimates <- readxl::read_xlsx("results_modelling_testing/modelling_results_peptide/foldchange_estimates.xlsx")


mapping <- map_ids_uniprot( fc_estimates, ID_col = "top_protein")
#head(mapping)
#mean(is.na(mapping$P_ENTREZGENEID))
writexl::write_xlsx(mapping, "results_modelling_testing/modelling_results_peptide/foldchange_estimates_Entrez.xlsx")



#Rscript C:\Users\wolski\prog\fgczgseaora\inst\run_scripts\lfq_multigroup_gsea.R .\results_modelling_intest\modelling_results_peptide\foldchange_estimates_entrez.xlsx -o sscrofa -i P_ENTREZGENEID -t entrezgene
#Rscript C:\Users\wolski\prog\fgczgseaora\inst\run_scripts\lfq_multigroup_ora.R .\results_modelling_intest\modelling_results_peptide\foldchange_estimates_entrez.xlsx -o sscrofa -i P_ENTREZGENEID -t entrezgene

#Rscript C:\Users\wolski\prog\fgczgseaora\inst\run_scripts\lfq_multigroup_gsea.R .\results_modelling_liver\modelling_results_peptide\foldchange_estimates_entrez.xlsx -o sscrofa -i P_ENTREZGENEID -t entrezgene
#Rscript C:\Users\wolski\prog\fgczgseaora\inst\run_scripts\lfq_multigroup_ora.R .\results_modelling_liver\modelling_results_peptide\foldchange_estimates_entrez.xlsx -o sscrofa -i P_ENTREZGENEID -t entrezgene
