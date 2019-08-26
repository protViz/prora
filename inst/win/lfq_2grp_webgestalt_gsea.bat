@echo off
:: start help with "lfq_2grp_webgestalt_gsea --help"

SET args=%*
  Rscript %~dp0/../run_scripts/lfq_2grp_webgestalt_gsea.R %args%


