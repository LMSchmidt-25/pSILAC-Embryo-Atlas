pSILAC Embryo Atlas v39 — integration candidate (not R-runtime-tested)

1. Keep your existing v38 directory unchanged.
2. Extract this ZIP into a NEW directory.
3. REQUIRED: Copy from your existing working v38 app into the new directory:
   data/shiny_data/protein_annotation.rds
   data/shiny_data/protein_localization.rds (if available)
   These files were not among the uploaded files, and are NOT included here.
4. Open the new directory as the R working directory and run shiny::runApp().
5. Validate startup, site counts, source filter, and a few known proteins before deployment.

This ZIP includes the uploaded enrichment and combined annotation RDS files and the original proteome and phosphoproteome data available in the conversation.
The v39 code checks for 12,002 combined annotation rows and 9,634 enrichment rows and tests unique site keys and mapping.
Enrichment fold changes and adjusted statistics are not SILAC incorporation.
The original v38 app is not modified.

IMPORTANT: The execution environment does not have R installed. R syntax, Shiny startup, and biological results have not been runtime-tested. This is a draft for local validation, not a verified deployment.


v39 readability/search revision (2026-11-05):
- Phosphosite table retains the eight v38 columns plus Measurement source; enrichment statistics are not displayed.
- Enrichment-only sites have NA maximum incorporation, never imputed from enrichment.
- 24 h/48 h boxplots and 5-day plots use larger graphics and labels.
- Protein search selects and loads the first ranked match automatically.
- No underlying dataset files or measurement calculations were changed.
