pSILAC Atlas v38 update (test locally first)

1. BACK UP your current app/app.R and app/data/shiny_data directory.
2. Copy app.R from this ZIP into <repo>/app/app.R.
3. Copy data/shiny_data/phosphosite_annotations_v38_CORRECTED.csv into <repo>/app/data/shiny_data/.
4. Keep your existing RDS data, other annotation files, and www assets unchanged.
5. In RStudio, set the working directory to <repo>/app and run shiny::runApp('.').
6. Optional graphical abstract: convert your original TIFF to PNG and place at
   <repo>/app/www/pSILAC_embryo_atlas.png.

IMPORTANT: R is not installed in the build environment, so the script has NOT
been parsed or run against your local RDS files. This is a candidate v38 build,
not a validated deployment. Test all tabs and CSV exports before committing.

The 95 sites with no stage are shown as 'Not assigned'.
The phosphosite annotation file contains mouse-only annotations supplied by you.
