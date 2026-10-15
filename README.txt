pSILAC Atlas v37 – visual and wording update

1. Back up your working psilac_shinylive_test/app/app.R.
2. Copy app.R from this package to psilac_shinylive_test/app/app.R.
3. Copy www/JVOlogo_large.png to psilac_shinylive_test/app/www/JVOlogo_large.png.
4. Keep the existing app/data directory unchanged.
5. Test with shiny::runApp(app_dir).
6. Re-export with shinylive::export(appdir=app_dir,destdir=export_dir).
7. Serve with httpuv::runStaticServer(export_dir).

Ethics text is not added yet because the full ethics statement could not be verified from the linked preprint. Please supply the exact text for faithful inclusion.
The minimum-replicate limit is computed from the maximum observed count in the 5-day daily summaries.
