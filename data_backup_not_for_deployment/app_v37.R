# pulsed SILAC in Embryos Atlas v1.5 — focal-protein correlation and STRING neighborhoods
# install.packages(c('shiny','dplyr','tidyr','stringr','ggplot2','DT','scales','visNetwork','igraph'))
library(shiny); library(dplyr); library(tidyr); library(stringr)
library(ggplot2); library(DT); library(scales); library(visNetwork); library(igraph)
root <- 'data'
read_data <- function(n) {
  f <- file.path(root,n)
  if(!file.exists(f)) stop('Missing: ',f)
  readRDS(f)
}
short <- read_data('atlas_pSILAC_proteome_public.rds')
long <- read_data('atlas_5day_proteome.rds')
hl <- read_data('halflife_protein.rds')
# Phosphosite tables are not exposed as separate tabs, but site counts remain
# available in the Protein Explorer overview.
phospho_files <- c('atlas_pSILAC_phosphoproteome_public.rds', 'atlas_5day_phosphoproteome.rds')
phospho_sites <- lapply(phospho_files, function(f) {
  path <- file.path(root, f)
  if(!file.exists(path)) return(NULL)
  z <- readRDS(path)
  if(!all(c('UniProtID','PhosphositeID') %in% names(z))) return(NULL)
  z %>% filter(!is.na(UniProtID), !is.na(PhosphositeID)) %>%
    distinct(UniProtID, PhosphositeID)
}) %>% bind_rows() %>% distinct(UniProtID, PhosphositeID)
phospho_counts <- phospho_sites %>% count(UniProtID, name='PhosphositeCount')
# Retain measured phosphosite rows for the Protein Explorer only.
phospho_detail <- bind_rows(lapply(phospho_files, function(f) {
  path <- file.path(root, f)
  if(!file.exists(path)) return(NULL)
  z <- readRDS(path)
  if(!all(c('UniProtID','PhosphositeID') %in% names(z))) return(NULL)
  if(!'Day' %in% names(z)) z$Day <- NA_real_
  if(!'SumIntensity' %in% names(z)) z$SumIntensity <- NA_real_
  if(!'Incorporation' %in% names(z)) z$Incorporation <- NA_real_
  if(!'Duration' %in% names(z)) z$Duration <- NA_character_
  z %>% transmute(UniProtID=as.character(UniProtID),PhosphositeID=as.character(PhosphositeID),
    Dataset=if(grepl('5day',f,fixed=TRUE)) '5-day' else as.character(Duration),
    Day=suppressWarnings(as.numeric(Day)),
    SumIntensity=suppressWarnings(as.numeric(SumIntensity)),
    Incorporation=suppressWarnings(as.numeric(Incorporation))) %>%
    filter(!is.na(UniProtID),!is.na(PhosphositeID),nzchar(PhosphositeID))
}))

# Mouse reference FASTA for site-centered +/-7 amino acid motifs.
# Set this to the location of the FASTA on your own computer if needed.
motif_candidates <- file.path(c(root, dirname(root)),
  rep(c('MOUSEc_UP000000589.fasta','MOUSEc_UP000000589(1).fasta'), each=2))
motif_candidates <- c(file.path(root,'MOUSEc_UP000000589.fasta'),
  file.path(dirname(root),'MOUSEc_UP000000589.fasta'),
  file.path(root,'MOUSEc_UP000000589(1).fasta'),
  file.path(dirname(root),'MOUSEc_UP000000589(1).fasta'))
motif_fasta <- if(any(file.exists(motif_candidates))) motif_candidates[which(file.exists(motif_candidates))[1]] else ''
read_motif_fasta <- function(path) {
  if(!file.exists(path)) return(setNames(character(), character()))
  lines <- readLines(path, warn=FALSE)
  headers <- grep('^>', lines)
  if(!length(headers)) return(setNames(character(), character()))
  ends <- c(headers[-1L]-1L,length(lines))
  accessions <- sub('^>[^|]*\\|([^|]+)\\|.*$', '\\1', lines[headers])
  seqs <- vapply(seq_along(headers), function(i)
    paste(lines[seq.int(headers[i]+1L, ends[i])], collapse=''), character(1))
  setNames(seqs[!duplicated(accessions)], accessions[!duplicated(accessions)])
}
motif_sequences <- read_motif_fasta(motif_fasta)
message('[Motifs] FASTA: ', if(nzchar(motif_fasta)) motif_fasta else 'NOT FOUND', '; sequences: ', length(motif_sequences))
site_motif <- function(site, accession) {
  # Only accept unambiguous residue-position IDs; never guess coordinates.
  site <- as.character(site)
  matches <- regmatches(site, gregexpr('(?<![A-Za-z])[STYsty]\\s*[0-9]+',site,perl=TRUE))[[1]]
  if(length(matches)!=1L) return(NA_character_)
  aa <- toupper(substr(matches,1,1)); pos <- suppressWarnings(as.integer(sub('^[STYsty]\\s*','',matches)))
  seq <- unname(motif_sequences[as.character(accession)])
  if(!length(seq) || is.na(seq) || !is.finite(pos) || pos<1 || pos>nchar(seq)) return(NA_character_)
  if(substr(seq,pos,pos)!=aa) return(NA_character_)
  positions <- seq.int(pos-7L,pos+7L)
  chars <- vapply(positions,function(j) if(j<1L || j>nchar(seq)) '-' else substr(seq,j,j),character(1))
  paste0(paste(chars[1:7],collapse=''),'[',chars[8],']',paste(chars[9:15],collapse=''))
}

ann <- readRDS(file.path(root,'shiny_data','protein_annotation.rds'))
locpath <- file.path(root,'shiny_data','protein_localization.rds')
if(file.exists(locpath)) {
  loc <- readRDS(locpath)
  ann <- ann %>% select(-any_of(c('UniProtLocation','GOCellularComponent','DeepLocLocation','DeepLocMembrane'))) %>%
    left_join(loc,by='UniProtID')
}
for(n in c('GeneSymbol','ProteinName','GeneSynonyms','UniProtLocation','GOCellularComponent','DeepLocLocation','DeepLocMembrane'))
  if(!n %in% names(ann)) ann[[n]] <- NA_character_
ann <- ann %>% distinct(UniProtID,.keep_all=TRUE) %>%
  mutate(SearchLabel=paste(coalesce(GeneSymbol,UniProtID),coalesce(ProteinName,''),paste0('[',UniProtID,']'),sep=' — '),
         SearchText=str_to_lower(paste(UniProtID,GeneSymbol,ProteinName,GeneSynonyms)))
# Contaminants/unmatched identifiers do not enter the public search.
ann <- ann %>% filter(!is.na(GeneSymbol) & nzchar(GeneSymbol))
mean_ok <- function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
sd_ok <- function(x) if(sum(is.finite(x))>1) sd(x[is.finite(x)]) else NA_real_
first_ok <- function(x) if(length(x)) min(x) else NA_real_
log2_ok <- function(x) ifelse(is.finite(x) & x>0,log2(x),NA_real_)

# ---------------------------------------------------------------------------
# Quantitative consistency QC
# Incorporation = H/(H+L), so a valid incorporation can only exist when a
# positive SumIntensity exists for the same replicate. Repair SumIntensity
# only from its defining channels when possible; never create zeros.
# ---------------------------------------------------------------------------
long <- long %>%
  mutate(
    SumIntensity_original = SumIntensity,
    SumIntensity = case_when(
      is.finite(SumIntensity) & SumIntensity > 0 ~ SumIntensity,
      is.finite(H) & is.finite(L) & (H + L) > 0 ~ H + L,
      is.finite(H) & !is.finite(L) & H > 0 ~ H,
      !is.finite(H) & is.finite(L) & L > 0 ~ L,
      TRUE ~ NA_real_
    ),
    inconsistent_incorp = is.finite(Incorporation) &
      (!is.finite(SumIntensity) | SumIntensity <= 0),
    Incorporation = ifelse(inconsistent_incorp, NA_real_, Incorporation),
    log2Intensity = log2_ok(SumIntensity)
  )

qc_repaired_sum <- sum(!is.finite(long$SumIntensity_original) & is.finite(long$SumIntensity))
qc_bad_incorp <- sum(long$inconsistent_incorp, na.rm=TRUE)
message('pSILAC app QC — SumIntensity repaired from H/L: ', qc_repaired_sum)
message('pSILAC app QC — incorporation rows without valid SumIntensity after repair: ', qc_bad_incorp)

daily <- long %>% filter(is.finite(Day)) %>% group_by(UniProtID,Day) %>%
  summarise(log2Abundance=mean_ok(log2Intensity),sdAbundance=sd_ok(log2Intensity),nAbundance=sum(is.finite(log2Intensity)),
            Incorporation=mean_ok(Incorporation),sdIncorp=sd_ok(Incorporation),nIncorp=sum(is.finite(Incorporation)),.groups='drop')
# 5 day summary metrics are only defined for proteins present in the 5 day experiment.
summary_5day <- daily %>% group_by(UniProtID) %>% summarise(
  FirstDetectedDay=first_ok(Day[nAbundance>0]),
  FirstIncorpDay=first_ok(Day[nIncorp>0]),
  DaysDetected=sum(nAbundance>0),
  MaxIncorp=if(any(is.finite(Incorporation))) max(Incorporation,na.rm=TRUE) else NA_real_,
  .groups='drop')

# Discovery is the atlas-wide protein catalogue: every annotated protein observed
# in the 24 h, 48 h, or 5 day proteome datasets is included.
short_presence <- short %>%
  mutate(
    Dataset = case_when(
      grepl('24', as.character(Duration), ignore.case=TRUE) ~ '24h',
      grepl('48', as.character(Duration), ignore.case=TRUE) ~ '48h',
      grepl("table\\s*0?5|^0?5$|si[_ ]?table\\s*0?5", str_to_lower(as.character(Table))) ~ '24h',
      grepl("table\\s*0?4|^0?4$|si[_ ]?table\\s*0?4", str_to_lower(as.character(Table))) ~ '48h',
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(UniProtID), !is.na(Dataset)) %>%
  distinct(UniProtID, Dataset)

long_presence <- long %>%
  filter(!is.na(UniProtID)) %>%
  distinct(UniProtID) %>%
  mutate(Dataset='5-day')

presence <- bind_rows(short_presence, long_presence) %>%
  distinct(UniProtID, Dataset) %>%
  group_by(UniProtID) %>%
  summarise(
    DetectedIn=paste(c('24 h','48 h','5 day')[c('24h','48h','5-day') %in% Dataset], collapse=' + '),
    .groups='drop'
  )

summary <- presence %>%
  left_join(summary_5day, by='UniProtID') %>%
  left_join(hl %>% select(any_of(c('UniProtID','half_life','k','R2','fit_status'))),by='UniProtID') %>%
  left_join(ann,by='UniProtID') %>%
  left_join(phospho_counts,by='UniProtID') %>%
  mutate(PhosphositeCount=coalesce(PhosphositeCount,0L)) %>%
  filter(!is.na(GeneSymbol) & nzchar(GeneSymbol))

# Turnover classes are relative to this dataset, not arbitrary biological cutoffs.
valid_half_lives <- summary$half_life[
  summary$fit_status == 'fit_ok' & is.finite(summary$half_life) & summary$half_life > 0
]
fast_cut <- if(length(valid_half_lives)) unname(quantile(valid_half_lives, 0.25, na.rm=TRUE)) else NA_real_
slow_cut <- if(length(valid_half_lives)) unname(quantile(valid_half_lives, 0.75, na.rm=TRUE)) else NA_real_

summary <- summary %>%
  mutate(
    TurnoverClass=case_when(
      fit_status != 'fit_ok' | !is.finite(half_life) | half_life <= 0 ~ 'Not estimable',
      is.finite(fast_cut) & half_life <= fast_cut ~ 'Fast',
      is.finite(slow_cut) & half_life >= slow_cut ~ 'Slow',
      TRUE ~ 'Intermediate'
    )
  )

days <- sort(unique(daily$Day)); days <- days[is.finite(days)]
max_replicates <- max(1L, as.integer(max(c(daily$nAbundance,daily$nIncorp), na.rm=TRUE)))
metric_choices <- c('Incorporation (%)'='Incorporation','Sum intensity (log2)'='SumIntensity',
                    'Heavy intensity (log2)'='H','Light intensity (log2)'='L','H/L ratio'='Ratio')
cols <- c('Oocyte'='#034EA2','Zygote'='#F4B978','2-cell'='#F4B978','4-cell'='#F6B1B4','6-cell'='#DB5962',
          '8-cell'='#A71930','HH'='#4477AA','HL'='#EE7733','LH'='#228833','LL'='#AA3377')
# Axis limits are user-defined, and applied with coord_cartesian (never delete observations).
axis_limit <- function(mode,lo,hi) {
  if(!identical(mode,'Custom') || !is.finite(lo) || !is.finite(hi) || lo>=hi) return(NULL)
  coord_cartesian(ylim=c(lo,hi))
}

# Robust source-table identification. This avoids relying on Medium text such as "5X".
table_text <- function(x) str_to_lower(as.character(x))
is_table04 <- function(x) grepl("table\\s*0?4|^0?4$|si[_ ]?table\\s*0?4", table_text(x))
is_table05 <- function(x) grepl("table\\s*0?5|^0?5$|si[_ ]?table\\s*0?5", table_text(x))

# Canonical stage parser used ONLY for plotting.
canonical_stage <- function(origin, stage, sample=NULL) {
  txt <- str_to_lower(paste(origin, stage, sample))
  case_when(
    str_detect(txt, "oocyte") ~ "Oocyte",
    str_detect(txt, "zygote") & !str_detect(txt, "[468][-_ ]?cell|[468]cell") ~ "Zygote",
    str_detect(txt, "(^|[^0-9])2[-_ ]?cell|2cell") ~ "2-cell",
    str_detect(txt, "(^|[^0-9])4[-_ ]?cell|4cell") ~ "4-cell",
    str_detect(txt, "(^|[^0-9])6[-_ ]?cell|6cell") ~ "6-cell",
    str_detect(txt, "(^|[^0-9])8[-_ ]?cell|8cell") ~ "8-cell",
    TRUE ~ NA_character_
  )
}

# Canonical sequential-label condition for 48 h data.
canonical_condition <- function(condition, label1, label2) {
  x <- str_to_upper(as.character(condition))
  fallback <- paste0(str_to_upper(as.character(label1)), str_to_upper(as.character(label2)))
  out <- case_when(
    str_detect(x, "(^|[^A-Z])HH([^A-Z]|$)") ~ "HH",
    str_detect(x, "(^|[^A-Z])HL([^A-Z]|$)") ~ "HL",
    str_detect(x, "(^|[^A-Z])LH([^A-Z]|$)") ~ "LH",
    str_detect(x, "(^|[^A-Z])LL([^A-Z]|$)") ~ "LL",
    fallback %in% c("HH","HL","LH","LL") ~ fallback,
    TRUE ~ NA_character_
  )
  out
}
ui <- navbarPage('pulsed SILAC in Embryos Atlas',id='nav',header=tags$head(tags$style(HTML('
body{background:#f7f8fb}.container-fluid{max-width:1700px}
.cardx{background:white;border:1px solid #e4e8ed;border-radius:10px;padding:18px;margin:12px 0}
.toprow{display:flex;align-items:stretch}
.toprow>[class*="col-"]{display:flex;align-items:stretch}
.toprow .shiny-html-output{width:100%;display:block}
.topcard{height:220px;width:100%;overflow:visible;margin:12px 0;box-sizing:border-box}
.metricgrid{display:grid;grid-template-columns:1fr 1fr;grid-template-rows:1fr 1fr;gap:10px;height:220px;width:100%;padding:12px 0;box-sizing:border-box}
.metricbox{background:white;border:1px solid #e4e8ed;border-radius:10px;padding:16px;display:flex;flex-direction:column;justify-content:center;min-width:0}
.metricvalue{font-size:25px;font-weight:700}.metriclabel{font-size:12px;color:#687487;margin-top:4px}
.smallnote{font-size:12px;color:#667085}.dataTables_wrapper{width:100%}.discovery-table{font-size:11px}.discovery-table table.dataTable tbody td{vertical-align:top;white-space:normal!important;word-break:normal}.discovery-table table.dataTable thead th{font-size:11px;white-space:normal!important}.selection-status{padding:8px 0;font-weight:600;color:#475467}
.searchhit{font-size:12px;color:#667085;margin-top:5px}
.topcard .form-group{margin-bottom:10px}
.topcard h2{margin-top:4px;margin-bottom:8px;font-size:30px}
.topcard h4{margin-top:2px}
.topcard p{margin-bottom:8px}
.metricbox{overflow:hidden}
.discovery-table table.dataTable thead input{width:100%!important;box-sizing:border-box;padding-right:20px!important;min-width:58px}
.network-note{background:#f8fafc;border-left:4px solid #98a2b3;padding:10px 12px;margin:8px 0 14px 0;font-size:12px;color:#475467}
.network-controls .form-group{margin-bottom:10px}
.network-legend{background:#fbfcfe;border:1px solid #dfe3e8;border-radius:8px;padding:10px 12px;margin-top:10px;font-size:12px}
.network-busy-overlay{display:none;position:fixed;z-index:99999;left:50%;top:50%;transform:translate(-50%,-50%);background:rgba(255,255,255,.96);border:1px solid #d0d5dd;border-radius:12px;padding:18px 24px;box-shadow:0 8px 28px rgba(16,24,40,.18);min-width:230px;text-align:center;color:#344054;font-size:14px;font-weight:600}
.shiny-busy .network-busy-overlay{display:block}
.network-spinner{width:30px;height:30px;margin:0 auto 10px auto;border:4px solid #e4e7ec;border-top-color:#667085;border-radius:50%;animation:network-spin .8s linear infinite}
@keyframes network-spin{to{transform:rotate(360deg)}}
.legend-dot{display:inline-block;width:12px;height:12px;border-radius:50%;margin-right:5px;border:1px solid #667085;vertical-align:-1px}
.legend-line{display:inline-block;width:30px;border-top:3px solid #667085;margin:0 6px 3px 0}
.legend-focal{display:inline-block;width:15px;height:15px;border-radius:50%;border:4px solid #344054;margin-right:5px;vertical-align:-3px}
.legend-row{display:flex;align-items:center;gap:7px;margin:4px 0}.legend-swatch{display:inline-block;width:14px;height:14px;border-radius:50%;border:1px solid #667085;flex:0 0 auto}.legend-size{display:inline-block;border-radius:50%;background:#D0D5DD;border:1px solid #667085;vertical-align:middle;flex:0 0 auto}.legend-section{margin-top:8px;padding-top:7px;border-top:1px solid #EAECF0}.legend-title{font-weight:700;margin-bottom:4px}.legend-edge{display:inline-block;width:38px;border-top-style:solid;border-top-color:#667085;margin-right:7px;vertical-align:middle}
')),
 tags$script(HTML("Shiny.addCustomMessageHandler('scrollToId', function(id){var el=document.getElementById(id); if(el){el.scrollIntoView({behavior:'smooth',block:'start'});}});"))),
 tabPanel('Home',
   div(class='cardx',tags$img(src='JVOlogo_large.png',alt='Olsen Group laboratory logo',style='float:right;width:190px;max-width:35%;height:auto;margin:0 0 12px 18px;'),h1('pulsed SILAC in Embryos Atlas'),p('Single-oocyte and embryo protein abundance, incorporation, turnover and developmental protein networks.'),
       actionButton('goexplore','Explore a protein'),actionButton('godiscover','Discover proteins'),actionButton('gonetwork','Explore networks'))
 ),
 tabPanel('Protein Explorer',
 fluidRow(class='toprow',
   column(4,
     div(class='cardx topcard',
       h4('Find a protein'),
       textInput('protein_search','Search gene / protein / UniProt',placeholder='e.g. Ezh2, EZHIP, Q61188...'),
       selectizeInput('protein_hits','Matching proteins',choices=NULL,options=list(dropdownParent='body',maxOptions=100)),
       p(class='smallnote','Search: gene symbol, protein name, synonym or UniProt ID.')
     )
   ),
   column(4,div(style='width:100%;',uiOutput('protein_identity'))),
   column(4,div(style='width:100%;',uiOutput('protein_metrics')))
 ),
 fluidRow(column(6,div(class='cardx',h3('5 day abundance'),plotOutput('abundance',height=320))),
          column(6,div(class='cardx',h3('5 day incorporation / turnover'),plotOutput('incorp',height=320)))),
 div(class='cardx',h3('24 h and 48 h pSILAC — 5X only'),
 fluidRow(column(3,selectInput('metric','Measurement',choices=metric_choices)),
          column(3,selectInput('axis_mode','Y-axis',c('Auto','Custom'))),
          column(3,numericInput('ymin','Y minimum',0)),column(3,numericInput('ymax','Y maximum',30))),
 p(class='smallnote','Intensities are log2-transformed and averaged per replicate. Pulse-Chase 48 h data are described according to the used Label: HH = H -> H, HL = H -> L, LH = L -> H. Labeling was applied for 24 h.'),
 fluidRow(column(6,h4('24 h: Oocyte vs 2-cell embryo'),plotOutput('short24',height=420)),
          column(6,h4('48 h: developmental stages and labeling'),plotOutput('short48',height=420)))),
 div(class='cardx',h3('Identified phosphosites'),
     p(class='smallnote','Motifs and exploratory sequence classes are precomputed from the mouse UniProt FASTA and read from phosphosite_annotations_with_motifs_fixed.csv. Reported-site, functional and kinase annotations are displayed only when supplied by a curated local annotation file; Not annotated does not mean novel. Motif classes are heuristic sequence patterns (S/T-P; R/K at -3/-2; D/E at +1 to +3; at least 3 prolines within ±7), can overlap, and do not establish kinase activity.'),
     selectInput('phospho_report_filter','UniProt annotation',choices=c('All sites','Reported in UniProt','Not annotated in UniProt'),selected='All sites'),
     DTOutput('protein_phosphosites'))),
 tabPanel('Protein Discovery',fluidRow(
  column(3,div(class='cardx',h3('Discovery filters'),
    textInput('discovery_search','Search within results',placeholder='Gene, protein name or UniProt ID'),
    selectInput('discovery_mode','Criterion',c('All proteins'='all','First detected (5 day)'='first',
      'First measurable incorporation (5 day)'='firstinc','Detected on day (5 day)'='detected',
      'Incorporation on day (5 day)'='inc','Persistent detection (5 day)'='persistent',
      'Fast turnover'='fast','Intermediate turnover'='intermediate','Slow turnover'='slow',
      'Turnover not estimable'='notestimable')),
    selectInput('day','Day',choices=days),numericInput('minrep','Minimum replicate measurements',1,min=1,max=max_replicates,step=1),
    sliderInput('threshold','Minimum incorporation',0,1,0,step=.01),
    textOutput('turnover_definition'),
    selectizeInput('localization','Localization',choices=c('All','Mitochondrion','Nucleus','Cytoplasm',
      'Endoplasmic reticulum','Golgi','Lysosome','Peroxisome','Plasma membrane','Ribosome','Extracellular'),selected='All'),
    p(class='smallnote','UniProt and GO annotations are used for the localization filter.'))),
  column(9,div(class='cardx',h3('Discovery results'),textOutput('resultcount'),
    p(class='smallnote','Select one row to explore the proteins network, or multiple rows to compare proteins. Clear selection removes selections retained from other table pages.'),
    fluidRow(
      column(3,actionButton('open_selected','Open selected protein',width='100%')),
      column(3,actionButton('compare_rows','Compare selected proteins',width='100%')),
      column(3,actionButton('network_rows','Explore selected network',width='100%')),
      column(3,actionButton('clear_selection','Clear selection',width='100%'))
    ),
    div(class='selection-status',textOutput('selection_status',inline=TRUE)),
    div(class='discovery-table',DTOutput('discovery'))),
    tags$div(id='compare_anchor'),
    div(class='cardx',h3('Compare proteins'),
      selectizeInput('compare_proteins','Selected proteins',choices=setNames(ann$UniProtID,ann$SearchLabel),
        multiple=TRUE,options=list(placeholder='Search gene, protein or UniProt ID...')),
      selectInput('compare_metric','Comparison measurement',c('log2 abundance'='log2Abundance','Incorporation'='Incorporation')),
      selectInput('compare_axis','Comparison y-axis',c('Auto','Custom')),
      fluidRow(column(6,numericInput('compare_min','Minimum',0)),column(6,numericInput('compare_max','Maximum',30))),
      fluidRow(
        column(6,actionButton('use_table_selection','Use selected table rows',width='100%')),
        column(6,actionButton('clear_compare','Clear comparison',width='100%'))
      ),
      textOutput('compare_status'),
      plotOutput('compareplot',height=440),
      p(class='smallnote','These are descriptive trajectories, not statistical correlations. Different proteins may have different missingness.'))))),
 tabPanel('Network',
   div(class='network-busy-overlay',div(class='network-spinner'),div('Calculating network…'),div(class='smallnote',style='margin-top:4px;font-weight:400;','Large protein groups can take a few seconds.')),
   tabsetPanel(id='network_workflow',
     tabPanel('Protein neighborhood',
       fluidRow(
         column(3,div(class='cardx network-controls',
           h3('Find a protein of interest'),
           selectizeInput('network_protein','Gene / protein / UniProt',
             choices=setNames(summary$UniProtID, paste0(summary$GeneSymbol,' · ',summary$UniProtID)),
             multiple=FALSE,selected=if(nrow(summary)) summary$UniProtID[1] else character(0),
             options=list(placeholder='Search one protein...')),
           actionButton('network_use_explorer','Use Protein Explorer protein',width='100%'),
           selectInput('network_type','Neighborhood',c('Correlated proteins'='correlation','STRING partners'='string','Combined evidence'='combined')),
           conditionalPanel("input.network_type == 'correlation' || input.network_type == 'combined'",
             selectInput('corr_metric','5 day trajectory',c('Abundance (log2)'='log2Abundance','Incorporation'='Incorporation')),
             selectInput('corr_method','Correlation',c('Pearson'='pearson','Spearman'='spearman')),
             sliderInput('corr_threshold','Minimum |r|',min=.5,max=.99,value=.85,step=.01),
             selectInput('corr_sign','Correlations',c('Positive only'='positive','Positive and negative'='both')),
             sliderInput('corr_minpoints','Minimum shared measured days',min=3,max=5,value=4,step=1)),
           conditionalPanel("input.network_type == 'string' || input.network_type == 'combined'",
             sliderInput('string_score','STRING confidence',min=0,max=1,value=.7,step=.05)),
           numericInput('max_edges','Maximum neighbors displayed',value=25,min=1,max=150,step=5),
           selectInput('node_color','Node color',c('Turnover class'='turnover','First detected day'='firstday','Localization'='localization','Maximum incorporation'='maxincorp','Detected dataset'='dataset')),
           selectInput('node_size','Node size',c('Uniform'='uniform','Mean 5 day abundance'='abundance','Maximum incorporation'='incorp')),
           p(class='smallnote','One protein of interest is compared with every eligible protein in the atlas. Correlations use observed 5 day means only. STRING partners are restricted to detected atlas proteins.'))),
         column(9,
           div(class='cardx',h3('Protein-centered developmental network'),
             div(class='network-note','Correlation is exploratory with only five developmental days and does not establish physical binding. STRING reports known/predicted functional associations, not necessarily direct binding. Combined shows partners supported by both.'),
             textOutput('network_status'),
             visNetworkOutput('protein_network',height='610px'),
             uiOutput('network_legend')),
           div(class='cardx',h3('Selected neighbor'),uiOutput('network_node_info')),
           div(class='cardx',h3('Ranked neighbors'),DTOutput('network_edges'))
         )
       )
     ),
     tabPanel('Explore protein groups',
       fluidRow(
         column(3,div(class='cardx network-controls',
           h3('Define a protein group'),
           selectizeInput('group_localization','Localization contains',choices=c('All proteins'='All'),selected='All',multiple=FALSE,
             options=list(create=TRUE,placeholder='e.g. Mitochondrial')),
           selectInput('group_firstday','First detected',c('Any day'='any','Day 1'='1','Day 2'='2','Day 3'='3','Day 4'='4','Day 5'='5','Day 2 or later'='2plus','Day 3 or later'='3plus')),
           selectInput('group_turnover','Turnover',c('Any'='any','Fast'='Fast','Intermediate'='Intermediate','Slow'='Slow','Not estimable'='Not estimable')),
           selectInput('group_metric','5 day trajectory',c('Abundance (log2)'='log2Abundance','Incorporation'='Incorporation')),
           selectInput('group_method','Correlation',c('Pearson'='pearson','Spearman'='spearman')),
           sliderInput('group_r','Minimum |r|',min=.5,max=.99,value=.85,step=.01),
           selectInput('group_sign','Correlations',c('Positive only'='positive','Positive and negative'='both')),
           sliderInput('group_minpoints','Minimum shared measured days',min=3,max=5,value=4,step=1),
           checkboxInput('group_match_turnover','Require similar half-life',FALSE),
           conditionalPanel("input.group_match_turnover == true",numericInput('group_halfdiff','Maximum half-life difference (days)',value=0.5,min=0.05,max=10,step=0.1)),
           numericInput('group_maxproteins','Maximum proteins in network',value=100,min=10,max=300,step=10),
           numericInput('group_maxedges','Maximum edges displayed',value=250,min=10,max=1500,step=25),
           selectInput('group_node_color','Node color',c('Module'='module','Turnover class'='turnover','First detected day'='firstday','Localization'='localization','Maximum incorporation'='maxincorp')),
           selectInput('group_node_size','Node size',c('Uniform'='uniform','Mean 5 day abundance'='abundance','Maximum incorporation'='incorp')),
           p(class='smallnote','This mode discovers coordinated modules within a biologically defined subset. Missing values are not imputed; correlations use only shared measured days.'))),
         column(9,
           div(class='cardx',h3('Data-driven protein-group network'),
             div(class='network-note','Edges represent similar developmental trajectories within the selected protein group. They are descriptive correlations, not evidence of physical interaction.'),
             textOutput('group_status'),
             tableOutput('group_diagnostics'),
             visNetworkOutput('group_network',height='610px'),
             uiOutput('group_legend')),
           div(class='cardx',h3('Selected protein'),uiOutput('group_node_info')),
           fluidRow(column(5,div(class='cardx',h3('Modules'),DTOutput('group_modules'))),
                    column(7,div(class='cardx',h3('Correlated pairs'),DTOutput('group_edges'))))
         )
       )
     )
   )
 ),
 tabPanel('About',div(class='cardx',h2('Interpretation'),
   p('Missing MS values do not establish absent expression. First measurable incorporation is not necessarily the onset of synthesis.'),
   p('Light amino acids may be recycled into newly synthesized proteins; heavy and light are not synonymous with new and old.'),
   p('Half-life model: I(t) = 1 - exp(-kt), with t=0 at the start of heavy labeling.'),
   p('App QC enforces that incorporation is only displayed when a positive SumIntensity exists for the same replicate.'),
   h3('Associated manuscript'),
   p(tags$a('Schmidt et al. (2026), Pulsed-SILAC in single mouse embryos reveals early embryonic protein synthesis dynamics and phosphosite regulation (bioRxiv)',href='https://www.biorxiv.org/content/10.64898/2026.09.24.754052v1',target='_blank',rel='noopener noreferrer')))))
server <- function(input,output,session) {
 observeEvent(input$goexplore,updateNavbarPage(session,'nav',selected='Protein Explorer'))
 observeEvent(input$godiscover,updateNavbarPage(session,'nav',selected='Protein Discovery'))
 observeEvent(input$gonetwork,updateNavbarPage(session,'nav',selected='Network'))
 # ---------------------------------------------------------------------------
 # Active protein state
 # ---------------------------------------------------------------------------
 # The plots use this reactiveVal directly. This guarantees that opening a
 # protein from Discovery changes the actual Explorer state, not only a widget.
 initial_id <- if(nrow(ann)) ann$UniProtID[1] else NA_character_
 active_protein <- reactiveVal(initial_id)
 comparison_ids <- reactiveVal(character(0))

 observeEvent(input$protein_search,{
   q <- str_trim(str_to_lower(input$protein_search))
   if(nchar(q) < 2) {
     updateSelectInput(session,'protein_hits',choices=character(0),selected=character(0))
     return()
   }
   hit <- ann %>%
     filter(grepl(q,SearchText,fixed=TRUE)) %>%
     slice_head(n=100)
   ch <- setNames(hit$UniProtID,ifelse(is.na(hit$GeneSymbol) | !nzchar(trimws(hit$GeneSymbol)),hit$UniProtID,hit$GeneSymbol))
   updateSelectInput(session,'protein_hits',choices=ch,
                     selected=if(length(ch)) unname(ch[1]) else character(0))
 },ignoreInit=TRUE)

 observeEvent(input$protein_hits,{
   req(input$protein_hits)
   active_protein(as.character(input$protein_hits))
 },ignoreInit=TRUE)

 chosen <- reactive({
   id <- active_protein()
   req(id,!is.na(id),nzchar(id))
   id
 })

 sel <- reactive(summary %>% filter(UniProtID==chosen()) %>% dplyr::slice(1))
 fmtday <- function(x) if(length(x) && is.finite(x)) paste('Day',x) else 'NA'

 output$protein_identity <- renderUI({
   a <- ann %>% filter(UniProtID==chosen()) %>% dplyr::slice(1)
   req(nrow(a))
   div(class='cardx topcard',
       h2(coalesce(a$GeneSymbol,a$UniProtID)),
       p(a$ProteinName),
       p(strong('UniProt: '),a$UniProtID),
       if(!is.na(a$GeneSynonyms) && nzchar(a$GeneSynonyms))
         p(class='smallnote',strong('Synonyms: '),a$GeneSynonyms))
 })

 output$protein_metrics <- renderUI({
   s0 <- sel()
   h <- if(nrow(s0) && 'half_life' %in% names(s0) && is.finite(s0$half_life[1]) &&
           (!'fit_status' %in% names(s0) || identical(as.character(s0$fit_status[1]),'fit_ok')))
          paste(round(s0$half_life[1],2),'d') else 'NA'
   det <- if(nrow(s0) && 'DetectedIn' %in% names(s0) && !is.na(s0$DetectedIn[1])) as.character(s0$DetectedIn[1]) else 'NA'
   sites <- if(nrow(s0) && 'PhosphositeCount' %in% names(s0)) as.character(s0$PhosphositeCount[1]) else '0'

   div(class='metricgrid',
     div(class='metricbox',div(class='metricvalue',fmtday(s0$FirstDetectedDay)),
         div(class='metriclabel','First detected')),
     div(class='metricbox',div(class='metricvalue',fmtday(s0$FirstIncorpDay)),
         div(class='metriclabel','First incorporation measured')),
     div(class='metricbox',div(class='metricvalue',h),
         div(class='metriclabel','Protein half-life')),
     div(class='metricbox',div(class='metricvalue',sites),
         div(class='metriclabel','Detected phosphosites'))
   )
 })

 # Optional, locally curated mouse phosphosite annotations.
 # No database evidence is inferred from the motif or from an absent record.
 annotation_path <- file.path(root,'shiny_data','phosphosite_annotations_with_motifs_fixed.csv')
 annotation_columns <- c('UniProtID','Residue','Position','Reported','FunctionalAnnotation','KinaseEvidence','SequenceMotif','MotifClass','MotifQC')
 if(file.exists(annotation_path)) {
   site_annotations <- tryCatch(read.csv(annotation_path,stringsAsFactors=FALSE,
                                         check.names=FALSE),error=function(e) {
     warning('Phosphosite annotation CSV could not be read: ',conditionMessage(e))
     data.frame()
   })
 } else site_annotations <- data.frame()
 if(!all(annotation_columns %in% names(site_annotations))) {
   if(nrow(site_annotations)) warning('Motif annotation CSV missing required columns; annotations disabled.')
   site_annotations <- as.data.frame(setNames(replicate(length(annotation_columns),
                                        character(0),simplify=FALSE),annotation_columns))
 }
 site_annotations <- site_annotations %>%
   mutate(UniProtID=as.character(UniProtID),
          Residue=toupper(trimws(as.character(Residue))),
          Position=suppressWarnings(as.integer(Position))) %>%
   filter(!is.na(Position),Residue %in% c('S','T','Y'),
          !is.na(UniProtID),nzchar(UniProtID)) %>%
   distinct(UniProtID,Residue,Position,.keep_all=TRUE)
 message('[Phosphosite annotations / precomputed motifs] ',nrow(site_annotations),
         ' unique accession/residue/position records loaded from ',
         if(file.exists(annotation_path)) annotation_path else 'no precomputed annotation file (copy the fixed CSV into shiny_data)')

 output$protein_phosphosites <- DT::renderDT({
   z <- phospho_detail %>% filter(UniProtID == chosen())
   validate(need(nrow(z)>0,'No identified phosphosites for this protein in the available phosphoproteomics datasets.'))
   z <- z %>% group_by(UniProtID,PhosphositeID) %>%
     summarise(`Maximum incorporation (%)`=
       if(any(is.finite(Incorporation))) 100*max(Incorporation[is.finite(Incorporation)])
       else NA_real_,.groups='drop')
   parsed <- regmatches(z$PhosphositeID,
     regexpr('(?<![A-Za-z])[STYsty]\\s*[0-9]+',z$PhosphositeID,perl=TRUE))
   z$Residue <- ifelse(nzchar(parsed),toupper(substr(parsed,1,1)),NA_character_)
   z$Position <- suppressWarnings(as.integer(sub('^[STYsty]\\s*','',parsed)))
   z <- left_join(z,site_annotations,by=c('UniProtID','Residue','Position'))
   z$`Motif (+/-7 AA)` <- as.character(z$SequenceMotif)
   z$`Motif class` <- as.character(z$MotifClass)
   z$`Motif class`[is.na(z$`Motif class`) | !nzchar(trimws(z$`Motif class`))] <- 'Not mapped'
   evidence <- function(x) {
     x <- trimws(as.character(x))
     ifelse(is.na(x)|!nzchar(x),'Not annotated',x)
   }
   z$`Previously reported` <- evidence(z$Reported)
   z$`Functional annotation` <- evidence(z$FunctionalAnnotation)
   z$`Kinase evidence` <- evidence(z$KinaseEvidence)
   if(identical(input$phospho_report_filter,'Reported in UniProt'))
     z <- z[!is.na(z$Reported) & nzchar(trimws(z$Reported)),,drop=FALSE]
   if(identical(input$phospho_report_filter,'Not annotated in UniProt'))
     z <- z[is.na(z$Reported) | !nzchar(trimws(z$Reported)),,drop=FALSE]
   z$`Maximum incorporation (%)` <- round(z$`Maximum incorporation (%)`,1)
   z <- z %>% select(`Phosphosite ID`=PhosphositeID,
     `Motif (+/-7 AA)`,`Motif class`,`Maximum incorporation (%)`,
     `Previously reported`,`Functional annotation`,`Kinase evidence`)
   DT::datatable(z,rownames=FALSE,filter='top',extensions='Buttons',
     options=list(pageLength=10,scrollX=TRUE,dom='Bftip',buttons=c('csv')))
 },server=TRUE)

 # Draw lines only between consecutive measured days. Missing days remain gaps.
 draw_consecutive <- function(day,value,pch=19,lwd=2,lty=1) {
   ord <- order(day); day <- day[ord]; value <- value[ord]
   points(day,value,pch=pch)
   if(length(day) >= 2) {
     for(i in seq_len(length(day)-1)) {
       if(is.finite(day[i]) && is.finite(day[i+1]) &&
          is.finite(value[i]) && is.finite(value[i+1]) &&
          (day[i+1]-day[i]) == 1) {
         segments(day[i],value[i],day[i+1],value[i+1],lwd=lwd,lty=lty)
       }
     }
   }
 }

 output$abundance <- renderPlot({
   d <- daily %>% filter(UniProtID==chosen(),is.finite(Day),is.finite(log2Abundance)) %>% arrange(Day)
   validate(need(nrow(d)>0,'No measured 5 day abundance.'))

   yl <- range(c(d$log2Abundance-d$sdAbundance,d$log2Abundance+d$sdAbundance),na.rm=TRUE)
   if(!all(is.finite(yl)) || diff(yl)==0) yl <- range(d$log2Abundance,na.rm=TRUE)+c(-1,1)

   plot(d$Day,d$log2Abundance,type='n',
        xlab='Day',ylab='Mean log2 SumIntensity',ylim=yl,
        xlim=c(1,5),xaxt='n')
   axis(1,at=1:5,labels=1:5)
   draw_consecutive(d$Day,d$log2Abundance,pch=19,lwd=2)
   ok <- is.finite(d$sdAbundance)
   if(any(ok)) arrows(d$Day[ok],d$log2Abundance[ok]-d$sdAbundance[ok],
                      d$Day[ok],d$log2Abundance[ok]+d$sdAbundance[ok],
                      angle=90,code=3,length=.05)
   grid()
 })

 output$incorp <- renderPlot({
   d <- daily %>%
     filter(UniProtID==chosen(),is.finite(Day),is.finite(Incorporation)) %>%
     arrange(Day)
   validate(need(nrow(d)>0,'No measured 5 day incorporation.'))

   # Measured mean incorporation trajectory.
   plot(d$Day,d$Incorporation,type='n',
        xlab='Day',ylab='Mean incorporation',
        xlim=c(1,5),ylim=c(0,1),
        xaxt='n',yaxt='n')
   axis(1,at=1:5,labels=1:5)
   draw_consecutive(d$Day,d$Incorporation,pch=19,lwd=2)
   axis(2,at=seq(0,1,.2),labels=paste0(seq(0,100,20),'%'),las=1)

   ok <- is.finite(d$sdIncorp)
   if(any(ok)) {
     arrows(d$Day[ok],
            pmax(0,d$Incorporation[ok]-d$sdIncorp[ok]),
            d$Day[ok],
            pmin(1,d$Incorporation[ok]+d$sdIncorp[ok]),
            angle=90,code=3,length=.05)
   }

   # Turnover model from the master analysis:
   # I(t) = 1 - exp(-k*t), with I(0)=0 fixed.
   h <- hl %>% filter(UniProtID==chosen()) %>% dplyr::slice(1)
   fit_drawn <- FALSE
   if(nrow(h) && is.finite(h$k[1]) && h$k[1] > 0 &&
      (!'fit_status' %in% names(h) ||
       identical(as.character(h$fit_status[1]),'fit_ok'))) {

     k_fit <- as.numeric(h$k[1])
     # Same fixed-zero turnover model used in the master analysis.
     # Only the experimentally displayed interval (days 1–5) is drawn.
     xx <- seq(1,5,length.out=500)
     yy <- 1 - exp(-k_fit * xx)

     lines(xx,yy,lwd=3,lty=2)
     fit_drawn <- TRUE
   }

   if(fit_drawn) {
     legend('bottomright',
            legend=c('Measured mean ± SD','Turnover fit: 1 - exp(-kt)'),
            lty=c(1,2),lwd=c(2,3),pch=c(19,NA),
            bty='n',cex=.8)
   }

   grid()
 })

 shortdata <- reactive({
   req(input$metric, chosen())
   d <- short %>% filter(UniProtID==chosen())
   if(!nrow(d)) return(d)

   # Table05 = 24 h 5X oocyte/zygote experiment.
   # Table04 = 48 h 5X sequential-label experiment.
   # Table03 (1X) is deliberately excluded from Protein Explorer.
   keep <- is_table04(d$Table) | is_table05(d$Table)
   d <- d[!is.na(keep) & keep,,drop=FALSE]
   if(!nrow(d)) return(d)

   d$Value <- suppressWarnings(as.numeric(d[[input$metric]]))
   if(input$metric %in% c('SumIntensity','H','L')) d$Value <- log2_ok(d$Value)

   d$PlotStage <- canonical_stage(
     d$Origin,
     d$Stage,
     if('Sample' %in% names(d)) d$Sample else NA_character_
   )

   d$PlotCondition <- canonical_condition(
     d$Condition,
     if('Label_0_24h' %in% names(d)) d$Label_0_24h else NA_character_,
     if('Label_24_48h' %in% names(d)) d$Label_24_48h else NA_character_
   )
   d
 })

 y_axis_label <- reactive({
   switch(input$metric,
     SumIntensity='log2 SumIntensity',
     H='log2 Heavy intensity',
     L='log2 Light intensity',
     Incorporation='Incorporation',
     Ratio='H/L ratio',
     input$metric
   )
 })

 output$short24 <- renderPlot({
   d <- shortdata() %>% filter(is_table05(Table),PlotStage %in% c('Oocyte','Zygote'))
   d$PlotStage <- ifelse(as.character(d$PlotStage)=='Zygote','2-cell',as.character(d$PlotStage))
   lev <- c('Oocyte','2-cell')
   d$PlotStage <- factor(d$PlotStage,levels=lev)
   validate(need(any(is.finite(d$Value)),'Protein was not identified in the 24 h dataset.'))

   if(input$axis_mode=='Custom' && is.finite(input$ymin) && is.finite(input$ymax) && input$ymin<input$ymax)
     yl <- c(input$ymin,input$ymax)
   else if(input$metric=='Incorporation') yl <- c(0,1)
   else {
     yl <- range(d$Value[is.finite(d$Value)],na.rm=TRUE)
     pad <- max(diff(yl)*.08,.25); yl <- yl+c(-pad,pad)
   }

   boxplot(Value ~ PlotStage,data=d,drop=FALSE,
           xlab='',ylab=y_axis_label(),ylim=yl,outline=FALSE,
           col=unname(cols[lev]),
           yaxt=if(input$metric=='Incorporation') 'n' else 's')
   if(input$metric=='Incorporation')
     axis(2,at=seq(0,1,.2),labels=paste0(seq(0,100,20),'%'),las=1)
   stripchart(Value ~ PlotStage,data=d,vertical=TRUE,method='jitter',
              pch=16,cex=.7,add=TRUE)
   grid()
 })

 output$short48 <- renderPlot({
   d <- shortdata() %>%
     filter(is_table04(Table),
            PlotStage %in% c('4-cell','6-cell','8-cell'),
            PlotCondition %in% c('HL','LH','HH'))

   validate(need(any(is.finite(d$Value)),'Protein was not identified in the 48 h dataset.'))

   # Requested order:
   # HL: 4-cell -> 6-cell -> 8-cell
   # LH: 4-cell -> 6-cell -> 8-cell
   # HH: 4-cell -> 6-cell -> 8-cell
   group_levels <- c(
     'HL | 4-cell','HL | 6-cell','HL | 8-cell',
     'LH | 4-cell','LH | 6-cell','LH | 8-cell',
     'HH | 4-cell','HH | 6-cell','HH | 8-cell'
   )
   d <- d %>% mutate(Group=factor(paste(PlotCondition,PlotStage,sep=' | '),
                                  levels=group_levels))

   if(input$axis_mode=='Custom' && is.finite(input$ymin) && is.finite(input$ymax) && input$ymin<input$ymax)
     yl <- c(input$ymin,input$ymax)
   else if(input$metric=='Incorporation') yl <- c(0,1)
   else {
     yl <- range(d$Value[is.finite(d$Value)],na.rm=TRUE)
     pad <- max(diff(yl)*.08,.25); yl <- yl+c(-pad,pad)
   }

   box_cols <- rep(c(cols['4-cell'],cols['6-cell'],cols['8-cell']),3)
   boxplot(Value ~ Group,data=d,drop=FALSE,xaxt='n',
           xlab='',ylab=y_axis_label(),ylim=yl,outline=FALSE,
           col=unname(box_cols),
           yaxt=if(input$metric=='Incorporation') 'n' else 's')
   if(input$metric=='Incorporation')
     axis(2,at=seq(0,1,.2),labels=paste0(seq(0,100,20),'%'),las=1)

   stripchart(Value ~ Group,data=d,vertical=TRUE,method='jitter',
              pch=16,cex=.55,add=TRUE)

   axis(1,at=1:9,
        labels=rep(c('4','6','8'),3),
        cex.axis=.82)
   mtext(c('HL','LH','HH'),side=1,at=c(2,5,8),line=2.1,font=2)
   mtext('Cell stage embryo',side=1,line=3.4)
   legend('topright',legend=c('4-cell','6-cell','8-cell'),
          fill=unname(cols[c('4-cell','6-cell','8-cell')]),
          bty='n',cex=.8)
   abline(v=c(3.5,6.5),lty=3)
   grid()
 })

 output$turnover_definition <- renderText({
   if(is.finite(fast_cut) && is.finite(slow_cut)) {
     paste0('Turnover classes use valid fitted half-lives: Fast = lowest quartile (≤ ',
            round(fast_cut,2),' d); Intermediate = middle 50%; Slow = highest quartile (≥ ',
            round(slow_cut,2),' d).')
   } else {
     'Turnover classes require valid fitted half-lives.'
   }
 })

 results <- reactive({
   d <- summary; mode <- input$discovery_mode; day <- as.numeric(input$day)
   if(mode=='first') d <- d %>% filter(FirstDetectedDay==day)
   if(mode=='firstinc') d <- d %>% filter(FirstIncorpDay==day,MaxIncorp>=input$threshold)
   if(mode %in% c('detected','inc')) {
     atday <- daily %>% filter(Day==day) %>% select(UniProtID,nAbundance,nIncorp,Incorporation)
     d <- d %>% left_join(atday,by='UniProtID')
     if(mode=='detected') d <- d %>% filter(nAbundance>=input$minrep)
     else d <- d %>% filter(nIncorp>=input$minrep,Incorporation>=input$threshold)
   }
   if(mode=='persistent') d <- d %>% filter(DaysDetected==length(days))
   if(mode=='fast') d <- d %>% filter(TurnoverClass=='Fast')
   if(mode=='intermediate') d <- d %>% filter(TurnoverClass=='Intermediate')
   if(mode=='slow') d <- d %>% filter(TurnoverClass=='Slow')
   if(mode=='notestimable') d <- d %>% filter(TurnoverClass=='Not estimable')
   if(!is.null(input$localization) && input$localization!='All') {
     term <- input$localization
     # Plasma membrane and cell membrane are treated as equivalent search terms.
     pattern <- if(term=='Plasma membrane') 'plasma membrane|cell membrane' else term
     matchtext <- paste(d$UniProtLocation,d$GOCellularComponent)
     d <- d[grepl(pattern,matchtext,ignore.case=TRUE),,drop=FALSE]
   }

   # Free-text search inside the currently filtered result set.
   q <- str_trim(str_to_lower(coalesce(input$discovery_search,'')))
   if(nzchar(q)) {
     txt <- str_to_lower(paste(d$GeneSymbol,d$ProteinName,d$UniProtID,d$GeneSynonyms,sep=' | '))
     d <- d[grepl(q,txt,fixed=TRUE),,drop=FALSE]
   }

   d %>% arrange(is.na(FirstDetectedDay),FirstDetectedDay,GeneSymbol)
 })
 output$resultcount <- renderText(paste(nrow(results()),'proteins match your filters'))
 output$discovery <- renderDT({
   d <- results() %>%
     transmute(
       Gene=GeneSymbol,
       `Protein name`=ProteinName,
       UniProt=UniProtID,
       `First day`=FirstDetectedDay,
       `First incorp.`=FirstIncorpDay,
       `Days seen`=DaysDetected,
       `Max incorp.`=MaxIncorp,
       `Half-life (d)`=half_life,
       Turnover=TurnoverClass,
       `Detected in`=DetectedIn,
       Localization=coalesce(UniProtLocation,GOCellularComponent)
     )

   # Keep the full protein name in a tooltip while showing a compact label.
   full_name <- ifelse(is.na(d$`Protein name`),'',d$`Protein name`)
   short_name <- ifelse(nchar(full_name)>58,paste0(substr(full_name,1,55),'...'),full_name)
   d$`Protein name` <- paste0('<span title="',
     htmltools::htmlEscape(full_name,attribute=TRUE),'">',
     htmltools::htmlEscape(short_name),'</span>')

   datatable(d,rownames=FALSE,selection='multiple',filter='top',escape=FALSE,
      extensions='Buttons',class='compact stripe hover',
      options=list(
       scrollX=FALSE,pageLength=15,autoWidth=FALSE,
       columnDefs=list(
         list(width='8%',targets=0),
         list(width='22%',targets=1),
         list(width='7%',targets=2),
         list(width='6%',targets=c(3,4,5)),
         list(width='7%',targets=c(6,7)),
         list(width='8%',targets=8),
         list(width='9%',targets=9),
         list(width='15%',targets=10)
       ),
       dom='Bfrtip',buttons=c('copy','csv')
      )) %>%
      formatPercentage('Max incorp.',1) %>%
      formatRound('Half-life (d)',2)
 },server=FALSE)

 # Track the most recently added table row, so Open selected protein never
 # silently opens an older selection from another DT page.
 previous_selection <- reactiveVal(integer(0))
 last_selected_id <- reactiveVal(NA_character_)
 observeEvent(input$discovery_rows_selected,{
   now <- input$discovery_rows_selected
   if(is.null(now)) now <- integer(0)
   added <- setdiff(now,previous_selection())
   if(length(added)) {
     d <- results()
     idx <- tail(added,1)
     if(idx >= 1 && idx <= nrow(d)) last_selected_id(as.character(d$UniProtID[idx]))
   }
   previous_selection(now)
 })

 output$selection_status <- renderText({
   n <- length(input$discovery_rows_selected)
   if(!n) 'No proteins selected.' else paste(n,'protein(s) selected.')
 })

 observeEvent(input$clear_selection,{
   dataTableProxy('discovery') %>% selectRows(NULL)
   previous_selection(integer(0))
   last_selected_id(NA_character_)
   comparison_ids(character(0))
   updateSelectizeInput(session,'compare_proteins',selected=character(0),server=FALSE)
 })

 observeEvent(input$open_selected,{
   id <- last_selected_id()
   if(is.na(id) || !nzchar(id)) {
     idx <- input$discovery_rows_selected
     req(length(idx)>=1)
     d <- results()
     id <- as.character(d$UniProtID[tail(idx,1)])
   }
   req(!is.na(id),nzchar(id))
   active_protein(id)
   a <- ann %>% filter(UniProtID==id) %>% dplyr::slice(1)
   if(nrow(a)) {
     updateTextInput(session,'protein_search',value=coalesce(a$GeneSymbol,id))
     updateSelectInput(session,'protein_hits',choices=setNames(a$UniProtID,ifelse(is.na(a$GeneSymbol) | !nzchar(trimws(a$GeneSymbol)),a$UniProtID,a$GeneSymbol)),selected=id)
   }
   updateNavbarPage(session,'nav',selected='Protein Explorer')
 })



 # ---------------------------------------------------------------------------
 # Protein comparison
 # ---------------------------------------------------------------------------
 # Selected rows in the Discovery table are the primary comparison state.
 # No intermediate widget update is required for the plot to react.
 table_selected_ids <- reactive({
   idx <- input$discovery_rows_selected
   d <- results()

   if(is.null(idx) || length(idx)==0 || nrow(d)==0)
     return(character(0))

   idx <- idx[idx >= 1 & idx <= nrow(d)]
   ids <- unique(as.character(d$UniProtID[idx]))
   ids[!is.na(ids) & nzchar(ids)]
 })

 # Comparison is explicit: table selection alone does not overwrite an existing comparison.
 # Clicking Compare selected proteins captures the current rows and scrolls to the plot.
 observeEvent(input$compare_rows,{
   ids <- table_selected_ids()
   req(length(ids) >= 1)
   comparison_ids(ids)
   lab <- ann %>% filter(UniProtID %in% ids)
   updateSelectizeInput(session,'compare_proteins',
     choices=setNames(ann$UniProtID,ann$SearchLabel),selected=ids,server=FALSE)
   session$sendCustomMessage('scrollToId','compare_anchor')
 })

 observeEvent(input$use_table_selection,{
   ids <- table_selected_ids()
   req(length(ids) >= 1)
   comparison_ids(ids)
   session$sendCustomMessage('scrollToId','compare_anchor')
 })

 # Manual multi-protein search remains available as an independent route.
 observeEvent(input$compare_proteins,{
   ids <- unique(as.character(input$compare_proteins))
   ids <- ids[!is.na(ids) & nzchar(ids)]
   if(length(ids)) comparison_ids(ids)
 },ignoreInit=TRUE)

 observeEvent(input$clear_compare,{
   comparison_ids(character(0))
   updateSelectizeInput(session,'compare_proteins',
                        selected=character(0),server=FALSE)
   proxy <- dataTableProxy('discovery')
   selectRows(proxy,NULL)
 })

 output$compare_status <- renderText({
   n <- length(comparison_ids())
   if(!n) {
     'No proteins selected for comparison.'
   } else {
     ids <- comparison_ids()
     labs <- ann %>%
       filter(UniProtID %in% ids) %>%
       mutate(lbl=coalesce(GeneSymbol,UniProtID)) %>%
       pull(lbl)
     paste0(n,' protein(s): ',paste(labs,collapse=', '))
   }
 })

 output$compareplot <- renderPlot({
   ids <- comparison_ids()
   validate(need(length(ids)>0,
     'Select multiple proteins in the table and click Compare selected proteins, or search them above.'))

   d <- daily %>%
     filter(UniProtID %in% ids) %>%
     left_join(ann %>% select(UniProtID,GeneSymbol),by='UniProtID') %>%
     mutate(Label=coalesce(GeneSymbol,UniProtID))

   metric <- input$compare_metric
   if(metric=='log2Abundance') {
     value_col <- 'log2Abundance'; sd_col <- 'sdAbundance'
     ylab <- 'Mean log2 SumIntensity'
   } else {
     value_col <- 'Incorporation'; sd_col <- 'sdIncorp'
     ylab <- 'Mean incorporation'
   }
   d <- d[is.finite(d[[value_col]]) & is.finite(d$Day),,drop=FALSE]
   validate(need(nrow(d)>0,'No measurements for the selected proteins.'))

   if(input$compare_axis=='Custom' && is.finite(input$compare_min) &&
      is.finite(input$compare_max) && input$compare_min<input$compare_max) {
     yl <- c(input$compare_min,input$compare_max)
   } else if(metric=='Incorporation') {
     yl <- c(0,1)
   } else {
     lo <- d[[value_col]] - d[[sd_col]]
     hi <- d[[value_col]] + d[[sd_col]]
     yl <- range(c(lo,hi,d[[value_col]]),na.rm=TRUE)
     if(!all(is.finite(yl)) || diff(yl)==0) yl <- range(d[[value_col]],na.rm=TRUE)+c(-1,1)
   }

   plot(NA,xlim=c(1,5),ylim=yl,xlab='Day',ylab=ylab,
        xaxt='n',yaxt=if(metric=='Incorporation') 'n' else 's')
   axis(1,at=1:5,labels=1:5)
   if(metric=='Incorporation')
     axis(2,at=seq(0,1,.2),labels=paste0(seq(0,100,20),'%'),las=1)

   labs <- unique(d$Label)
   comparison_cols <- setNames(grDevices::hcl.colors(max(3,length(labs)), palette='Dark 3')[seq_along(labs)], labs)

   for(i in seq_along(labs)) {
     z <- d[d$Label==labs[i],,drop=FALSE]
     z <- z[order(z$Day),,drop=FALSE]
     # All proteins use the same solid line; protein identity is encoded by color.
     # draw_consecutive handles missing days without bridging gaps.
     zz <- z[is.finite(z$Day) & is.finite(z[[value_col]]),,drop=FALSE]
     if(nrow(zz)) {
       points(zz$Day, zz[[value_col]], pch=16, col=comparison_cols[labs[i]])
       if(nrow(zz) >= 2) {
         for(j in seq_len(nrow(zz)-1)) {
           if(zz$Day[j+1] - zz$Day[j] == 1)
             segments(zz$Day[j],zz[[value_col]][j],
                      zz$Day[j+1],zz[[value_col]][j+1],
                      lwd=2,lty=1,col=comparison_cols[labs[i]])
         }
       }
     }
     ok <- is.finite(z[[sd_col]]) & is.finite(z[[value_col]]) & is.finite(z$Day)
     if(any(ok)) arrows(z$Day[ok],z[[value_col]][ok]-z[[sd_col]][ok],
                        z$Day[ok],z[[value_col]][ok]+z[[sd_col]][ok],
                        angle=90,code=3,length=.035,col=comparison_cols[labs[i]])
   }
   legend('topright',legend=labs,pch=16,lty=1,lwd=2,
          col=unname(comparison_cols[labs]),bty='n',cex=.8)
   grid()
 })


 # ---------------------------------------------------------------------------
 # Focal-protein network: search all 5 day trajectories; STRING partners in atlas
 # ---------------------------------------------------------------------------
 observeEvent(input$network_rows,{
   ids <- table_selected_ids()
   req(length(ids)>=1)
   updateSelectizeInput(session,'network_protein',selected=tail(ids,1))
   updateNavbarPage(session,'nav',selected='Network')
 })
 observeEvent(input$network_use_explorer,{
   req(active_protein())
   updateSelectizeInput(session,'network_protein',selected=active_protein())
 })

 # Compact, precomputed day-by-protein matrices avoid repeated atlas-wide joins.
 make_network_matrix <- function(metric) {
   d <- daily %>% filter(Day %in% 1:5) %>%
     select(UniProtID,Day,all_of(metric)) %>%
     pivot_wider(names_from=Day,values_from=all_of(metric),names_prefix='D')
   for(nm in paste0('D',1:5)) if(!nm %in% names(d)) d[[nm]] <- NA_real_
   m <- as.matrix(d[,paste0('D',1:5),drop=FALSE])
   storage.mode(m) <- 'double'; rownames(m) <- d$UniProtID
   m
 }
 abundance_mat <- make_network_matrix('log2Abundance')
 incorporation_mat <- make_network_matrix('Incorporation')

 network_focus <- reactive({
   id <- as.character(input$network_protein)
   validate(need(length(id)==1 && nzchar(id) && id %in% summary$UniProtID,
                 'Select a protein from the atlas.'))
   id
 })

 correlation_neighbors <- reactive({
   id <- network_focus()
   mat <- if(input$corr_metric=='Incorporation') incorporation_mat else abundance_mat
   empty <- tibble(from=character(),to=character(),correlation=double(),shared_days=integer())
   if(!id %in% rownames(mat)) return(empty)
   x <- as.numeric(mat[id,]); okx <- is.finite(x)
   if(sum(okx)<input$corr_minpoints) return(empty)
   n <- rowSums(is.finite(mat) & matrix(okx,nrow(mat),ncol(mat),byrow=TRUE))
   candidates <- which(n>=input$corr_minpoints & rownames(mat)!=id)
   if(!length(candidates)) return(empty)
   rr <- vapply(candidates,function(i) {
     ok <- is.finite(mat[i,]) & okx
     if(sum(ok)<input$corr_minpoints) return(NA_real_)
     suppressWarnings(cor(x[ok],mat[i,ok],method=input$corr_method))
   },numeric(1))
   out <- tibble(from=id,to=rownames(mat)[candidates],correlation=rr,shared_days=n[candidates]) %>%
     filter(is.finite(correlation),abs(correlation)>=input$corr_threshold)
   if(input$corr_sign=='positive') out <- out %>% filter(correlation>0)
   out %>% arrange(desc(abs(correlation)))
 })

 # Fetch once per protein of interest, then filter by confidence locally. An explicit
 # partner API call is essential: /network with multiple IDs only links inputs.
 string_cache <- reactiveVal(list())
 fetch_string_partners <- function(id) {
   empty <- tibble(from=character(),to=character(),string_score=double())
   endpoint <- 'https://version-12-5.string-db.org/api/tsv/interaction_partners'
   focal_gene <- ann$GeneSymbol[match(id, ann$UniProtID)]
   # STRING does not consistently resolve UniProt accessions. Try both
   # accession and the annotated mouse gene symbol, and use the first
   # successful nonempty response. Query all partners before atlas filtering.
   candidates <- unique(c(id, focal_gene))
   candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
   query_one <- function(query) {
     q <- paste0('?identifiers=', URLencode(query,reserved=TRUE),
       '&species=10090&required_score=0&limit=5000&caller_identity=pSILAC_atlas')
     tryCatch({
       con <- url(paste0(endpoint,q),open='rt'); on.exit(close(con),add=TRUE)
       read.delim(con,stringsAsFactors=FALSE,check.names=FALSE)
     },error=function(e) NULL)
   }
   responses <- lapply(candidates, query_one)
   responses <- Filter(function(z) !is.null(z) && nrow(z)>0 &&
     all(c('preferredName_A','preferredName_B','score') %in% names(z)), responses)
   if(!length(responses)) {
     warning('STRING returned no interaction partners for ',id,
       ' (mouse accession/gene symbol). Check connectivity and STRING mapping.')
     return(empty)
   }
   z <- bind_rows(responses)
   # STRING aliases are not guaranteed to match UniProt preferred gene names.
   # Map preferred names and synonyms, but exclude ambiguous atlas mappings.
   symbols <- ann %>% transmute(UniProtID,
     aliases=paste(coalesce(GeneSymbol,''),coalesce(GeneSynonyms,''),sep=';'))
   amap <- symbols %>% mutate(alias=strsplit(aliases,'[;,[:space:]]+')) %>%
     tidyr::unnest(alias) %>% mutate(alias=toupper(trimws(alias))) %>%
     filter(nzchar(alias)) %>% select(UniProtID,alias) %>% distinct() %>%
     add_count(alias) %>% filter(n==1) %>% select(-n)
   mapping <- setNames(amap$UniProtID,amap$alias)
   from_a <- toupper(as.character(z$preferredName_A))
   from_b <- toupper(as.character(z$preferredName_B))
   a_id <- unname(mapping[from_a]); b_id <- unname(mapping[from_b])
   focal_aliases <- toupper(c(id,focal_gene))
   a_id[from_a %in% focal_aliases] <- id
   b_id[from_b %in% focal_aliases] <- id
   partner <- rep(NA_character_,length(a_id))
   is_a <- !is.na(a_id) & a_id==id
   is_b <- !is.na(b_id) & b_id==id
   partner[is_a] <- b_id[is_a]
   partner[is_b] <- a_id[is_b]
   out <- tibble(from=id,to=partner,string_score=as.numeric(z$score)) %>%
     filter(!is.na(to),to!=id,to %in% summary$UniProtID,
       is.finite(string_score)) %>%
     group_by(from,to) %>% summarise(string_score=max(string_score),.groups='drop') %>%
     arrange(desc(string_score))
   if(!nrow(out)) message('STRING: partners found but none mapped to atlas for ',id,
     '. Check species, annotations and confidence threshold.')
   out
 }
 string_neighbors <- reactive({
   id <- network_focus()
   cache <- string_cache()
   if(is.null(cache[[id]])) {
     cache[[id]] <- fetch_string_partners(id)
     string_cache(cache)
   }
   cache[[id]] %>% filter(string_score>=input$string_score)
 })

 active_edges <- reactive({
   mode <- input$network_type
   if(mode=='correlation') {
     e <- correlation_neighbors() %>% mutate(string_score=NA_real_,source='pSILAC correlation')
     if(nrow(e)) e <- e %>% arrange(desc(abs(correlation)))
   } else if(mode=='string') {
     e <- string_neighbors() %>% mutate(correlation=NA_real_,shared_days=NA_integer_,source='STRING')
     if(nrow(e)) e <- e %>% arrange(desc(string_score))
   } else {
     c0 <- correlation_neighbors(); s0 <- string_neighbors()
     e <- inner_join(c0,s0,by=c('from','to')) %>%
       mutate(source='Correlation + STRING')
     if(nrow(e)) e <- e %>% arrange(desc(string_score),desc(abs(correlation)))
   }
   if(nrow(e)>input$max_edges) e <- head(e,input$max_edges)
   e
 })
 # Explicit, stable node colours so the graph and legend always use the same mapping.
 # Broad, interpretable display categories; membership filtering remains separate.
 location_class <- function(x) {
   x <- as.character(x); x[is.na(x)] <- ''
   vapply(x,function(z) {
     if(!nzchar(trimws(z))) return('Unknown')
     flags <- c(Mitochondrial=grepl('mitochondri|mitochon|\\bMICOS\\b',z,ignore.case=TRUE),
       Nuclear=grepl('nucle|chromatin|chromosom',z,ignore.case=TRUE),
       Cytoplasmic=grepl('cytoplas|cytosol',z,ignore.case=TRUE),
       Membrane=grepl('membrane|plasma membrane',z,ignore.case=TRUE),
       Extracellular=grepl('extracellular|secreted',z,ignore.case=TRUE),
       Golgi=grepl('golgi',z,ignore.case=TRUE),
       `Endoplasmic reticulum`=grepl('endoplasmic reticulum',z,ignore.case=TRUE))
     if(any(flags)) return(names(flags)[which(flags)[1]])
     'Other'
   },character(1))
 }
 network_palette <- function(groups, kind='generic') {
   g <- as.character(groups); lev <- unique(g[!is.na(g) & nzchar(g)])
   if(kind=='localization') {
     fixed <- c(Mitochondrial='#C65B66',Nuclear='#7B61A8',Cytoplasmic='#2D9C9B',Membrane='#D59C37',Extracellular='#549B52',Golgi='#4D87C7',`Endoplasmic reticulum`='#A46B42',Other='#89939F',Unknown='#C5CBD3')
     out <- fixed[g]; out[is.na(out)] <- '#89939F'; return(unname(out))
   }
   if(kind=='turnover') {
     fixed <- c('Fast'='#2E86AB','Intermediate'='#F2C14E','Slow'='#D1495B','Not estimable'='#98A2B3','NA'='#D0D5DD')
     out <- fixed[g]; out[is.na(out)] <- '#D0D5DD'; return(unname(out))
   }
   if(kind=='firstday') {
     fixed <- c('Day 1'='#3B82F6','Day 2'='#22C55E','Day 3'='#F2C14E','Day 4'='#F97316','Day 5'='#A855F7','NA'='#D0D5DD')
     out <- fixed[g]; out[is.na(out)] <- '#D0D5DD'; return(unname(out))
   }
   if(kind=='maxincorp') {
     vals <- suppressWarnings(as.numeric(sub('%','',g)))
     pal <- grDevices::colorRampPalette(c('#E8F1FA','#2E86AB'))(11)
     idx <- pmax(1,pmin(11,round(vals/10)+1)); out <- pal[idx]; out[!is.finite(vals)] <- '#D0D5DD'; return(out)
   }
   if(!length(lev)) return(rep('#98A2B3',length(g)))
   cols <- grDevices::hcl.colors(max(length(lev),3),'Dark 3')[seq_along(lev)]
   mp <- setNames(cols,lev); out <- unname(mp[g]); out[is.na(out)] <- '#D0D5DD'; out
 }
 legend_rows <- function(labels, colors) {
   Map(function(lab,col) tags$div(class='legend-row',tags$span(class='legend-swatch',style=paste0('background:',col,';')),tags$span(lab)), labels, colors)
 }
 size_legend_rows <- function(values, sizes, formatter=function(x) format(round(x,2),trim=TRUE)) {
   Map(function(v,sz) tags$div(class='legend-row',tags$span(class='legend-size',style=paste0('width:',sz,'px;height:',sz,'px;')),tags$span(formatter(v))),values,sizes)
 }
 network_nodes <- reactive({
   focus <- network_focus(); e <- active_edges()
   ids <- unique(c(focus,e$to))
   m <- summary %>% filter(UniProtID %in% ids) %>%
     mutate(Label=coalesce(GeneSymbol,UniProtID),
       Localization=coalesce(UniProtLocation,GOCellularComponent,'Unknown'))
   mean_ab <- daily %>% filter(UniProtID %in% ids) %>% group_by(UniProtID) %>%
     summarise(MeanAbundance=mean_ok(log2Abundance),.groups='drop')
   m <- left_join(m,mean_ab,by='UniProtID')
   groups <- switch(input$node_color,
     turnover=m$TurnoverClass,
     firstday=ifelse(is.finite(m$FirstDetectedDay),paste('Day',m$FirstDetectedDay),'NA'),
     localization=location_class(m$Localization),
     maxincorp=ifelse(is.finite(m$MaxIncorp),paste0(round(m$MaxIncorp*10)*10,'%'),'NA'),
     dataset=m$DetectedIn)
   sizes <- rep(22,nrow(m))
   if(input$node_size=='abundance' && any(is.finite(m$MeanAbundance))) {
     x <- m$MeanAbundance; x[!is.finite(x)] <- min(x[is.finite(x)])
     sizes <- if(length(unique(x))>1) rescale(x,to=c(15,38)) else sizes
   }
   if(input$node_size=='incorp' && any(is.finite(m$MaxIncorp))) {
     x <- m$MaxIncorp; x[!is.finite(x)] <- 0
     sizes <- if(length(unique(x))>1) rescale(x,to=c(15,38)) else sizes
   }
   sizes[m$UniProtID==focus] <- max(38,sizes[m$UniProtID==focus])
   node_cols <- network_palette(groups,input$node_color)
   m %>% transmute(id=UniProtID,label=Label,group=as.character(groups),color=node_cols,value=sizes,
     borderWidth=ifelse(UniProtID==focus,5,1),
     title=paste0('<b>',Label,'</b><br>',coalesce(ProteinName,''),
       '<br>UniProt: ',UniProtID,'<br>Detected: ',DetectedIn,
       '<br>Turnover: ',TurnoverClass))
 })
 output$network_status <- renderText({
   id <- network_focus(); e <- active_edges()
   gene <- summary$GeneSymbol[match(id,summary$UniProtID)]
   paste0(coalesce(gene,id),' · ',nrow(e),' neighbors shown · ',input$network_type,
     if(input$network_type!='string') paste0(' · ',nrow(correlation_neighbors()),' qualifying correlations') else '')
 })
 output$protein_network <- renderVisNetwork({
   focus <- network_focus(); e <- active_edges(); n <- network_nodes()
   validate(need(nrow(n)>0,'No protein available.'))
   ve <- if(nrow(e)) e %>% transmute(from,to,
     value=ifelse(is.finite(string_score),pmax(string_score,.01)*5,pmax(abs(correlation),.01)*5),
     title=case_when(source=='pSILAC correlation' ~ paste0('r = ',round(correlation,3),' · shared days = ',shared_days),
       source=='STRING' ~ paste0('STRING score = ',round(string_score,3)),
       TRUE ~ paste0('r = ',round(correlation,3),' · STRING = ',round(string_score,3))),
     dashes=source=='STRING') else data.frame(from=character(),to=character())
   visNetwork(n,ve,width='100%',height='640px') %>%
     visNodes(shape='dot',font=list(size=18)) %>% visEdges(smooth=FALSE) %>%
     visOptions(highlightNearest=list(enabled=TRUE,degree=1,hover=TRUE),nodesIdSelection=TRUE) %>%
     visPhysics(stabilization=TRUE,solver='forceAtlas2Based')
 })
 output$network_node_info <- renderUI({
   id <- input$protein_network_selected
   if(is.null(id) || !length(id)) return(p(class='smallnote','Click a neighbor to inspect it.'))
   id <- as.character(id[1]); a <- summary %>% filter(UniProtID==id) %>% dplyr::slice(1)
   if(!nrow(a)) return(NULL)
   div(h4(coalesce(a$GeneSymbol,a$UniProtID)),p(a$ProteinName),
       p(strong('UniProt: '),a$UniProtID,' · ',strong('Detected in: '),a$DetectedIn),
       p(strong('First detected: '),fmtday(a$FirstDetectedDay),' · ',strong('Turnover: '),a$TurnoverClass),
       actionButton('network_refocus','Make protein of interest'),
       actionButton('open_network_node','Open in Protein Explorer'))
 })
 observeEvent(input$network_refocus,{
   id <- input$protein_network_selected; req(length(id)>=1)
   updateSelectizeInput(session,'network_protein',selected=as.character(id[1]))
 })
 observeEvent(input$open_network_node,{
   id <- input$protein_network_selected; req(length(id)>=1); id <- as.character(id[1])
   active_protein(id)
   a <- ann %>% filter(UniProtID==id) %>% dplyr::slice(1)
   if(nrow(a)) {
     updateTextInput(session,'protein_search',value=coalesce(a$GeneSymbol,id))
     updateSelectInput(session,'protein_hits',choices=setNames(a$UniProtID,ifelse(is.na(a$GeneSymbol) | !nzchar(trimws(a$GeneSymbol)),a$UniProtID,a$GeneSymbol)),selected=id)
   }
   updateNavbarPage(session,'nav',selected='Protein Explorer')
 })
 output$network_edges <- renderDT({
   e <- active_edges()
   validate(need(nrow(e)>0,'No neighbors pass the current settings (or STRING is unavailable).'))
   labels <- ann %>% select(UniProtID,GeneSymbol) %>% distinct(UniProtID,.keep_all=TRUE)
   e %>% left_join(labels,by=c('to'='UniProtID')) %>%
     mutate(Protein=coalesce(GeneSymbol,to)) %>%
     select(Protein,UniProt=to,any_of(c('correlation','shared_days','string_score','source'))) %>%
     datatable(rownames=FALSE,options=list(pageLength=10,scrollX=FALSE))
 })


 # Dynamic legend for the focal-protein network.
 output$network_legend <- renderUI({
   n <- network_nodes(); req(nrow(n)>0)
   # actual colour categories currently visible
   cats <- n %>% distinct(group,color) %>% filter(!is.na(group),nzchar(group)) %>% arrange(group)
   color_label <- switch(input$node_color,turnover='Turnover class',firstday='First detected day',localization='Localization',maxincorp='Maximum incorporation',dataset='Detected dataset')
   color_part <- tagList(tags$div(class='legend-title',paste0('Node colour — ',color_label)),legend_rows(cats$group,cats$color))
   size_part <- if(input$node_size=='uniform') tagList(tags$div(class='legend-title','Node size'),tags$div(class='legend-row','Uniform')) else {
     vals <- if(input$node_size=='abundance') {
       ids<-n$id; daily %>% filter(UniProtID %in% ids) %>% group_by(UniProtID) %>% summarise(v=mean_ok(log2Abundance),.groups='drop') %>% pull(v)
     } else summary %>% filter(UniProtID %in% n$id) %>% pull(MaxIncorp)
     vals<-vals[is.finite(vals)]
     if(length(vals)) {q<-as.numeric(quantile(vals,c(.1,.5,.9),na.rm=TRUE,names=FALSE)); labsz<-c(12,20,30); fmt<-if(input$node_size=='incorp') function(x) percent(x,accuracy=1) else function(x) sprintf('%.1f log2 intensity',x); tagList(tags$div(class='legend-title',paste0('Node size — ',ifelse(input$node_size=='incorp','Maximum incorporation','Mean 5 day abundance'))),size_legend_rows(q,labsz,fmt))} else NULL
   }
   edge_txt <- switch(input$network_type,correlation='pSILAC trajectory correlation; width ∝ |r|',string='STRING functional association; width ∝ STRING confidence',combined='Both pSILAC correlation + STRING; width ∝ STRING confidence')
   settings <- if(input$network_type=='string') paste0('STRING ≥ ',sprintf('%.2f',input$string_score)) else paste0(ifelse(input$corr_metric=='Incorporation','Incorporation','Abundance'),' · ',tools::toTitleCase(input$corr_method),' · |r| ≥ ',sprintf('%.2f',input$corr_threshold),' · ≥',input$corr_minpoints,' shared days',if(input$network_type=='combined') paste0(' · STRING ≥ ',sprintf('%.2f',input$string_score)) else '')
   div(class='network-legend',tags$b('Legend'),
     tags$div(class='legend-section',color_part),
     tags$div(class='legend-section',size_part),
     tags$div(class='legend-section',tags$div(class='legend-title','Special node'),HTML('<span class="legend-focal"></span>Protein of Interest')),
     tags$div(class='legend-section',tags$div(class='legend-title','Edge'),HTML('<span class="legend-edge" style="border-top-width:3px"></span>'),edge_txt),
     tags$div(class='legend-section',style='color:#667085;',paste0('Current filters: ',settings)))
 })

 # -------------------------------------------------------------------------
 # Exploratory protein-group network
 # -------------------------------------------------------------------------
 observe({
   common <- c('All','Mitochondrial','Nuclear','Cytoplasmic','Membrane','Extracellular','Golgi','Endoplasmic reticulum')
   updateSelectizeInput(session,'group_localization',choices=setNames(common,common),selected=isolate(input$group_localization),server=TRUE)
 })

 # Membership uses ALL available location fields, not the first nonmissing one.
 group_universe <- reactive({
   d <- summary
   fields <- intersect(c('UniProtLocation','GOCellularComponent'),names(d))
   if(length(fields)) {
     d$LocationSearch <- apply(d[,fields,drop=FALSE],1,function(v) paste(v[!is.na(v)],collapse='; '))
   } else d$LocationSearch <- ''
   loc <- input$group_localization
   pattern <- switch(loc,
     Mitochondrial='mitochondri|mitochon|\\bMICOS\\b', Nuclear='nucle|chromatin',
     Cytoplasmic='cytoplas|cytosol', Membrane='membrane',
     Extracellular='extracellular|secreted', Golgi='golgi',
     `Endoplasmic reticulum`='endoplasmic reticulum',
     NULL)
   if(!is.null(loc) && nzchar(loc) && loc!='All') {
     if(is.null(pattern)) pattern <- loc
     d <- d %>% filter(str_detect(LocationSearch,regex(pattern,ignore_case=TRUE)))
   }
   fd <- input$group_firstday
   if(fd %in% as.character(1:5)) d <- d %>% filter(FirstDetectedDay==as.numeric(fd))
   if(fd=='2plus') d <- d %>% filter(is.finite(FirstDetectedDay),FirstDetectedDay>=2)
   if(fd=='3plus') d <- d %>% filter(is.finite(FirstDetectedDay),FirstDetectedDay>=3)
   if(input$group_turnover!='any') d <- d %>% filter(TurnoverClass==input$group_turnover)
   mat <- if(input$group_metric=='Incorporation') incorporation_mat else abundance_mat
   # Safe accession lookup: summary may include proteins absent from the
   # selected 5 day matrix. Missing matrix rows mean no measured days.
   idx <- match(as.character(d$UniProtID), rownames(mat))
   d$n_measured <- integer(nrow(d))
   present <- which(!is.na(idx))
   if(length(present)) {
     observed <- mat[idx[present], , drop=FALSE]
     d$n_measured[present] <- rowSums(is.finite(observed))
   }
   d
 })
 group_candidates <- reactive({
   d <- group_universe() %>% filter(n_measured>=input$group_minpoints) %>%
     arrange(desc(n_measured),desc(MaxIncorp))
   if(nrow(d)>input$group_maxproteins) d <- head(d,input$group_maxproteins)
   d
 })
 output$group_diagnostics <- renderTable({
   d <- group_universe()
   data.frame(`Group proteins`=nrow(d),
     `With trajectory measurements`=sum(d$n_measured>0),
     `With 3+ days`=sum(d$n_measured>=3),
     `Passing selected day threshold`=sum(d$n_measured>=input$group_minpoints),
     check.names=FALSE)
 },striped=TRUE,bordered=FALSE,spacing='xs')

 group_edges_raw <- reactive({
   d <- group_candidates(); empty <- tibble(from=character(),to=character(),correlation=double(),shared_days=integer())
   if(nrow(d)<2) return(empty)
   mat0 <- if(input$group_metric=='Incorporation') incorporation_mat else abundance_mat
   mat <- mat0[d$UniProtID,,drop=FALSE]
   out <- vector('list',nrow(mat)-1); k <- 0L
   for(i in seq_len(nrow(mat)-1)) {
     js <- (i+1):nrow(mat); xi <- mat[i,]
     vals <- lapply(js,function(j){
       ok <- is.finite(xi) & is.finite(mat[j,]); n <- sum(ok)
       if(n<input$group_minpoints) return(NULL)
       r <- suppressWarnings(cor(xi[ok],mat[j,ok],method=input$group_method))
       if(!is.finite(r) || abs(r)<input$group_r || (input$group_sign=='positive' && r<=0)) return(NULL)
       if(isTRUE(input$group_match_turnover)) {
         h1 <- d$half_life[match(rownames(mat)[i],d$UniProtID)]; h2 <- d$half_life[match(rownames(mat)[j],d$UniProtID)]
         if(!is.finite(h1) || !is.finite(h2) || abs(h1-h2)>input$group_halfdiff) return(NULL)
       }
       tibble(from=rownames(mat)[i],to=rownames(mat)[j],correlation=r,shared_days=n)
     })
     z <- bind_rows(vals); if(nrow(z)){k<-k+1L; out[[k]]<-z}
   }
   e <- bind_rows(out[seq_len(k)])
   if(!nrow(e)) return(empty)
   e %>% arrange(desc(abs(correlation)))
 })

 group_graph_data <- reactive({
   e <- group_edges_raw()
   if(nrow(e)>input$group_maxedges) e <- head(e,input$group_maxedges)
   ids <- unique(c(e$from,e$to))
   if(!length(ids)) return(list(edges=e,nodes=tibble(),modules=tibble()))
   g <- igraph::graph_from_data_frame(e %>% select(from,to),directed=FALSE,vertices=ids)
   comp <- igraph::components(g)$membership
   module <- tibble(UniProtID=names(comp),Module=paste0('Module ',as.integer(comp)))
   d <- summary %>% filter(UniProtID %in% ids) %>%
     mutate(Localization=coalesce(UniProtLocation,GOCellularComponent,'Unknown')) %>% left_join(module,by='UniProtID')
   mean_ab <- daily %>% filter(UniProtID %in% ids) %>% group_by(UniProtID) %>% summarise(MeanAbundance=mean_ok(log2Abundance),.groups='drop')
   d <- left_join(d,mean_ab,by='UniProtID')
   groups <- switch(input$group_node_color,module=d$Module,turnover=d$TurnoverClass,
     firstday=ifelse(is.finite(d$FirstDetectedDay),paste('Day',d$FirstDetectedDay),'NA'),
     localization=location_class(d$Localization),maxincorp=ifelse(is.finite(d$MaxIncorp),paste0(round(d$MaxIncorp*10)*10,'%'),'NA'))
   sizes <- rep(20,nrow(d))
   if(input$group_node_size=='abundance' && any(is.finite(d$MeanAbundance))) {x<-d$MeanAbundance;x[!is.finite(x)]<-min(x[is.finite(x)]);if(length(unique(x))>1)sizes<-rescale(x,to=c(13,34))}
   if(input$group_node_size=='incorp' && any(is.finite(d$MaxIncorp))) {x<-d$MaxIncorp;x[!is.finite(x)]<-0;if(length(unique(x))>1)sizes<-rescale(x,to=c(13,34))}
   node_cols <- network_palette(groups,input$group_node_color)
   nodes <- d %>% transmute(id=UniProtID,label=coalesce(GeneSymbol,UniProtID),group=as.character(groups),color=node_cols,value=sizes,
     title=paste0('<b>',coalesce(GeneSymbol,UniProtID),'</b><br>',coalesce(ProteinName,''),'<br>UniProt: ',UniProtID,
       '<br>Module: ',Module,'<br>First detected: ',ifelse(is.finite(FirstDetectedDay),paste0('Day ',FirstDetectedDay),'NA'),
       '<br>Turnover: ',TurnoverClass,'<br>Half-life: ',ifelse(is.finite(half_life),paste0(round(half_life,2),' d'),'NA'),
       '<br>Max incorporation: ',ifelse(is.finite(MaxIncorp),percent(MaxIncorp,accuracy=.1),'NA'),'<br>Localization: ',Localization))
   mods <- nodes %>% transmute(UniProtID=id,Module=d$Module) %>% count(Module,name='Proteins') %>% arrange(desc(Proteins))
   list(edges=e,nodes=nodes,modules=mods)
 })

 output$group_status <- renderText({
   d<-group_candidates(); z<-group_graph_data(); loc<-input$group_localization
   paste0(ifelse(loc=='All','All localizations',loc),' · ',nrow(d),' eligible proteins · ',nrow(z$nodes),' connected proteins · ',nrow(z$edges),' edges · ',nrow(z$modules),' modules')
 })
 output$group_network <- renderVisNetwork({
   z<-group_graph_data(); validate(need(nrow(z$nodes)>0,'No correlated protein pairs pass the current filters. Try lowering minimum |r| or broadening the protein group.'))
   ve<-z$edges %>% transmute(from,to,value=pmax(abs(correlation),.01)*5,title=paste0('r = ',round(correlation,3),' · shared days = ',shared_days),color=ifelse(correlation<0,'#D97736','#168579'))
   visNetwork(z$nodes,ve,width='100%',height='610px') %>% visNodes(shape='dot',font=list(size=16)) %>%
     visEdges(smooth=FALSE) %>% visOptions(highlightNearest=list(enabled=TRUE,degree=1,hover=TRUE),nodesIdSelection=TRUE) %>%
     visPhysics(stabilization=TRUE,solver='forceAtlas2Based')
 })
 output$group_legend <- renderUI({
   z<-group_graph_data(); req(nrow(z$nodes)>0); n<-z$nodes
   cats<-n %>% distinct(group,color) %>% filter(!is.na(group),nzchar(group)) %>% arrange(group)
   color_label<-switch(input$group_node_color,module='Connected module',turnover='Turnover class',firstday='First detected day',localization='Localization',maxincorp='Maximum incorporation')
   color_part<-tagList(tags$div(class='legend-title',paste0('Node colour — ',color_label)),legend_rows(cats$group,cats$color))
   size_part<-if(input$group_node_size=='uniform') tagList(tags$div(class='legend-title','Node size'),tags$div(class='legend-row','Uniform')) else {
     vals<-if(input$group_node_size=='abundance') daily %>% filter(UniProtID %in% n$id) %>% group_by(UniProtID) %>% summarise(v=mean_ok(log2Abundance),.groups='drop') %>% pull(v) else summary %>% filter(UniProtID %in% n$id) %>% pull(MaxIncorp)
     vals<-vals[is.finite(vals)]
     if(length(vals)){q<-as.numeric(quantile(vals,c(.1,.5,.9),na.rm=TRUE,names=FALSE));fmt<-if(input$group_node_size=='incorp') function(x) percent(x,accuracy=1) else function(x) sprintf('%.1f log2 intensity',x);tagList(tags$div(class='legend-title',paste0('Node size — ',ifelse(input$group_node_size=='incorp','Maximum incorporation','Mean 5 day abundance'))),size_legend_rows(q,c(12,20,30),fmt))} else NULL
   }
   settings<-paste0(ifelse(input$group_metric=='Incorporation','Incorporation','Abundance'),' · ',tools::toTitleCase(input$group_method),' · |r| ≥ ',sprintf('%.2f',input$group_r),' · ≥',input$group_minpoints,' shared days',if(isTRUE(input$group_match_turnover)) paste0(' · Δ half-life ≤ ',input$group_halfdiff,' d') else '')
   div(class='network-legend',tags$b('Legend'),tags$div(class='legend-section',color_part),tags$div(class='legend-section',size_part),
     tags$div(class='legend-section',tags$div(class='legend-title','Edge'),HTML('<div class="legend-row"><span class="legend-edge" style="border-top-color:#168579;border-top-width:3px"></span>Positive correlation (r &gt; 0)</div><div class="legend-row"><span class="legend-edge" style="border-top-color:#D97736;border-top-width:3px"></span>Negative correlation (r &lt; 0; shown when enabled)</div><div>Edge width ∝ |r| (correlation strength).</div>')),
     tags$div(class='legend-section',style='color:#667085;',paste0('Current filters: ',settings)))
 })
 output$group_node_info <- renderUI({
   id<-input$group_network_selected
   if(is.null(id)||!length(id)) return(p(class='smallnote','Click a protein node to inspect it.'))
   id<-as.character(id[1]); a<-summary %>% filter(UniProtID==id) %>% dplyr::slice(1); if(!nrow(a)) return(NULL)
   div(h4(coalesce(a$GeneSymbol,a$UniProtID)),p(a$ProteinName),p(strong('UniProt: '),a$UniProtID,' · ',strong('Detected in: '),a$DetectedIn),
     p(strong('First detected: '),fmtday(a$FirstDetectedDay),' · ',strong('Half-life: '),ifelse(is.finite(a$half_life), sprintf('%.2f', a$half_life), 'NA'),' d · ',strong('Turnover: '),a$TurnoverClass),
     actionButton('open_group_node','Open in Protein Explorer'))
 })
 observeEvent(input$open_group_node,{id<-input$group_network_selected;req(length(id)>=1);id<-as.character(id[1]);active_protein(id);a<-ann%>%filter(UniProtID==id)%>%dplyr::slice(1);if(nrow(a)){updateTextInput(session,'protein_search',value=coalesce(a$GeneSymbol,id));updateSelectInput(session,'protein_hits',choices=setNames(a$UniProtID,ifelse(is.na(a$GeneSymbol) | !nzchar(trimws(a$GeneSymbol)),a$UniProtID,a$GeneSymbol)),selected=id)};updateNavbarPage(session,'nav',selected='Protein Explorer')})
 output$group_modules <- renderDT({z<-group_graph_data();validate(need(nrow(z$modules)>0,'No modules.'));datatable(z$modules,rownames=FALSE,options=list(pageLength=8,dom='tip'))})
 output$group_edges <- renderDT({z<-group_graph_data();validate(need(nrow(z$edges)>0,'No correlated pairs.'));labs<-ann%>%select(UniProtID,GeneSymbol)%>%distinct(UniProtID,.keep_all=TRUE);x<-z$edges%>%left_join(labs,by=c('from'='UniProtID'))%>%dplyr::rename(Protein_A=GeneSymbol)%>%left_join(labs,by=c('to'='UniProtID'))%>%dplyr::rename(Protein_B=GeneSymbol)%>%transmute(Protein_A=coalesce(Protein_A,from),Protein_B=coalesce(Protein_B,to),Correlation=round(correlation,3),`Shared days`=shared_days);datatable(x,rownames=FALSE,options=list(pageLength=10,scrollX=FALSE))})

}
shinyApp(ui,server)
