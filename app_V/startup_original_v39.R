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

# Motif sequences and classes are read exclusively from the corrected CSV.

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
   mutate(GeneSymbol=trimws(sub('[; ,].*$', '', as.character(GeneSymbol)))) %>%
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
  # Stage is authoritative for cell number; Origin describes experimental origin.
  x <- str_to_lower(as.character(stage))
  o <- str_to_lower(as.character(origin))
  out <- case_when(
    grepl('8\\s*cell|8cell', x) ~ '8-cell',
    grepl('6\\s*cell|6cell', x) ~ '6-cell',
    grepl('4\\s*cell|4cell', x) ~ '4-cell',
    grepl('2\\s*cell|2cell', x) ~ '2-cell',
    x == 'oocyte' ~ 'Oocyte',
    o == 'oocyte' ~ 'Oocyte',
    o == 'zygote' ~ 'Zygote',
    TRUE ~ NA_character_
  )
  out
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
