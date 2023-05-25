# Aggregate microdata on fetal and infant death to feto-infant lifetable

# Init ------------------------------------------------------------

library(tidyverse)

paths <- list()
paths$input <- list(
  fetoinfant = 'tmp/21-fetoinfant.RData',
  lifetable_functions = 'src/00-fnct-feto_infant_lt.R'
)
paths$output <- list(
  fetoinfant_lt = 'out/30-fetoinfant_lifetables.RData'
)

# multistate aggregation of events and exposure times
source(paths$input$lifetable_functions)

cnst <- list()
# aggregate over each single week LMP from 24 to 77
cnst$lifetable_breaks <- 24:77

# Load data -------------------------------------------------------

# individual level event history data on feto-infant survival
load(paths$input$fetoinfant)

# Aggregate observations into multistate lt -----------------------

filt <- list()

# total (cohort 2009)
filt$total09 <-
  FILT(
    fetoinfant_event_histories %>%
      filter(date_of_conception_y == 2009),
    breaks = cnst$lifetable_breaks
  )

# by sex (cohort 2009)
filt$sex09 <-
  FILT(
    fetoinfant_event_histories %>%
      filter(date_of_conception_y == 2009),
    breaks = cnst$lifetable_breaks,
    stratum = sex
  )

# by cohort
filt$cohort <-
  FILT(
    fetoinfant_event_histories,
    breaks = cnst$lifetable_breaks,
    stratum = date_of_conception_y
  )

# by origin of mother (cohort 2009)
filt$origin09 <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2009) %>%
  mutate(
    maternal_origin = recode_factor(
      race_and_hispanic_orig_of_mother,
      'Mexican' = 'Hispanic',
      'Puerto Rican' = 'Hispanic',
      'Cuban' = 'Hispanic',
      'Central or South American' = 'Hispanic',
      'Other and unknown Hispanic' = 'Hispanic',
      'Non-Hispanic other races' = 'Other',
    ) %>%
      forcats::fct_explicit_na()
  ) %>%
  filter(
    maternal_origin %in%
      c('Hispanic', 'Non-Hispanic Black', 'Non-Hispanic White')
  ) %>%
  FILT(
    df = .,
    breaks = cnst$lifetable_breaks,
    stratum = maternal_origin
  )

# Cause of death --------------------------------------------------

# cohort 2014

filt$maternal <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    d_out = 'destination2', death_state_name = 'Maternal',
    breaks = cnst$lifetable_breaks
  )
filt$neoplasm <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    d_out = 'destination2', death_state_name = 'Neoplasms',
    breaks = cnst$lifetable_breaks
  )
filt$external <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    d_out = 'destination2', death_state_name = 'External',
    breaks = cnst$lifetable_breaks
  )
filt$PCML <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    d_out = 'destination2', death_state_name = 'PCML',
    breaks = cnst$lifetable_breaks
  )
filt$Prematurity <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    d_out = 'destination2', death_state_name = 'Prematurity',
    breaks = cnst$lifetable_breaks
  )
filt$Other <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    d_out = 'destination2', death_state_name = 'Other',
    breaks = cnst$lifetable_breaks
  )
filt$Respiratory <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    d_out = 'destination2', death_state_name = 'Respiratory',
    breaks = cnst$lifetable_breaks
  )
filt$Unspecific <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    d_out = 'destination2', death_state_name = 'Unspecific',
    breaks = cnst$lifetable_breaks
  )

# Export ----------------------------------------------------------

save(filt, file = paths$output$fetoinfant_lt)
