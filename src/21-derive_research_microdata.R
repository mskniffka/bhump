# Derive microdata in event-history format for analysis of
# fetal-infant transition

# Init ------------------------------------------------------------

library(data.table)
library(lubridate)
library(magrittr)
library(readxl)
library(tidyverse)

# use all available CPUs
setDTthreads(0)

paths <- list()
paths$input <- list(
  fetoinfant = 'tmp/20-fetoinfant.RData'
)
paths$output <- list(
  fetoinfant = 'tmp/21-fetoinfant.RData'
)

# constants
cnst <-
  list(
    left_truncation_gestage = 24,
    right_censoring_gestage = 76.99
  )

# Load data -------------------------------------------------------

# US fetal- infant deaths and births individual level data
load(paths$input$fetoinfant)
setDT(fetoinfant)

# Subset ----------------------------------------------------------

# only consider cases which survived to
# week 24, i.e. were delivered (alive or dead)
# at week 24 or later and were conceived during
# the years 1989, 1999, 2009, and 2014.
# to that end we only consider deliveries during
# 1989 & 1990 -> 1989,
# 1999 & 2000 -> 1999,
# 2009 & 2010 -> 2009
# 2014 & 2015 -> 2014
fetoinfant <-
  fetoinfant[
    gestation_at_delivery_w >= cnst$left_truncation_gestage &
      date_of_delivery_y %in%
      c(1989, 1990, 1999, 2000, 2009, 2010, 2014, 2015),
    ]

# Calculate date of conception ------------------------------------

# calculate date of conception as
# date of delivery - (weeks of gestation at delivery - 2 weeks)
# the second term is the ferilization age at delivery which
# is shorter than the gestational age as fertilization on average
# happens 2 weeks following the last period
fetoinfant[
  ,
  date_of_delivery_ym :=
    paste(date_of_delivery_y,
          date_of_delivery_m,
          '15', sep = '-') %>%
    ymd()
  ]
fetoinfant[
  ,
  date_of_conception_ym :=
    (date_of_delivery_ym - weeks(gestation_at_delivery_w - 2)) %>%
    round_date(unit = 'month')
  ]
fetoinfant[
  ,
  date_of_conception_y :=
    year(date_of_conception_ym)
  ]

# subset to conception cohorts 1989, 1999, 2009, and 2014
fetoinfant <-
  fetoinfant[
    date_of_conception_y %in% c(1989, 1999, 2009, 2014)
    ]

# delete space after icd code

fetoinfant <-
  fetoinfant[, 
             cod_icd10 := gsub(" ", "", cod_icd10)]

# Recode cause of death categories --------------------------------

cod_codes <- read_xlsx("./dat/10-cod-list/cod.xlsx")[-c(3,4,5,6)] 


# cod_codes <- list(
#   Maternal = c('P00', 'P01'),
#   PCML = c('P02', 'P03', 'P04', 'P10', 'P11',
#            'P12', 'P13', 'P14', 'P15'),
#   Prematurity = c('P05', 'P07', 'P08'),
#   Respiratory = c('P20', 'P21', 'P22', 'P23', 'P24', 'P25',
#                   'P26', 'P27', 'P28'),
#   External = c('P35','P36','P37','P38','P39',
#                paste0('A', formatC(0:99, width = 2, flag = '0')),
#                paste0('B', formatC(0:99, width = 2, flag = '0')),
#                paste0('V', formatC(0:99, width = 2, flag = '0')),
#                paste0('W', formatC(0:99, width = 2, flag = '0')),
#                paste0('X', formatC(c(0:59, 85:99), width = 2, flag = '0')),
#                paste0('Y', formatC(c(0:9, 10:36, 40:84), width = 2, flag = '0'))
#   ),
#   Unspecific = c('P95', paste0('R', formatC(0:99, width = 2, flag = '0'))),
#   Neoplasms = c(paste0('Q', formatC(0:99, width = 2, flag = '0')),
#                 paste0('C', formatC(0:99, width = 2, flag = '0')),
#                 paste0('D', formatC(0:49, width = 2, flag = '0')))
# )


fetoinfant <- merge(fetoinfant, cod_codes, by = "cod_icd10", all = T)
fetoinfant <- fetoinfant[ !is.na(type)] 

fetoinfant <- fetoinfant %>% 
  mutate(cod_cat = case_when(
    (date_of_delivery_y == 2014 & type == "fetus" & cod_icd10 == "") ~ "Other",
    (date_of_delivery_y == 2015 & type == "fetus" & cod_icd10 == "") ~ "Other",
    (cod_cat == "Maternal Complications" & age_at_death_d >= 100) ~ "Other",
    (cod_cat == "Labour, Cord, Membrane and Labour Complications" & age_at_death_d >= 100) ~ "Other",
    TRUE ~ cod_cat
  ))



fetoinfant[
  ,
  cod_cat := fcase(
    # substr(cod_icd10,1,3)%in%cod_codes$Maternal, 'Maternal',
    # substr(cod_icd10,1,3)%in%cod_codes$PCML, 'PCML',
    # substr(cod_icd10,1,3)%in%cod_codes$Prematurity, 'Prematurity',
    # substr(cod_icd10,1,3)%in%cod_codes$InfectionsParacitesOperations, 'InfectionsParacitesOperations',
    # substr(cod_icd10,1,3)%in%cod_codes$ViolenceAccidents, 'ViolenceAccidents',
    # substr(cod_icd10,1,3)%in%cod_codes$UnspecificStillbirth, 'UnspecificStillbirth',
    # substr(cod_icd10,1,3)%in%cod_codes$SuddenInfantDeath, 'SuddenInfantDeath',
    # substr(cod_icd10,1,3)%in%cod_codes$TreatableNeoplasms, 'TreatableNeoplasms',
    # substr(cod_icd10,1,3)%in%cod_codes$UntreatableNeoplasms, 'UntreatableNeoplasms',
    !is.na(cod_icd10), cod_cat,
    default = NA
  )
]


# Add flags for vital events --------------------------------------

# add flags for vital events fetal-death, life-birth,
# neonatal and post-neonatal death and survival
fetoinfant[
  ,
  `:=`(
    # fetal deaths are all cases from the
    # fetal death file
    fetal_death =
      type == 'fetus',
    # life births are all cases from the
    # birth registry file
    life_birth =
      type == 'infant',
    # neonatal deaths are all cases with an entry
    # for the age at death in days which is < 7
    # FALSE & NA = FALSE
    neonatal_death =
      (!is.na(age_at_death_d)) & (age_at_death_d < 7),
    # neonatal survivors are all cases from the infant file
    # without an age at death entry or with age at death >= 7
    neonatal_survivor =
      type == 'infant' & (is.na(age_at_death_d) | (age_at_death_d >= 7)),
    # postneonatal deaths are all cases with an
    # age at death >= 7
    postneonatal_death =
      (!is.na(age_at_death_d)) & (age_at_death_d >= 7),
    # postneonatal survivors are all cases from the infant file
    # without an age at death
    postneonatal_survivor =
      type == 'infant' & is.na(age_at_death_d)
  )
]

# Add timing of vital events --------------------------------------

# assume uniform distribution of events over week of
# gestation for fetal deaths and life-births
fetoinfant[
  ,
  `:=`(
    # gestational age at fetal death in weeks
    # assume uniform distribution of fetal deaths
    # within week
    gestage_at_fetal_death_w =
      fifelse(
        fetal_death,
        gestation_at_delivery_w + runif(.N),
        as.numeric(NA)
      ),
    # gestational age at life-birth
    # assume uniform distribution of life-births
    # within week
    gestage_at_life_birth_w =
      fifelse(
        life_birth,
        gestation_at_delivery_w + runif(.N),
        as.numeric(NA)
      ),
    # chronological age at infant death in (fractional) weeks
    # assume that death occurs uniformly throughout a day
    age_at_death_w =
      (age_at_death_d + runif(.N))/7
  )]


fetoinfant[
  ,
  `:=`(
    # gestational age at neonatal death in weeks
    gestage_at_neonatal_death_w =
      fifelse(
        neonatal_death,
        gestage_at_life_birth_w + age_at_death_w,
        as.numeric(NA)
      ),
    # gestational age at postneonatality in weeks
    gestage_at_postneonatality_w =
      fifelse(
        neonatal_survivor,
        gestage_at_life_birth_w + 1,
        as.numeric(NA)
      ),
    # gestation at postneonatal death in weeks
    gestage_at_postneonatal_death_w =
      fifelse(
        postneonatal_death,
        gestage_at_life_birth_w + age_at_death_w,
        as.numeric(NA)
      ),
    # censoring at end of week 76 = 1 year post viability
    gestage_at_censoring_w =
      cnst$right_censoring_gestage
  )]

# add id variable
fetoinfant[,id := 1:.N]

# select variables
fetoinfant <-
  fetoinfant[,
          .(
            id, sex,
            race_and_hispanic_orig_of_mother,
            date_of_conception_y, date_of_conception_ym,
            fetal_death, life_birth,
            neonatal_death, neonatal_survivor,
            postneonatal_death, postneonatal_survivor,
            gestage_at_fetal_death_w,
            gestage_at_life_birth_w,
            gestage_at_neonatal_death_w,
            gestage_at_postneonatality_w,
            gestage_at_postneonatal_death_w,
            gestage_at_censoring_w,
            cod_cat
          )]

# Convert data to event-history format ----------------------------

# population strata
strata <-
  c(
    'id', 'sex',
    'race_and_hispanic_orig_of_mother',
    'date_of_conception_y'
  )

# Fetal death event histories -------------------------------------

# fetal death event histories
# fetus -> death
fetal_deaths <- fetoinfant[fetal_death == TRUE]
fetal_deaths_histories <-
  cbind(
    fetal_deaths[,..strata],
    data.table(
      origin = 'fetus',
      entry_time = cnst$left_truncation_gestage,
      destination = 'death',
      destination2 = fetal_deaths$cod_cat,
      exit_time = fetal_deaths$gestage_at_fetal_death_w
    )
  )

# Neonatal death event histories ----------------------------------

# neonatal death event histories
neonatal_deaths <- fetoinfant[neonatal_death == TRUE]
neonatal_deaths_histories <-
  rbind(
    # fetus -> neonate
    data.table(
      origin = 'fetus',
      entry_time = cnst$left_truncation_gestage,
      destination = 'neonate',
      destination2 = 'neonate',
      exit_time = neonatal_deaths$gestage_at_life_birth_w
    ),
    # neonate -> death
    data.table(
      origin = 'neonate',
      entry_time = neonatal_deaths$gestage_at_life_birth_w,
      destination = 'death',
      destination2 = neonatal_deaths$cod_cat,
      exit_time = neonatal_deaths$gestage_at_neonatal_death_w
    )
  )
# why does this work?
# the number of rows in neonatal_deaths_histories is an integer
# multiple of the number of rows in neonatal_deaths, the latter
# rows being recycled to the length of the neonatal_deaths_histories.
# as all the individuals in neonatal_deaths_histories experience the
# same number of events/rows and the order of events is encoded in the
# order of rows, the recycled id's actually belong to the same individual.
# the same logic applies to the other event_histories.
neonatal_deaths_histories <-
  cbind(
    neonatal_deaths[,..strata],
    neonatal_deaths_histories
  )

# Postneonatal death event histories ------------------------------

# postneonatal death event histories
# fetus -> neonatal -> postneonatal -> (death|censored)
postneonatal_deaths <- fetoinfant[postneonatal_death == TRUE]
postneonatal_deaths[
  ,
  censored := gestage_at_postneonatal_death_w > gestage_at_censoring_w
  ]
postneonatal_deaths_histories <-
  rbind(
    # fetus -> neonate
    data.table(
      origin = 'fetus',
      entry_time = cnst$left_truncation_gestage,
      destination = 'neonate',
      destination2 = 'neonate',
      exit_time = postneonatal_deaths$gestage_at_life_birth_w
    ),
    # neonate -> post-neonate
    data.table(
      origin = 'neonate',
      entry_time = postneonatal_deaths$gestage_at_life_birth_w,
      destination = 'postneonate',
      destination2 = 'postneonate',
      exit_time = postneonatal_deaths$gestage_at_postneonatality_w
    ),
    # post-neonate -> (death | censoring)
    # infant -> (death|censoring)
    data.table(
      origin = 'postneonate',
      entry_time = postneonatal_deaths$gestage_at_postneonatality_w,
      destination =
        fifelse(
          postneonatal_deaths$censored == FALSE,
          'death',
          'censored'
        ),
      destination2 =
        fifelse(
          postneonatal_deaths$censored == FALSE,
          postneonatal_deaths$cod_cat,
          'censored'
        ),
      exit_time =
        fifelse(
          postneonatal_deaths$censored == FALSE,
          postneonatal_deaths$gestage_at_postneonatal_death_w,
          postneonatal_deaths$gestage_at_censoring_w
        )
    )
  )
postneonatal_deaths_histories <-
  cbind(
    postneonatal_deaths[,..strata],
    postneonatal_deaths_histories
  )

# Postneonatal survivor event histories ---------------------------

# infant survivor event histories
# fetus -> neonatal -> postneonatal -> censored
postneonatal_survivors <- fetoinfant[postneonatal_survivor == TRUE]
postneonatal_survivors_histories <-
  rbind(
    # fetus -> neonate
    data.table(
      origin = 'fetus',
      entry_time = cnst$left_truncation_gestage,
      destination = 'neonate',
      destination2 = 'neonate',
      exit_time = postneonatal_survivors$gestage_at_life_birth_w
    ),
    # neonate -> post-neonate
    data.table(
      origin = 'neonate',
      entry_time = postneonatal_survivors$gestage_at_life_birth_w,
      destination = 'postneonate',
      destination2 = 'postneonate',
      exit_time = postneonatal_survivors$gestage_at_postneonatality_w
    ),
    # post-neonate -> (death | censoring)
    # infant -> (death|censoring)
    data.table(
      origin = 'postneonate',
      entry_time = postneonatal_survivors$gestage_at_postneonatality_w,
      destination = 'censored',
      destination2 = 'censored',
      exit_time = postneonatal_survivors$gestage_at_censoring_w
    )
  )
postneonatal_survivors_histories <-
  cbind(
    postneonatal_survivors[,..strata],
    postneonatal_survivors_histories
  )

# Complete feto-infant event histories ----------------------------

# event histories for each subject during
# the feto-infant period
fetoinfant_event_histories <-
  rbind(
    fetal_deaths_histories,
    neonatal_deaths_histories,
    postneonatal_deaths_histories,
    postneonatal_survivors_histories
  )
fetoinfant_event_histories <-
  fetoinfant_event_histories[order(id)]

# Consistency checks --------------------------------------------

# cohort size matches
fetoinfant[,.N] == fetoinfant_event_histories[,length(unique(id))]

# number of fetal deaths matches
fetoinfant_event_histories[
  origin == 'fetus' & destination == 'death',
  .N
  ] ==
  fetoinfant[
    fetal_death == TRUE,
    .N
    ]
# number of neonatal deaths matches
fetoinfant_event_histories[
  origin == 'neonate' & destination == 'death',
  .N
  ] ==
  fetoinfant[
    neonatal_death == TRUE,
    .N
    ]
# number of post-neonatal deaths matches
fetoinfant_event_histories[
  origin == 'postneonate' & destination == 'death',
  .N
  ] ==
  fetoinfant[
    postneonatal_death == TRUE &
      (gestage_at_postneonatal_death_w < gestage_at_censoring_w),
    .N
    ]
# number of censorings matches
fetoinfant_event_histories[
  origin == 'postneonate' & destination == 'censored',
  .N
  ] ==
  fetoinfant[
    postneonatal_survivor == TRUE |
      (gestage_at_postneonatal_death_w >= gestage_at_censoring_w),
    .N
    ]

# is entry time always smaller than exit time?
all(fetoinfant_event_histories$entry_time <
  fetoinfant_event_histories$exit_time)

# is the time of exit from a spell identical to the time
# of entry into the next spell?
fetoinfant_event_histories[
  .N > 1,
  .(same = entry_time[-1] == exit_time[-length(exit_time)]),
  by = id
  ][,
    all(same)
  ]

# save the processed microdata
save(
  fetoinfant_event_histories,
  file = paths$output$fetoinfant
)
