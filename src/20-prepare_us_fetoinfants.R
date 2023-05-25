# Read, label and concatenate data on US births, infant- and fetal deaths
#
# (1) Read NCHS data into R applying variable specifications stored
#     in custom codebook
# (2) Concatenate NCHS data across multiple years

# Init ------------------------------------------------------------

library(yaml); library(dplyr)

#memory.limit(64000)

# Constants -------------------------------------------------------

paths <- list()
paths$input <- list(
  codebook_functions = 'src/00-codebook.R',
  codebook_infant = 'cfg/codebook-us_cohort_infant_births_deaths_minimal.yaml',
  codebook_fetus = 'cfg/codebook-us_fetal_deaths_minimal.yaml',
  zip_infant = 'dat/10-nchs-us_cohort_linked_infant_deaths_births/',
  zip_fetus = 'dat/10-nchs-us_fetal_deaths/'
)
paths$output <- list(
  fetoinfant = 'tmp/20-fetoinfant.RData'
)

# codebook function
source(paths$input$codebook_functions)

# Read codebook ---------------------------------------------------

codebook <- list()

codebook$fetus <- ReadCodebook(paths$input$codebook_fetus)
codebook$infant <- ReadCodebook(paths$input$codebook_infant)

# Read data into R and apply varspecs -----------------------------

infant <- ReadFromZip(
  codebook$infant, paths$input$zip_infant, subset = c(
    'Cohort1989','Cohort1990',
    'Cohort1999','Cohort2000',
    'Cohort2009','Cohort2010',
    'Cohort2014','Cohort2015'
  ))
fetus <- ReadFromZip(
  codebook$fetus, paths$input$zip_fetus, subset = c(
    'Period1989','Period1990',
    'Period1999','Period2000',
    'Period2009','Period2010',
    'Period2014','Period2015'
  ))

# Concatenate data ------------------------------------------------

# merge data on births, fetal- and infant deaths cross years
fetoinfant <-
  bind_rows(
    infant = bind_rows(infant),
    fetus = bind_rows(fetus),
    .id = 'type'
  )

# Save ------------------------------------------------------------

save(
  fetoinfant,
  file = paths$output$fetoinfant
)
