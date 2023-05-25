# Download data on US births, infant- and fetal deaths
#
# Jonas Schöley
# 2020-06-10
#
# Download data on US births, fetal and infant deaths from the web.

# Init ------------------------------------------------------------

library(yaml)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  codebook_infant = 'cfg/codebook-us_cohort_infant_births_deaths_minimal.yaml',
  codebook_fetus = 'cfg/codebook-us_fetal_deaths_minimal.yaml'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  infant = 'dat/10-nchs-us_cohort_linked_infant_deaths_births/',
  fetus = 'dat/10-nchs-us_fetal_deaths/'
)

# codebooks
codebook <- list()
codebook$infant <- read_yaml(paths$input$codebook_infant)
codebook$fetus <- read_yaml(paths$input$codebook_fetus)

# Functions -------------------------------------------------------

# download and save files
DownloadFiles <- function (url, save_path) {
  options(timeout = max(300, getOption('timeout')))
  for (i in url) {
    cat("Download", i)
    download.file(
      url = i,
      destfile =
        # save using the filename on server
        paste0(save_path, rev(strsplit(i, '/')[[1]])[1]),
      mode = 'wb'
    )
  }
}

# Download files --------------------------------------------------

files <- list()

# extract file urls from codebook
files$infant_url <- lapply(codebook$infant$files, function(x) x$url)
files$fetus_url <- lapply(codebook$fetus$files, function(x) x$url)

# download US cohort linked infant birth / death data and guides
DownloadFiles(
  url = files$infant_url, save_path = paths$output$infant
)
# DownloadFiles(
#   url = files$infantguide_url, save_path = paths$output$infant
# )

# download US period fetal death data and guides
DownloadFiles(
  url = files$fetus_url, save_path = paths$output$fetus
)
# DownloadFiles(
#   url = files$fetalguide_urls, save_path = paths$output$fetus
# )
