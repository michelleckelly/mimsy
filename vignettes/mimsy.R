## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", 
                      attr.source='.numberLines'
)
library(mimsy)
library(knitr)
library(xfun)
xfun::pkg_load2(c("base64enc", "htmltools", "mime"))
library(kableExtra)

## ----install.packages, eval = FALSE-------------------------------------------
#  # Install the package
#  install.packages("mimsy")
#  
#  # Load the package
#  library(mimsy)

## ----install.developers, eval = FALSE-----------------------------------------
#  # Install and load devtools
#  install.packages("devtools")
#  library(devtools)
#  
#  # Download mimsy from Github
#  install_github("michelleckelly/mimsy")
#  
#  # Load the package
#  library(mimsy)

## ----load.csv_view, eval = FALSE----------------------------------------------
#  # Load data into R
#  data <- read.csv(file = "data_twoTemp.csv", header = TRUE, stringsAsFactors = FALSE)
#  
#  # Check it out
#  data

## ----load.csv_hide, eval = TRUE, echo = FALSE---------------------------------
# Load data into R
data <- read.csv(file = "data_twoTemp.csv", header = TRUE, stringsAsFactors = FALSE)
# Check it out
data %>%
  kable() %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "500px")

## ----run, eval = TRUE, warning=FALSE------------------------------------------
# Run the function
results <- mimsy(data, baromet.press = 977.2, units = "hPa")

## ----results, eval=FALSE------------------------------------------------------
#  # Check out the summarized results
#  results$results

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Check out the summarized results
results$results %>%
  kable() %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "500px")

## ----solub, eval=FALSE--------------------------------------------------------
#  # Check out the solubility concentrations
#  results$solubility.Concentrations

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Check out the solubility concentrations
results$solubility.Concentrations %>%
  kable() %>%
  kable_styling()

## ----save, eval = FALSE-------------------------------------------------------
#  # Save output to an Excel workbook
#  mimsy.save(results, file = "results.xlsx")
#  
#  # Save output to an RData file
#  save(results, file = "results.RData")

## ----fullScript, eval=FALSE---------------------------------------------------
#  # Install mimsy
#  install.packages("mimsy")
#  
#  # Load mimsy
#  library(mimsy)
#  
#  # Load data into R
#  data <- read.csv(file = "data.csv", header = TRUE, stringsAsFactors = FALSE)
#  
#  # Run the mimsy function
#  results <- mimsy(data, baromet.press = 977.2, units = "hPa")
#  
#  # Save the results
#  mimsy.save(results, file = "results.xlsx") # To Excel file
#  save(results, file = "results.RData") # To RData file
#  
#  # Done! :)

