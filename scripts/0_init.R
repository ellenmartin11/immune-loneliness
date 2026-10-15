

install.packages("renv")
library(renv)
source("renv/activate.R")
renv::restore() #restores my env with all my packages



library(data.table)
library(dplyr)
library(biomaRt)
library(httr2)
library(jsonlite)
library(ieugwasr)
library(MungeSumstats)
