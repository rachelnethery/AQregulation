library(mvtnorm)
library(dplyr)
library(splines)
library(BayesTree)
library(xtable)
library(ggplot2)

## create dataset ##
source('create_dataset.R')

rm(list=ls())

## run the analysis ##
source('analysis_wholeUS.R')

rm(list=ls())

## make table 3, export to file 'table_3.txt' ##
source('table_3.R')

rm(list=ls())

## make figure 3, export to file 'figure_3.pdf' ##
source('figure_3.R')
