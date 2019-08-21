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

## make table 2, export to file 'table_2.txt' ##
source('table_2.R')

rm(list=ls())

## make figure 2, export to file 'figure_2.pdf' ##
source('figure_2.R')
