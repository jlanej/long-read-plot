print("loading libraries")
print("loading optparse")
library(optparse)
print("loading crayon")
if(!require(crayon)) install.packages("crayon")
library(crayon)

source("bamRUtils.R")
print("hi")