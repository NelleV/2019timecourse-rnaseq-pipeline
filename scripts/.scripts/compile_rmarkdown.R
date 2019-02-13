library(rmarkdown)

args = commandArgs(trailingOnly=TRUE)
filename = args[1]
render(filename)
