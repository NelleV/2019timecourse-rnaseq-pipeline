library(rmarkdown)

args = commandArgs(trailingOnly=TRUE)
filename = args[1]
outname = args[2]

render(filename, output_format="html_document", output_file=outname)
