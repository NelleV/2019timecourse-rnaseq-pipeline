library(rmarkdown)

args = commandArgs(trailingOnly=TRUE)
filename = args[1]
outname = args[2]

render(filename, output_file=outname)
draft("MyArticle.Rmd",
      template="f1000_article",
      package="BiocWorkflowTools",
      edit=FALSE)
