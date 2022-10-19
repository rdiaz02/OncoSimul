#!/bin/bash
## simplify my workflow.


texi2pdf -q relfunct.tex
pdfcrop --noverbose relfunct.pdf relfunct-c.pdf
mv relfunct-c.pdf relfunct.pdf


## Beware that the first table, the one based on GSR, should
## not be caught by this grep. In that table there is no
## space after ":" --- could use a different convention
## and name all others as "tab:tabtab" or something.
grep -v -F "Table: (\#tab:" OncoSimulR.Rmd > OncoSimulR-tex.Rmd
sed -i 's/## caption = "\\\\label/caption = "\\\\label/' OncoSimulR-tex.Rmd
sed -i 's/## panderOptions("table.split.cells", 8)/panderOptions("table.split.cells", 8)/' OncoSimulR-tex.Rmd
sed -i 's/## panderOptions("table.split.cells", 12)/panderOptions("table.split.cells", 12)/' OncoSimulR-tex.Rmd
sed -i 's/## panderOptions("table.split.cells", 15)/panderOptions("table.split.cells", 15)/' OncoSimulR-tex.Rmd


sed -i 's/relfunct.png/relfunct.pdf/' OncoSimulR-tex.Rmd
sed -i 's/hurlbut.png/hurlbut.pdf/' OncoSimulR-tex.Rmd

## Using lualatex allows us to deal with some symbols, but messes other things
## such as a simple balbalba$\alpha$
## and if the symbols are used in a figure, then we get the PDF with
## incorrect symbols a warnings.
## So avoid it
## Using xelatex lead to problems in other places, like tables: \blandscape
## Rscript -e 'library(rmarkdown); library(BiocStyle); library(bookdown); render("OncoSimulR-tex.Rmd", output_format = bookdown::pdf_document2(toc = TRUE, toc_depth = 4, keep_tex = TRUE, latex_engine = "lualatex"))'
## Rscript -e 'library(rmarkdown); library(BiocStyle); library(bookdown); render("OncoSimulR-tex.Rmd", output_format = bookdown::pdf_document2(toc = TRUE, toc_depth = 4, keep_tex = TRUE, latex_engine = "xelatex"))'
Rscript -e 'library(rmarkdown); library(BiocStyle); library(bookdown); render("OncoSimulR-tex.Rmd", output_format = bookdown::pdf_document2(toc = TRUE, toc_depth = 4, keep_tex = TRUE))'

cp OncoSimulR-tex.pdf OncoSimulR.pdf


