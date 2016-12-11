#!/bin/bash
## simplify my workflow.
cd ./OncoSimulR/vignettes

grep -v -F "Table: (\#tab:" OncoSimulR.Rmd > OncoSimulR-tex.Rmd
sed -i 's/## caption = "\\\\label/caption = "\\\\label/' OncoSimulR-tex.Rmd
sed -i 's/## panderOptions("table.split.cells", 8)/panderOptions("table.split.cells", 8)/' OncoSimulR-tex.Rmd

Rscript -e 'library(rmarkdown); library(BiocStyle); library(bookdown); render("OncoSimulR-tex.Rmd", output_format = bookdown::pdf_document2(toc = TRUE, toc_depth = 4, keep_tex = TRUE))'

cp OncoSimulR-tex.pdf OncoSimulR.pdf

cd ../..
