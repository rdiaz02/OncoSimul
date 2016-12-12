#!/bin/bash
## Create html and pdf vignettes and place in github pages directory.

## cd ./OncoSimulR/vignettes
## Rscript -e 'library(rmarkdown); library(BiocStyle); render("OncoSimulR.Rmd")'
./html-vignette.sh
cd ./OncoSimulR/vignettes
mv OncoSimulR.html ../../../oncosimul-gh-pages/.
cd ../../
./pdf-vignette.sh 
## Rscript -e 'library(rmarkdown); library(BiocStyle); library(bookdown); render("OncoSimulR.Rmd", output_format = bookdown::pdf_document2(toc = TRUE, toc_depth = 4, keep_tex = TRUE))'
cd ./OncoSimulR/vignettes
mv OncoSimulR.pdf ../../../oncosimul-gh-pages/pdfs/.
rm -r -f OncoSimulR_files
rm -r -f OncoSimulR-tex_files
rm OncoSimulR.tex
rm OncoSimulR.synctex.gz
rm OncoSimulR-tex.tex
rm OncoSimulR-tex.Rmd
rm OncoSimulR-tex.pdf
rm OncoSimulR-tex.synctex.gz
rm rl1.txt
gv0=$(git rev-parse --short HEAD)
cd ../../../oncosimul-gh-pages
gv=$(git rev-parse --short HEAD)
echo "   ****   Now, Commit and push, if appropriate, including gh-pages. We are at " $gv " at " $gv0 " original"








