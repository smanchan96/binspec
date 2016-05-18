R CMD BATCH build.R
cat build.Rout
cat NAMESPACE
R CMD Rd2pdf --pdf --title=binspec -o binspec.pdf man/*.Rd
