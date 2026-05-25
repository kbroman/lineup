all: doc docs/lineup.html

# build package documentation
doc:
	R -e 'devtools::document()'

docs/lineup.html: vignettes/lineup.Rmd
	cd $(<D);R -e "rmarkdown::render('$(<F)')"
	mv $(<D)/$(@F) $@
