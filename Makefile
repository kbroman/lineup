all: doc docs/lineup.html

# build package documentation
doc:
	R -e 'devtools::document()'

docs/lineup.html: vignettes/lineup.Rmd docs/badges.html docs/paste_badges.R
	cd $(<D);R -e "rmarkdown::render('$(<F)')"
	mv $(<D)/$(@F) $@
	cd $(@D);paste_badges.R $(@F)
