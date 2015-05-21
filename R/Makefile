# Makefile for generating .Rmd file by inserting code snippets from external .R files
#
# Author: dkulp 5/18/2015
#
# Just run "make" to generate the .Rmd file.
# Run "make tutorial.html" to run knitr and pandoc from command line
# Run "make test" to execute each code snippet independently.

PERL = perl
RSCRIPT = rscript
R = R
RBATCH = $(R) CMD BATCH

SRC=tutorial.Rmdp
RMD=$(SRC:.Rmdp=.Rmd)
HTML=$(SRC:.Rmdp=.html)

RCODE=code0.r code1.r code2.r code3.r code4.r code5.r code6.r code7.r code8.r code9.r
RTEST=$(RCODE:.r=.r.Rout)

.PHONY: all
all:	$(RMD)
	@echo "Now process $(RMD) in Rstudio or run 'make html' to approximate Rstudio's output"

%.Rmd : %.Rmdp
	$(PERL) rmd-inline.pl $< > $@

%.md : %.Rmd
	$(RSCRIPT) --vanilla -e 'library(knitr); knit(commandArgs(TRUE)[1])' $<

%.html : %.md
	/opt/homebrew-cask/Caskroom/rstudio/0.98.1102/RStudio.app/Contents/MacOS/pandoc/pandoc \
	 $< --to html --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash-implicit_figures \
	--output $@ --smart --email-obfuscation none --self-contained --standalone --section-divs \
	--template /usr/local/Cellar/r/3.1.2_1/R.framework/Versions/3.1/Resources/library/rmarkdown/rmd/h/default.html \
	--variable 'theme:bootstrap' --mathjax --variable 'mathjax-url:https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' \
	--no-highlight --variable highlightjs=/usr/local/Cellar/r/3.1.2_1/R.framework/Versions/3.1/Resources/library/rmarkdown/rmd/h/highlight 

%.r.Rout : %.r
	$(RBATCH) $<

clean:
	rm $(RMD) $(HTML)


test: $(RTEST)

