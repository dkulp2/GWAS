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

SRC=tut_intro.Rmdp tut_preproc.Rmdp tut_generation.Rmdp tut_analysis.Rmdp tut_post.Rmdp
RMD=$(SRC:.Rmdp=.Rmd)
HTML=$(SRC:.Rmdp=.html)
ERB=$(SRC:.Rmdp=.html.erb)

RCODE=code0.r code1.r code2.r code3.r code4.r code5.r code6.r code7.r code8.r code9.r
RTEST=$(RCODE:.r=.r.Rout)

DEST=../site/source/tut

.PHONY: all tutorial install
all:	$(RMD)
	@echo "Now process $(RMD) in Rstudio or run 'make tutorial.html' to approximate Rstudio's output"

tutorial:	$(HTML)

install:	$(ERB)
	mkdir -p $(DEST)
	cp $(ERB) $(DEST)

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

%.html : %.Rmd
	$(RSCRIPT) --vanilla -e 'library(knitr); knitr::knit2html(commandArgs(TRUE)[1])' $<

# Strip out the header tags from the HTML and remove the padding from the body.
# This 1-liner is dependent on the way that knit2html generates the HTML file.
%.html.erb : %.html

	$(PERL) -ne 'print unless /DOCTYPE|<html|html>|head>|<meta|title>|body>|padding: 1em|data:text\/css|max-width: 800px/' $< | sed -e 's/^<pre><code/<pre class="boot-clear"><code/' > $@


%.r.Rout : %.r
	$(RBATCH) $<

clean:
	rm -f $(RMD) $(HTML) $(RTEST) $(ERB)

realclean:	clean
	rm -rf cache

test: $(RTEST)

