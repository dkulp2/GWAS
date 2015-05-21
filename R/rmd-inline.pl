#!/usr/bin/perl
#
# Author: dkulp, 5/18/2015
#
# Quick code to insert code snippets from .R files into an Rmd file.
# This serves the same purpose as code externalization using knitr:read_chunk(),
# but that is broken.  So instead, I extend the syntax where code snippets
# of the form
#
# ```{r "code0.r"/code0, eval=FALSE}
# ```
#
# are replaced with the code0 snippet in the file code0.r.  The code0 snippet
# is marked by a comment line with 4 dashes and the snippet name, i.e.
#
# ---- code0 
#
# and ended by the next occurrence of a 4 dash comment.

use FileHandle;
use strict;

while (<>) {
    if (my ($pre, $fname, $snippet, $rest) = (/^(\`\`\`\{r )\"([^\"]+)\"\/([^ ,}]+)(.*)/)) {
        my $code = "";
        my $fn = new FileHandle $fname;
        my $line;

        while ($line = <$fn> and $line !~ /^\#\s+\-{4,}\s+${snippet}/) {}

        while ($line = <$fn> and $line !~ /^\#\s+\-{4,}/) {
            $code .= $line;
        }
        print "$pre$snippet$rest\n";
        print $code;
        $fn->close();
    } else {
        print $_;
    }
}
