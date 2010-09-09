SHELL = /bin/bash

.PHONY: clean new full

poster.pdf: poster.tex
	pdflatex ${<:.tex=}

full: poster.tex poster.bib 
	pdflatex ${<:.tex=}
	bibtex ${<:.tex=}
	pdflatex ${<:.tex=}
	pdflatex ${<:.tex=}

clean:
	rm -f poster.aux poster.log poster.nav poster.out poster.pdf poster.snm poster.toc poster.blg poster.bbl texput.log

new: clean poster.pdf
