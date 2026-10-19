main.pdf: main.tex sections/*.tex tables/assessment/*.tex references.bib
	pdflatex -interaction=nonstopmode main.tex
	bibtex main
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex

clean:
	rm -f *.aux *.bbl *.blg *.lof *.log *.lot *.out *.toc
