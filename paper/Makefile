#
# Makefile for plots $Id: Makefile,v 1.11 2013/12/16 20:11:39 spoel Exp $
#

EPLOTS  = sasa.pdf rmsd.pdf hb.pdf msd.pdf rmsf.pdf diff_cyto_over_singles.pdf #rg.pdf 
PLOTS = rdf_RNA_metabolites.pdf avclust.pdf rdf_mg.pdf
SPLOTS  = fig.pdf
TABS    = tbl.tex # tbl_liquids_corr.tex 
MKEPS	= xmgrace -hardcopy -hdevice EPS -printfile 
MKPDF	= xmgrace -hardcopy -hdevice PDF -printfile 
MKPNG	= xmgrace -hardcopy -hdevice PNG -printfile
XVG = 	xvg
DAT=	*.dat
DB	= $(LIQUIDS)/DATABASE/molecules.sqlite3.dat 
MYVI	= ../bibtex/myvi
RM	= /bin/rm -f
PDF	= virtual_cytoplasm_v2.pdf $(EPLOTS)  #virtual_cytoplasm-supporting.pdf

all:	$(PDF)

virtual_cytoplasm_v2.pdf: virtual_cytoplasm_v2.tex $(EPLOTS) $(PLOTS)  #tbl.tex
			pdflatex virtual_cytoplasm_v2
			pdflatex virtual_cytoplasm_v2
			bibtex virtual_cytoplasm_v2
			pdflatex virtual_cytoplasm_v2
virtual_cytoplasm-supporting.pdf:	virtual_cytoplasm-supporting.tex $(SPLOTS)


sasa.pdf: $(DAT) barplots.py
	python3 barplots.py -prop sasa

rmsd.pdf: $(DAT) barplots.py
	python3 barplots.py -prop rmsd

hb.pdf: $(DAT) barplots.py
	python3 barplots.py -prop hb

rg.pdf: $(DAT) barplots.py
	python3 barplots.py -prop rg

msd.pdf: $(DAT) barplots.py
	python3 barplots.py -prop msd

diff_cyto_over_singles.pdf: $(DAT) diff_plots.py
	python3 diff_plots.py	


%.pdf:	%.tex
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

%.dvi:	%.tex
	latex $<
	bibtex $*
	latex $<
	latex $<
	latex $<

%.ps:	%.dvi
	dvips -t letter -o $@ $<

%.pdf:	%.eps
	ps2pdf $<

%.pdf:  %.xvg
	$(MKPDF) $@ $<

# msd1.pdf:  msd1.xvg msd.par
# 	$(MKPDF) $@ $< -param msd.par
# msd2.pdf:  msd2.xvg msd.par
# 	$(MKPDF) $@ $< -param msd.par
# msd4.pdf:  msd4.xvg msd.par
# 	$(MKPDF) $@ $< -param msd.par
# msd8.pdf:  msd8.xvg msd.par
# 	$(MKPDF) $@ $< -param msd.par
# msd64.pdf:  msd64.xvg msd.par
# 	$(MKPDF) $@ $< -param msd.par

%.eps:	%.xvg
	$(MKEPS) $@ $<
	chmod +x $(MYVI)
	$(MYVI) $@
%.png:	%.xvg
	$(MKPNG) $@ $<


clean:
	$(RM) *.o  *.log *.aux *.bbl *.blg *~ $(PDF) *.dvi *.eps *.lot *.lof *.toc  $(TABS) $(EPLOTS)

