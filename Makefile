LATEX := pdflatex
BIBTEX := bibtex

PRODUCTION := paper

open: $(addsuffix .pdf.open,$(PRODUCTION))

all: $(addsuffix .pdf,$(PRODUCTION)) latex.zip

up:
	cvs up -dPR

TEXFILES := paper.tex paper.bbl
latex.zip: $(TEXFILES)
	zip $@ $(TEXFILES)

%.clean:
	rm $*.aux $*.bbl $*.blg $*.log $*.pdf

%.pdf: %.tex
	test -e $*.aux && rm $*.aux || eval
	$(LATEX) $*
	-$(BIBTEX) $*
	$(LATEX) $*
	$(LATEX) $*

%.open: %
	open $<

.SECONDARY:

t/refg-orthologs-mouse-human.txt:
	mysql go_latest_lite -e "SELECT xref_dbname1,xref_key1, symbol1, xref_dbname2,xref_key2, symbol2 FROM orthologous_to_J_xref WHERE species1_id IN (SELECT id FROM species WHERE ncbi_taxa_id=9606) AND  species2_id IN (SELECT id FROM species WHERE ncbi_taxa_id=10090) " > $@
.PRECIOUS: t/refg-orthologs-mouse-human.txt

t/refg-assocs-mouse.txt: t/refg-orthologs-mouse-human.txt
	blip -r go -r go_assoc_local/mgi -i $< -f "tbl(ortho)" -u ontol_db -u curation_db findall "ortho(_,_,_,A,B,_),concat_atom([A,B],':',X),curation_statement(_,X,_,T),class(T,TN)" -select "assoc(X,T,TN)" > $@.tmp && sort -u $@.tmp > $@

t/refg-assocs-human.txt: t/refg-orthologs-mouse-human.txt
	blip -r go -r go_assoc_local/goa_human -i $< -f "tbl(ortho)" -u ontol_db -u curation_db findall "ortho(A,B,_,_,_,_),concat_atom([A,B],':',X),curation_statement(_,X,_,T),class(T,TN)" -select "assoc(X,T,TN)" > $@.tmp && sort -u $@.tmp > $@
