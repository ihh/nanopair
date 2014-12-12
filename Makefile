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


HDF5 = /usr/local/hdf5

TEST5 = ../../nanopore/loman/LomanLabz_PC_E.coli_MG1655_ONI_3058_1_ch101_file20_strand.fast5

ifdef HDF5
CPPFLAGS += -I${HDF5}/include -W -Wall -ansi
LDFLAGS += -L${HDF5}/lib -Wl,-rpath -Wl,${HDF5}/lib
endif

LIBS += -lhdf5_hl -lhdf5
CFLAGS += -O2 -g

SRCFILES = $(wildcard src/*.c)


bin:
	mkdir $@

bin/dump_fast5events: $(SRCFILES) Makefile bin
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} $(SRCFILES) -o $@ ${LIBS}

dump:
	h5dump -d /Analyses/Basecall_2D_000/BaseCalled_template/Events $(TEST5) | less

ptdump:
	poretools events $(TEST5) | less

test: bin/dump_fast5events
	bin/dump_fast5events $(TEST5)
