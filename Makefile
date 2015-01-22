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
ifeq "$(wildcard $(HDF5))" ""
HDF5 = /usr/local
endif

TEST_FAST5 = nanopair-data/LomanLabz_PC_E.coli_MG1655_ONI_3058_1_ch101_file20_strand.fast5
TEST_FASTA = nanopair-data/U00096.2.fas

ifdef HDF5
CPPFLAGS += -I${HDF5}/include -W -Wall -Wno-unused-function -Wno-unused-parameter -std=c99
LDFLAGS += -L${HDF5}/lib -Wl,-rpath -Wl,${HDF5}/lib -lz
endif

LIBS += -lhdf5_hl -lhdf5
# CFLAGS += -O2 -g
CFLAGS += -g

SRCFILES = $(wildcard src/*.c)


bin:
	mkdir $@

bin/%: $(SRCFILES) t/%.c Makefile bin
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} $(SRCFILES) t/$*.c -o $@ ${LIBS}

dump:
	h5dump -d /Analyses/Basecall_2D_000/BaseCalled_template/Events $(TEST_FAST5) | less

ptdump:
	poretools events $(TEST_FAST5) | less

nanopair: bin/nanopair

test: bin/dump_fast5events
	bin/dump_fast5events $(TEST_FAST5)
