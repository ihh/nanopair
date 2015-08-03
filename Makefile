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
TINY_FAST5 = data/tiny.fast5

TEST_FASTA = nanopair-data/U00096.2.fas
SHORT_FASTA = nanopair-data/ch101_file20_strand.fasta
TINY_FASTA = data/tiny.fasta

ifdef HDF5
CPPFLAGS += -I${HDF5}/include -W -Wall -Wno-unused-function -Wno-unused-parameter -std=c99
LDFLAGS += -L${HDF5}/lib -Wl,-rpath -Wl,${HDF5}/lib -lz
endif

LIBS += -lhdf5_hl -lhdf5 -lgsl
# CFLAGS += -O2 -g
CFLAGS += -g

SRCFILES = $(wildcard src/*.c)


bin:
	test -e $@ || mkdir $@

bin/%: $(SRCFILES) t/%.c Makefile bin
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} $(SRCFILES) t/$*.c -o $@ ${LIBS}

dump:
	h5dump -d /Analyses/Basecall_2D_000/BaseCalled_template/Events $(TEST_FAST5) | less

ptdump:
	poretools events $(TEST_FAST5) | less

nanopair: bin/nanopair submodule

submodule: $(TEST_FAST5)

${TEST_FAST5} ${TEST_FASTA}:
	git submodule init
	git submodule update

test: bin/dump_fast5events $(TEST_FAST5)
	bin/dump_fast5events $(TEST_FAST5)

copytest: bin/copy_fast5events $(TEST_FAST5)
	bin/copy_fast5events $(TEST_FAST5) copy.fast5

tinytest: $(TINY_FAST5)

tinybug: bin/nanopair
	bin/nanopair train -eventseed $(TINY_FASTA) $(TINY_FAST5)

rep-tinybug: bin/train-tinyfast5
	bin/train-tinyfast5

$(TINY_FAST5): bin/tinyfast5
	bin/tinyfast5 $@

smalltrain: bin/nanopair submodule
	bin/nanopair train -modelseed $(SHORT_FASTA) $(TEST_FAST5)

eventseed.xml:
	bin/nanopair eventseed ~/nanopore/loman/LomanLabz_PC_E.coli_MG1655_ONI_3058_1_ch3*.fast5 >$@

