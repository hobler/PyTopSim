# ----------------------------------------------------------------------
# Makefile for the generation of manual.pdf
#
# manual.pdf and auxiliary files are generated in the work directory
#
# "make" |
# "make all"            ... make manual.pdf
#
# "make clean"          ... delete the work manual directory 
#
# 29-Jan-2016  G.Hobler  creation
# 14-Oct-2016  G.Hobler  adjusted for topsim
# ----------------------------------------------------------------------


# Macros ---------------------------------------------------------------

#BRANCH=$(shell git describe --tags --exact-match || git symbolic-ref --short HEAD)
#BRANCH=$(shell if test "$(shell git status -z)"; then git symbolic-ref --short HEAD; else git describe --tags --exact-match || git symbolic-ref --short HEAD; fi)
BRANCH=$(shell git symbolic-ref --short HEAD)
OUTDIR=/home/hobler/work/topsim/${BRANCH}/manual
PDF=${OUTDIR}/manual.pdf

# Targets --------------------------------------------------------------

all : ${PDF}

${PDF} : algo.tex beam.tex diffusion.tex examples.tex files.tex grid.tex imsil.sty include.tex initial.tex input.tex internal.tex ion.tex mode.tex numerics.tex output.tex overview.tex parameters.tex physics.tex program.tex promis12.sty region.tex regions.tex scan.tex setup.tex time.tex times.sty titlepag.tex versions.tex makefile
	if test ! -d ${OUTDIR}; then mkdir -p ${OUTDIR}; fi
	pdflatex -output-directory=${OUTDIR} manual.tex; pdflatex -output-directory=${OUTDIR} manual.tex
	
clean :
	-rm -rf ${OUTDIR}
