.SUFFIXES:.f90

PROJ=cmft

PROGS=cmft xcmft
VERSION=0.1

FC=gfortran
CC=gcc

FFLAGS=
CFLAGS=

EXTRAS=zgin.dat zg.dat

cmft_SOURCES=CMFT.f90
xcmft_SOURCES=CMFT_gnuplot.f90
nspcg_SOURCES=nspcg.f90
gnuplot_SOURCES=gnuplot.c

cmft_OBJS=${cmft_SOURCES:.f90=.o}
xcmft_OBJS=${xcmft_SOURCES:.f90=.o}
nspcg_OBJS=${nspcg_SOURCES:.f90=.o}
gnuplot_OBJS=${gnuplot_SOURCES:.c=.o}

SOURCES=${cmft_SOURCES} ${xcmft_SOURCES} ${nspcg_SOURCES} ${gnuplot_SOURCES}
OBJS=${cmft_OBJS} ${xcmft_OBJS} ${nspcg_OBJS} ${gnuplot_OBJS}

cmft_OBJS+=${nspcg_OBJS}
xcmft_OBJS+=${nspcg_OBJS} ${gnuplot_OBJS}

all: ${PROGS}

nspcg.o: ${nspcg_SOURCES}
	${FC} -std=legacy ${FFLAGS} -c $<

cmft: ${cmft_OBJS}
	${FC} ${FFLAGS} -o $@ $?
	@echo "---> Created executable $@"

xcmft: ${xcmft_OBJS}
	${FC} ${FFLAGS} -o $@ $?
	@echo "---> Created executable $@"

dist:
	@mkdir ${PROJ}-${VERSION}
	@cp ${SOURCES} ${EXTRAS} ${PROJ}-${VERSION}
	@cp Makefile ${PROJ}-${VERSION}
	@tar cvfz ${PROJ}-${VERSION}.tar.gz ${PROJ}-${VERSION} 
	@rm -rf ${PROJ}-${VERSION}
	@echo "---> Created distribution ${PROJ}-${VERSION}.tar.gz"

clean:
	rm -f ${PROGS} ${OBJS} core

.f90.o:
	${FC} ${FFLAGS} -c $<

