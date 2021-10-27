#
# GENESIS makefile
#
#
SOURCE=  main.f check.f diagno.f esource.f field.f incoherent.f \
 math.f partsim.f pushp.f loadbeam.f loadrad.f magfield.f \
 tdepend.f track.f string.f rpos.f scan.f source.f stepz.f \
 timerec.f initrun.f  input.f output.f mpi.f
# output.f  \
#
# LIB=  -lmdbmth -lSDDS1 -lmdblib -lz 
#
#OBJECTS =$(SOURCE:.f=.o)
#
#COMPILER = gfortran
COMPILER = mpif77
#
#  executable name
#
EXECUTABLE = genesis 
#EXECUTABLE = genesis_mpi 
#
# targets
#
genesis:	$(SOURCE)
	$(COMPILER) -g -w -O -Wall -o $(EXECUTABLE) $(SOURCE) $(LIB)
	@echo ' ******* end of linking ****** '

clean:
	rm -f *~
	rm -f *.o

install:
	cp ./$(EXECUTABLE) ~/bin/

single:
	cp mpi.f.single mpi.f
	cp mpif.h.single mpif.h

multi:
	cp mpi.f.multi mpi.f
	rm mpif.h










