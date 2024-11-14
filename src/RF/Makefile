INST_DIR = $(HOME)/bin
BIN = 
#FC = /usr/bin/gfortran
#CC = /usr/bin/g++
FC = /opt/homebrew/bin/gfortran
CC = /usr/bin/g++
FFLAGS = -O2 -Wall -ffixed-line-length-none
CFLAGS = -O2 -Wall
LDLIBS =

FOBJS = four1.o qlayer.o

theo.o: theo.f $(FOBJS)
	$(FC) $(FFLAGS) -c -o theo.o theo.f $(FOBJS) -lgfortran -L/opt/homebrew/Cellar/gcc/13.1.0/lib/gcc/13

clean : 
	rm -f $(BIN) core $(FOBJS)

