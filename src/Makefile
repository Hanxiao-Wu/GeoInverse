INST_DIR = $(HOME)/bin
BIN = MC_main
BIN1 = MC_Posterior
FC = /opt/homebrew/bin/gfortran
CC = /usr/bin/g++
#CC = /opt/homebrew/bin/g++-13
FFLAGS = -O2 -Wall -ffixed-line-length-none
CFLAGS = -std=c++20 -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ 
LDLIBS = -L/opt/homebrew/Cellar/boost/1.82.0_1/lib/ 
FOBJS = MC_main.o \
        ./DISP/surfa.o ./DISP/flat1.o ./DISP/init.o ./DISP/calcul.o ./DISP/fast_surf.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o
FOBJS1 = MC_Posterior.o \
	 ./DISP/surfa.o ./DISP/flat1.o ./DISP/init.o ./DISP/calcul.o ./DISP/fast_surf.o \
	 ./RF/theo.o ./RF/qlayer.o ./RF/four1.o

$(BIN) : $(FOBJS)
	$(CC) $(CFLAGS) $(FOBJS) -o $(BIN) ${LDLIBS}  -lgfortran -L/opt/homebrew/Cellar/gcc/13.2.0/lib/gcc/13

$(BIN1) : $(FOBJS1)
	    $(CC) $(CFLAGS) $(FOBJS1) -o $(BIN1) ${LDLIBS}  -lgfortran -L/opt/homebrew/Cellar/gcc/13.2.0/lib/gcc/13

MC_main.o : MC_main.cpp 
	$(CC)  -O2 -c $^ -I/opt/homebrew/Cellar/boost/1.82.0_1/include/  -std=c++20

MC_Posterior.o : MC_Posterior.cpp
	$(CC) -O2 -c $^ -I/opt/homebrew/Cellar/boost/1.82.0_1/include/  -std=c++20
clean ::
	rm -f $(BIN) core.* $(FOBJS) 
	rm -f $(BIN1) core.* $(FOBJS1)
