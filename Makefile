# executable name
PROG = /home/cgiocoli/bin/MapSim.CoDECS-2.0_swap

MAIN = main.cpp readSUBFIND.cpp

LDFLAGS +=-Wl,-rpath  -lstdc++ -DSWAP #-DCOSMOLIB 

NR_DIR = /home/cgiocoli/lib/NR
COSMOLIB_DIR = /home/cgiocoli/lib/CosmoLib


# gsl, cfitsio, CCfits, fftw
LIBS = -L/home/cgiocoli/lib/gsl-1.13/lib/  -lgsl -lgslcblas  \
       -L/home/cgiocoli/lib/cfitsio/lib/ \
       -L/home/cgiocoli/lib/CCfits/lib/ -lCCfits -lcfitsio \
       -L/home/cgiocoli/lib/fftw-3.2.2/lib/ -lfftw3 -lm \
       -L$(COSMOLIB_DIR)/ -lCosmoLib \
       -L$(NR_DIR)/ -lNR \

# gsl, cfitsio, CCfits, fftw  
ALLFLAGS = -I/home/cgiocoli/lib/gsl-1.13/include/gsl/ \
	   -I/home/cgiocoli/lib/gsl-1.13/include/ \
           -I/home/cgiocoli/lib/cfitsio/include/ \
           -I/home/cgiocoli/lib/CCfits/include/ \
           -I/home/cgiocoli/lib/fftw-3.2.2/include/ \
           -I./ \
           -I$(NR_DIR)/include \
           -I$(COSMOLIB_DIR)/include \

# 
DEBUG = -g -O2 
CC = g++ 
#
RM = rm -f -r
#
OBJ = $(SOURCES:.cpp=.o)
#

CFLAGS=-O2  -g -fPIC 

default: main
main: 
	$(CC) $(CFLAGS) ${ALLFLAGS} $(MAIN) ${LIBS} ${LDFLAGS} -o ${PROG}
clean:
	$(RM) $(PROG) $(OBJ) *~

