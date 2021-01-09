BIN   = ../bin/cross_spec_beam
CC    = icc
OBJ   = cross_spec_beam.o geo_seis.o
FLAGS = -lm -lfftw3

$(BIN) : $(OBJ)
	$(CC) -o $(BIN) $(OBJ) $(FLAGS)

$(OBJ) : geo_seis.h

clean : 
	rm -f $(OBJ)
