CC = g++
LIBS = -lfftw3
FLAGS = -O3 -Wall

all: FNG_iso.x FNG_elastic.x

FNG_iso.x: FNG_iso.o
	$(CC) -lm FNG_iso.o -o $@ $(FLAGS) $(LIBS)

FNG_iso.o: FNG_iso.cpp
	$(CC) -c FNG_iso.cpp -o $@ $(FLAGS) $(LIBS)

FNG_elastic.x: FNG_elastic.o
	$(CC) -lm FNG_elastic.o -o $@ $(FLAGS) $(LIBS)

FNG_elastic.o: FNG_elastic.cpp
	$(CC) -c FNG_elastic.cpp -o $@ $(FLAGS) $(LIBS)

clean:
	rm -rf *.o *.x
