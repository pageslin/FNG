CC = g++
LIBS = -lfftw3
FLAGS = -O3 -Wall

all: noise_generator.x

noise_generator.x: noise_generator.o
	$(CC) -lm noise_generator.o -o $@ $(FLAGS) $(LIBS)

noise_generator.o: noise_generator.cpp
	$(CC) -c noise_generator.cpp -o $@ $(FLAGS) $(LIBS)

clean:
	rm -rf *.o *.x
