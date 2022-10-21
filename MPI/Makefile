CC	= mpicc
CFLAGS	= -O0 -g -Wall
LDFLAGS = -lm -lpthread -lX11 -fopenmp
VERBOSE	=
TARGET	= nbody_brute_force nbody_barnes_hut nbody_brute_force_async nbody_barnes_hut_v2
OBJS	= ui.o xstuff.o nbody_tools.o nbody_alloc.o

#DISPLAY = -DDISPLAY
DISPLAY =

#DUMP = -DDUMP_RESULT
DUMP =

all: $(TARGET)

nbody_brute_force: nbody_brute_force.o $(OBJS)
	$(CC) $(VERBOSE) -o $@ $< $(OBJS) $(LDFLAGS)

nbody_brute_force_async: nbody_brute_force_async.o $(OBJS)
	$(CC) $(VERBOSE) -o $@ $< $(OBJS) $(LDFLAGS)

nbody_barnes_hut: nbody_barnes_hut.o $(OBJS)
	$(CC) $(VERBOSE) -o $@ $< $(OBJS)  $(LDFLAGS)

nbody_barnes_hut_v2: nbody_barnes_hut_v2.o $(OBJS)
	$(CC) $(VERBOSE) -o $@ $< $(OBJS)  $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< $(VERBOSE) $(DISPLAY) $(DUMP)
clean:
	rm -f *.o $(TARGET)
