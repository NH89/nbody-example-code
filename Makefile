
SYSTYPE = "GNU"
ZOLTANDIR = ../zoltan
SPRNGDIR = ../sprng

ifeq ($(SYSTYPE),"GNU")
CC = mpicc
OPT = -O3 -ffast-math -ftree-vectorize -fopenmp 
INCLUDES = -I. -I$(ZOLTANDIR)/include -I$(SPRNGDIR)/include
LIBS = -lm -L$(ZOLTANDIR)/lib -lzoltan -lsprng -lgmp
endif

EXEC = particles

OPTIONS = $(OPT)

OBJS = main.o init.o decompose.o io.o dynamics.o

INCL = particles.h

CFLAGS = $(OPTIONS) $(INCLUDES) 

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC) results/*


