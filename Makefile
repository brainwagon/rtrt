LIBS= -lSDL2 -lm
CFLAGS=-g -O -D_REENTRANT -I/usr/include/SDL2
CC=cc

rtrt:	rtrt.o
	$(CC) -o rtrt $(CFLAGS) rtrt.o $(LIBS)
