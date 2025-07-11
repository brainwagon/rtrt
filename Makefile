LIBS= -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
CFLAGS=-g -O
CC=cc

rtrt3:	rtrt3.o
	$(CC) -o rtrt3 $(CFLAGS) rtrt3.o $(LIBS)
