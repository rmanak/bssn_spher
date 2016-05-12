.IGNORE:
SHELL = /bin/sh
LIBS = -lm -lbbhutil -lvutil -lifcore -lbssn -lutils
LDFLAGS = -O3 -L/usr/loca/lib -L. -L./lib
CC = icc
CFLAGS = -O3
CPPFLAGS = -I. -I/usr/local/include -I./FD -I./human
CC_COMP = $(CC) -c $(CFLAGS) $(CPPFLAGS)
CC_LOAD = $(CC) $(LDFLAGS)
F77 = ifort
F77FLAGS = -fast -w -w90 -w95 -cm -Vaxlib
F77_COMP = $(F77) -c $(F77FLAGS)
F77_LOAD = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS)

EXECUTABLE = main 

.f.o:
	$(F77_COMP) $*.f

.c.o:
	$(CC_COMP) $*.c

default: all

all: FDlibs humanlibs main

FDlibs:
	cd FD; make

humanlibs:
	cd human; make

OBJECTS = main.o sdf_read.o

main: $(OBJECTS)
	$(CC_LOAD) $(OBJECTS) $(LIBS) -o main
	/bin/cp main ../bin/
	/bin/cp main ../runs/

clean:
	/bin/rm -rf $(EXECUTABLE) $(OBJECTS)
