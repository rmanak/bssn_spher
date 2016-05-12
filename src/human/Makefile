.IGNORE:

AR = xiar
SHELL = /bin/sh
LIBS = -lm -lbbhutil -lvutil -lifcore 
LDFLAGS = -O3 -L/usr/loca/lib -L.
CC = icc
CFLAGS = -O3
CPPFLAGS = -I. -I/usr/local/include
CC_COMP = $(CC) -c $(CFLAGS) $(CPPFLAGS)
CC_LOAD = $(CC) $(LDFLAGS)
F77 = ifort
F77FLAGS = -fast -w -w90 -w95 -cm -Vaxlib
F77_COMP = $(F77) -c $(F77FLAGS)
F77_LOAD = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS)


.f.o:
	$(F77_COMP) $*.f
.c.o:
	$(CC_COMP) $*.c

OBJECTS =  Divcalc.o	 compute_b.o  compute_mass.o  interpolator.o  utils.o compute_athth.o  compute_dmdr.o  compute_tt.o	 rel_diff.o

all: libutils.a install

libutils.a: $(OBJECTS)
	$(AR) r libutils.a $(OBJECTS)

install: libutils.a
	/bin/cp libutils.a ../lib/

clean:
	/bin/rm -rf $(OBJECTS) libutils.a
