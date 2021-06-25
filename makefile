#Version 1. This is my own makefile, it will grow.
#For now: a makefile goes like this:
#target:prereq1 prereq2
#	commands
#'all' is target, executed by make because it's the first one, 
#it says that you have to take into account main.c and cash2.c 
#all the times you change something.
#All commands must stat with a tab

PROJECT = rp_system

OBJCASH = arithmetic.o basic.o color.o filter.o io.o logical.o \
margolus.o movie.o neighbors.o noise.o png.o ps.o random.o shift.o\
x11.o

OBJCASH2 = cash2.o mersenne.o

OBJMAIN = main.o parameters.o initialisation.o graphical_output.o replication.o complex_formation.o

CC=/usr/bin/gcc

#CFLAGS=-I.
#DEPS = my.h
#CCOPT = -O0 -ggdb -Wall # For debugging.
CCOPT = -O3 -Wall #Most optimised

#all: $(OBJCASH) $(OBJCASH2) $(OBJMAIN) source.tar.gz
#	$(CPP) $(OBJCASH) $(OBJCASH2) $(OBJMAIN) $(CCOPT) -o $(PROJECT) $(LIBS) $(LDIR)

LDIR = -L/usr/X11R6/lib
LIBS = -lpng -lz -lX11 -lm #-lgrace_np -lRNA
#IDIR = -I/usr/include/ViennaRNA

all: $(OBJCASH) $(OBJCASH2) $(OBJMAIN) 
	$(CC) $(OBJCASH) $(OBJCASH2) $(OBJMAIN) $(CCOPT) -o $(PROJECT) $(LIBS) $(LDIR) 
	
main.o: my.h cash.h cash2.h mersenne.h constants.h parameters.h makefile initialisation.h replication.h complex_formation.h

$(OBJCASH): cash.h makefile
$(OBJCASH2): cash2.h mersenne.h makefile
parameters.o: parameters.h constants.h makefile
intialisation.o: initialisation.h cash.h cash2.h parameters.h constants.h makefile
graphical_output.o: graphical_output.h cash2.h cash.h parameters.h makefile
replication.o: replication.h cash2.h makefile
complex_formation.o: complex_formation.h cash2.h parameters.h constants.h

source.tar.gz: $(OBJCASH) $(OBJCASH2) $(OBJMAIN)
	tar -zcf source.tar.gz *.c *.h makefile

.c.o:
	$(CC) $(CCOPT) -c $< -o $@

#.c.o: 
#	$(CC) -c $(CCOPT) $(IDIR) $<
	
#you have to make things like arithmetic.o and I don't know if the thingy up here does it!

