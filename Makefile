PROGNAME = assign3_1
SRCFILES = assign3_1.c file.c timer.c simulate.c
TARNAME = assign3_1.tgz

RUNARGS = 100 1000 # i_max t_max, increase this when testing on the DAS4!
NODES = 8 # How many DAS4 nodes
PROCESSESPERNODE = 8 # How many processes to spawn on one machine.

# if macos then
ifeq ($(shell uname), Darwin)
	PRUNARGS = -v -np $(NODES) -$(PROCESSESPERNODE) \
			   -sge-script $$PRUN_ETC/prun-openmpi
	IMAGEVIEW = open -a Preview
	CC = mpicc

	WARNFLAGS = -Wall -Werror-implicit-function-declaration -Wshadow \
			  -Wstrict-prototypes -pedantic-errors
	CFLAGS = -std=c99 -ggdb -O2 $(WARNFLAGS) -D_POSIX_C_SOURCE=200112
	LFLAGS = -lm
else
	PRUNARGS = -v -np $(NODES) -$(PROCESSESPERNODE) \
			   -sge-script $$PRUN_ETC/prun-openmpi
	IMAGEVIEW = display
	CC = mpicc

	WARNFLAGS = -Wall -Werror-implicit-function-declaration -Wshadow \
			  -Wstrict-prototypes -pedantic-errors
	CFLAGS = -std=c99 -ggdb -O2 $(WARNFLAGS) -D_POSIX_C_SOURCE=200112
	LFLAGS = -lm -lrt
endif

# Do some substitution to get a list of .o files from the given .c files.
OBJFILES = $(patsubst %.c,%.o,$(SRCFILES))

.PHONY: all run runlocal plot clean dist todo

all: $(PROGNAME)

$(PROGNAME): $(OBJFILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

run: $(PROGNAME)
	prun $(PRUNARGS) $(PROGNAME) $(RUNARGS)

runlocal: $(PROGNAME)
	mpirun -n $(PROCESSESPERNODE) $(PROGNAME) $(RUNARGS)

plot: result.txt
	gnuplot plot.gnp
	$(IMAGEVIEW) plot.png

dist:
	tar cvzf $(TARNAME) Makefile *.c *.h data/

clean:
	rm -fv $(PROGNAME) $(OBJFILES) $(TARNAME) result.txt plot.png
