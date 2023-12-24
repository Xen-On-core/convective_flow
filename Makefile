TESTS = ./tests/test_default
CFLAGS = -Wall
BINDIR = ./bin
DATADIR = ./output_data
BINFILE = convective_flow
OBJECTS = main.o \
		  src/utils/array.o	
INCLUDES = ./src/include/

all: utils
ifeq ($(wildcard $(BINDIR)/$(BINFILE)),)
	rm -rf $(BINDIR) $(DATADIR)
	mkdir $(BINDIR) $(DATADIR)
endif


	gcc $(CFLAGS) -I $(INCLUDES) -c main.c 
	gcc $(CFLAGS) -I $(INCLUDES) $(OBJECTS) -o $(BINDIR)/$(BINFILE)  -lm


utils: 
	gcc $(CFLAGS) -c src/utils/array.c -I $(INCLUDES) -o src/utils/array.o



clean:
	rm -rf $(BINDIR) $(DATADIR)

check:
	$(TESTS)
  
vg: all
	valgrind ./bin/convective_flow --time-points 10000 --x-points 20 --y-points 20 --x0 1.0 --y0 1.0 --t0 1.0 --reynolds 1.0 --grashof 10000.0 --prandtl 1.0

.PHONY: all clean check vg
