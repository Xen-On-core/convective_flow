TESTS = ./tests/test_default
CFLAGS = -O2 -Wall
OBJECTS = ./src/main.o \
		  ./src/utils/array.o	
INCLUDES = ./src/include/
DATADIR = ./output_data
BINFILE = ./bin/convective_flow
BINDIR = ./bin

all: utils
ifeq ($(wildcard $(BINFILE)),)
	rm -rf $(BINDIR) $(DATADIR)
	mkdir $(BINDIR) $(DATADIR)
endif
	gcc $(CFLAGS) -I $(INCLUDES) -c ./src/main.c -o ./src/main.o -lm
	gcc $(CFLAGS) -I $(INCLUDES) $(OBJECTS) -o $(BINFILE) -lm

utils: 
	gcc $(CFLAGS) -c src/utils/array.c -I $(INCLUDES) -o src/utils/array.o

clean:
	find -name '*.o' -delete
	rm -rf $(BINDIR) $(DATADIR)

check:
	$(TESTS)
  
vg: all
	valgrind $(BINFILE) --time-points 10000 --x-points 20 --y-points 20 --x_1 1.0 --y_1 1.0 --t_1 1.0 --reynolds 1.0 --grashof 10000.0 --prandtl 1.0

.PHONY: all clean check vg
