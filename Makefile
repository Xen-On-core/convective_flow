TESTS = ./tests/test_default
CFLAGS = -Wall
BINDIR = ./bin
DATADIR = ./output_data
BINFILE = convective_flow

all:
ifeq ($(wildcard $(BINDIR)/$(BINFILE)),)
	rm -rf $(BINDIR) $(DATADIR)
	mkdir $(BINDIR) $(DATADIR)
endif
	gcc $(CFLAGS) -o $(BINDIR)/$(BINFILE) main.c -lm

clean:
	rm -rf $(BINDIR) $(DATADIR)

check:
	$(TESTS)
  
vg: all
	valgrind ./convective_flow

.PHONY: all clean check vg
