TESTS = ./tests/test_default
CFLAGS = -O2 -Wall -lm
BINDIR = ./bin
SRCDIR = ./src
DATADIR = ./output_data
BINFILE = convective_flow

all:
ifeq ($(wildcard $(BINDIR)/$(BINFILE)),)
	rm -rf $(BINDIR) $(DATADIR)
	mkdir $(BINDIR) $(DATADIR)
endif
	gcc $(SRCDIR)/main.c $(SRCDIR)/utils.c -o $(BINDIR)/$(BINFILE) $(CFLAGS)

clean:
	rm -rf $(BINDIR) $(DATADIR)

check:
	$(TESTS)
  
vg: all
	valgrind ./convective_flow

.PHONY: all clean check vg
