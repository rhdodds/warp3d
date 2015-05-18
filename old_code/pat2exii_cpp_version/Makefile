include Makefile.inc

all: 
	cd src; make all

check: all
	./$(EXE) < test/testcase

clean:
	rm -f $(EXE)
	rm -f *.exo
	cd src; make clean

install:
	mkdir -p $(INSTDIR)
	install $(EXE) $(INSTDIR)
