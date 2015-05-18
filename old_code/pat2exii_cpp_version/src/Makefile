include ../Makefile.inc

OBJS = Main.o FortranIO.o Translator.o Helpers.o ModelClasses.o

all: $(EXE)

$(EXE): $(OBJS)
	$(LD) $(OBJS) $(LIBS) -o $@
	mv $(EXE) ..

%.o : %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

clean: 
	rm -f *.o
