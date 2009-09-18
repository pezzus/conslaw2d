CXXCOMPILER = g++

WORKDIR = .

PROGRAM = shockbubble
CONSLAW2DDIR = ../../src/
EIGENDIR = ../../external/eigen2/
IFLAGS = -I$(CONSLAW2DDIR) -I$(EIGENDIR)

FLAGS = -O3 -msse2 -std=c++0x -pedantic -Wall

all:
	$(CXXCOMPILER) $(FLAGS) -o $(WORKDIR)/$(PROGRAM) $(PROGRAM).cpp $(IFLAGS)

clean:
	rm -f $(PROGRAM)
