DOXYGEN = doxygen
DOCSRC = $(wildcard doc/*)


all:
	@echo " === Compilazione test ==="
	$(MAKE) -C test/testmesh
	$(MAKE) -C test/sodproblem
	$(MAKE) -C test/shockbubble
	$(MAKE) -C test/shockreflection
	$(MAKE) -C test/dambreak2d
	$(MAKE) -C test/acousticwave

doc: Doxyfile $(DOCSRC)
	@echo " === Compilazione documentazione ==="
	$(DOXYGEN) Doxyfile
	$(MAKE) -C doc/report