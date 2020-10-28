#
# Makefile

CWD             = `pwd`
CONFIGURE_FLAGS = --bindir=$(CWD)/../../bin --datadir=$(CWD)/../../data

all: install 

install: 
	@(cd ..; pip install -e .)

clean:
	@(cd latch; make clean)
	@(cd PolymerBuilder; make clean)

distclean:
	@(cd latch; make distclean)
	@(cd PolymerBuilder; make distclean)

