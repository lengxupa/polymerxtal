#
# Makefile

CWD             = `pwd`
CONFIGURE_FLAGS = --bindir=$(CWD)/../../bin --datadir=$(CWD)/../../data

all: install 

install: 
	@(cd ..; pip install -e .)
